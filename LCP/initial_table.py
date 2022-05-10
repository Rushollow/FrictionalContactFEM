# This .py we need to form Contact Stiffness Matrix and Contact Load Vector
from FEM.scheme import LoadVector  # to create lv for csm
import numpy as np
from input_data import FRICTION_COEFFICIENT, INITIALIZE_GAP


class InitialTable:
    """
    Forming Contact Stiffness Matrix. Used for LCP (Linear Complementary Problem) for initial table.
    Have methods for creating CSM.
    """
    def __init__(self, element_null, stiffness_matrix, lv_const, u_linear_const,
                 u_linear_variable, lv_variable):
        """
        Initialize contact stiffness matrix
        :param element_null: FEM.element_null.ElementNullContainer
        :param stiffness_matrix: FEM.scheme.StiffnessMatrix
        :param lv_const: FEM.scheme.LoadVector (constant part of load)
        :param u_linear_const: result of solving SLAE for linear formulation with const load
        :param u_linear_variable: result of solving SLAE for linear formulation with variable load
        :param lv_variable: variable load vector for force increment algorithm
        """
        self.unit_u = np.array([], dtype=float)  # displacements from unit dislocation
        self.csm = np.array([], dtype=float)  # contact stiffness matrix
        self.clv_const = np.array([], dtype=float)  # contact load vector
        self.clv_variable = np.array([], dtype=float)
        self.table = np.array([], dtype=float)  # initial table for LCP
        self.rf = np.array([], dtype=float)  # react vector in initial table
        self.eta = np.array([], dtype=float)  # vector of gaps

        element_null.order_elements_null()  # order the elements with only correct order
        self.element_null = element_null
        self.stiffness_matrix = stiffness_matrix
        self.lv_const = lv_const
        self.u_linear_const = u_linear_const
        self.u_linear_variable = u_linear_variable
        self.lv_variable = lv_variable
        self.force_inc = False
        self.force_int_step = 0  # step number in force increment algorithm
        self.rf_const = np.zeros(1)  # to know in Lemke's algorithm if lv_const exists
        # check if load vector was set for using force incrementation algorithm
        if np.any(lv_variable.rf):
            self.force_inc = True
        self.n_amount, self.t_amount = 0, 0  # amount of normal and tangent unilateral connections in scheme
        for null_el in self.element_null:
            if null_el.orientation == 'n':
                self.n_amount += 1  # we need to remember how much normal
            else:
                self.t_amount += 1  # and how much tangential joints in contact pairs

    def __str__(self):
        """
        What we are returning if the access to our object of the class was as to string
        :return: info
        """
        return 'This is initial table for Lemke''s algorithm'

    def clear(self):
        """
        clear the memory from CSM, CLV, Initial Table itself
        :return: None
        """
        self.unit_u = np.array([], dtype=float)
        self.csm = np.array([], dtype=float)
        self.clv_const = np.array([], dtype=float)
        self.clv_variable = np.array([], dtype=float)
        self.table = np.array([], dtype=float)
        self.rf = np.array([], dtype=float)

    def _form_unit_u(self):
        """
        Forming global displacements matrix from unit dislocation of the nodes which are connected by null-elements
        :return: None
        """
        load_vector_for_csm = LoadVector(element_null=self.element_null)
        k = 0
        for element in self.element_null:
            Qe = element.form_Qe()
            load_vector_for_csm._add_Qe_to_RF(element, Qe, for_csm=True, el_null_num=k)
            k += 1
        self.unit_u = np.linalg.solve(self.stiffness_matrix.r, load_vector_for_csm.rf)

    def _form_csm(self):
        """
        Forming Contact Stiffness Matrix
        :return: None
        """
        csm = np.zeros([len(self.element_null), len(self.element_null)], dtype=float)
        # iterate over unit dislocation result of each null-element
        i, j = 0, 0
        for u_col_unit in self.unit_u.T:  # iterate over columns of array unit_u
            # iterate over each null-element
            j = 0
            for el_null in self.element_null:
                csm[j, i] = el_null.get_strain_effort(u_col_unit)
                j += 1
            i += 1
        self.csm = csm

    def _form_clv(self):
        """
        Contact Load Vector. Used for LCP (Linear Complementary Problem) initial table.
        :return: None
        """
        self.clv_const = np.zeros(len(self.element_null), dtype=float)
        if self.u_linear_const is not None:  # if we have constant load
            k = 0
            for element in self.element_null:
                self.clv_const[k] = element.get_strain_effort(self.u_linear_const)
                k += 1
        # if we have global displacements for variable load and using force increment algorithm
        if self.force_inc:
            self.clv_variable = []  # create list
            for i, u_linear_variable in enumerate(self.u_linear_variable):  # iterate over all loads
                self.clv_variable.append(np.zeros(len(self.element_null), dtype=float))  # create zeros vec in list
                for j, element in enumerate(self.element_null):  # add forces to the clv
                    self.clv_variable[i][j] = element.get_strain_effort(u_linear_variable)

    def _form_eta(self, n_amount):
        """
        Forming a vector of gaps in contact pairs which have null-elements
        :param n_amount: amount of normal unilateral contact pairs in scheme
        :return: None
        """
        self.eta = np.zeros([n_amount, 1], dtype=float)
        for i in range(n_amount):
            self.eta[i] = self.element_null.elements_list[i].gap_length

    def form_initial_table(self):
        """
        Forming an initial table for LCP (linear complementary problem)
        General view of the table in matrices:
        e|r|p|rf : where e - identity matrix, r - modified CSM, p - matrix with ones (1), rf - modified CLV
        p - coefficients of tightening weight in each contact pair on directions where null_elements were added
        f - coefficient of friction
        How its forming:
        ------------------------------------------------------------------------------------------
        |    | -r_nn                 -r_nt               r_nt            |    |  rf_n             |
        | e  | -r_tn - f * r_nn      -r_tt - f * r_nt    r_tt + f * r_nt | -p |  rf_t + f * rf_n  |
        |    |  r_tn - f * r_nn       r_tt - f * r_nt   -r_tt + f * r_nt |    |  -rf_t + f * rf_n |
        ------------------------------------------------------------------------------------------
        If we take into account gaps:
        ----------------------------------------------------------------------------------------------------------------------
        |    | -r_nn                 -r_nt               r_nt            |    |  rf_n - r_nn * eta                            |
        | e  | -r_tn - f * r_nn      -r_tt - f * r_nt    r_tt + f * r_nt | -p |  rf_t + f * (rf_n - r_nn * eta) - r_tn * eta  |
        |    |  r_tn - f * r_nn       r_tt - f * r_nt   -r_tt + f * r_nt |    |  -rf_t + f * (rf_n - r_nn * eta) + r_tn * eta |
        ------------------------------------------------------------------------------------------ ---------------------------
                If we are using force incrementation method:
        ----- -------------------csm---------------------------------- ---- -clv variables(multiple)--------------------clv constant---------------
        |    | -r_nn              -r_nt               r_nt            |    | rf_n_v                | rf_n - r_nn * eta                            |
        | e  | -r_tn - f * r_nn   -r_tt - f * r_nt    r_tt + f * r_nt | -p | rf_t_v + f * rf_n_v   | rf_t + f * (rf_n - r_nn * eta) - r_tn * eta  |
        |    |  r_tn - f * r_nn    r_tt - f * r_nt   -r_tt + f * r_nt |    | rf_t_v - f * rf_n_v   | -rf_t + f * (rf_n - r_nn * eta) + r_tn * eta |
        --------------------------------------------------------------------------------------------------------------------------------------------
                If we are using force incrementation method and there is no constant load or gaps (no constant clv):
        ----- -------------------csm------------------------------------- ---- clv variables(multiple)- -ones-
        |    | -r_nn                 -r_nt               r_nt            |    | rf_n_v                 |     |
        | e  | -r_tn - f * r_nn      -r_tt - f * r_nt    r_tt + f * r_nt | -p | rf_t_v + f * rf_n_v    |  p  |
        |    |  r_tn - f * r_nn       r_tt - f * r_nt   -r_tt + f * r_nt |    | rf_t_v - f * rf_n_v    |     |
        -----------------------------------------------------------------------------------------------

        where 'eta' - vector of gaps in contact pairs connected by null-elements
        eta_i = gap_length in i contact pair
        (len(eta) == number of contact pairs connected 'normal' to the contact zone by null elements)
        :return: None
        """
        n_amount = self.n_amount
        t_amount = self.t_amount
        self._form_unit_u()  # forming data for creating table for Lemke's algorithm
        self._form_csm()  # form contact stiffness matrix
        self._form_clv()  # form contact load vector
        self._form_eta(n_amount)  # form gaps vector
        f = FRICTION_COEFFICIENT
        r_nn = self.csm[:n_amount, :n_amount]  # CSM in normal (n) directions from unit normal (n) dislocations
        r_nt = self.csm[:n_amount, n_amount:]  # CSM in normal (n) directions from unit tangential (t) dislocations
        r_tn = self.csm[n_amount:, :n_amount]  # CSM in tangential (t) directions from unit normal (n) dislocations
        r_tt = self.csm[n_amount:, n_amount:]  # CSM in tangential (t) directions from unit tangential (t) dislocations
        rf_n = np.array([self.clv_const[:n_amount]]).T  # CLV in normal (n) directions
        rf_t = np.array([self.clv_const[n_amount:]]).T  # CLV in tangential (t) directions
        if self.force_inc:  # if using force incrementation algorithm
            rf_n_v, rf_t_v = [], []
            for i, vec in enumerate(self.clv_variable):
                rf_n_v.append(vec[:n_amount])  # CLV in normal (n) directions for variable load
                rf_t_v.append(vec[n_amount:])  # CLV in tangential (t) directions for variable load
            if t_amount == 0:  # if there is no tangent null-elements
                self._concatenate_table_n_force_inc(r_nn, rf_n, rf_n_v)
            else:
                self._concatenate_table_nt_force_inc(f, r_nn, r_nt, r_tn, r_tt, rf_n, rf_t, rf_n_v, rf_t_v)
        else:
            if t_amount == 0:
                self._concatenate_table_n(r_nn, rf_n)
            else:
                self._concatenate_table_nt(f, r_nn, r_nt, r_tn, r_tt, rf_n, rf_t)

    def _concatenate_table_n(self, *args):
        """
        Make table for unilateral frictionless contact with force incrementation algorithm
        :param args: params from self.form_initial_table() function
        :return:
        """
        r_nn, rf_n = args
        n_amount = self.n_amount
        r = r_nn  # get modified CSM
        self.rf = rf_n
        e = np.identity(n_amount, dtype=float)
        p = np.ones(shape=(n_amount, 1), dtype=float)
        self.table = np.concatenate((e, r, -p, self.rf), axis=1)

    def _concatenate_table_nt(self, *args):
        """
        Make table for unilateral frictional contact
        :param args: params from self.form_initial_table() function
        :return:
        """
        f, r_nn, r_nt, r_tn, r_tt, rf_n, rf_t = args
        n_amount = self.n_amount
        t_amount = self.t_amount
        r_row1 = np.concatenate((-r_nn, -r_nt, r_nt), axis=1)
        r_row2 = np.concatenate((np.subtract(-r_tn, f * r_nn[:t_amount, :n_amount]),
                                 np.subtract(-r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        r_row3 = np.concatenate((np.subtract(r_tn, f * r_nn[:t_amount, :n_amount]),
                                 np.subtract(r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(-r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        r = np.concatenate((r_row1, r_row2, r_row3), axis=0)  # get modified CSM
        rf_row1 = np.subtract(rf_n, r_nn.dot(self.eta))
        rf_row2 = np.add(rf_t,
                         np.subtract(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                     r_tn.dot(self.eta)))
        rf_row3 = np.add(-rf_t,
                         np.add(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                r_tn.dot(self.eta)))
        self.rf = np.concatenate((rf_row1, rf_row2, rf_row3))
        e = np.identity(n_amount + t_amount * 2, dtype=float)
        p = np.ones(shape=(n_amount + t_amount * 2, 1), dtype=float)
        self.table = np.concatenate((e, r, -p, self.rf), axis=1)

    def _concatenate_table_n_force_inc(self, *args):
        """
        Make table for unilateral frictionless contact with force incrementation algorithm
        :param args: params from self.form_initial_table() function
        :return:
        """
        r_nn, rf_n, rf_n_v = args
        n_amount = self.n_amount

        r = r_nn  # get modified CSM
        # forming constant react vector for table
        self.rf_const = np.subtract(rf_n, r_nn.dot(self.eta))
        # forming variable react vector for table
        rf_n_v = np.array(rf_n_v).T
        rf_variable = rf_n_v
        e = np.identity(n_amount, dtype=float)
        p = np.ones(shape=(n_amount, 1), dtype=float)
        if not np.any(self.rf_const):  # if there is no constant load (gaps or forces) than change table:
            self.table = np.concatenate((e, r, -p, -rf_variable, p), axis=1)
        else:  # general case
            self.table = np.concatenate((e, r, -p, -rf_variable, self.rf_const), axis=1)

    def _concatenate_table_nt_force_inc(self, *args):
        """
        Make table for unilateral frictional contact with force incrementation algorithm
        :param args: params from self.form_initial_table() function
        :return:
        """
        f, r_nn, r_nt, r_tn, r_tt, rf_n, rf_t, rf_n_v, rf_t_v = args
        n_amount = self.n_amount
        t_amount = self.t_amount

        r_row1 = np.concatenate((-r_nn, -r_nt, r_nt), axis=1)
        r_row2 = np.concatenate((np.subtract(-r_tn, f * r_nn[:t_amount, :n_amount]),
                                 np.subtract(-r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        r_row3 = np.concatenate((np.subtract(r_tn, f * r_nn[:t_amount, :n_amount]),
                                 np.subtract(r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(-r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        r = np.concatenate((r_row1, r_row2, r_row3), axis=0)  # get modified CSM
        # forming constant react vector for table
        rf_row1 = np.subtract(rf_n, r_nn.dot(self.eta))
        rf_row2 = np.add(rf_t,
                         np.subtract(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                     r_tn.dot(self.eta)))
        rf_row3 = np.add(-rf_t,
                         np.add(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                r_tn.dot(self.eta)))
        self.rf_const = np.concatenate((rf_row1, rf_row2, rf_row3))
        # forming variable react vector for table
        # rf_var_row1 = rf_n_v  # np.expand_dims(arr2, 1)
        # rf_var_row2 = np.add(rf_t_v, f * rf_n_v[:t_amount])
        # rf_var_row3 = np.add(-rf_t_v, f * rf_n_v[:t_amount])
        rf_n_v = np.array(rf_n_v).T
        rf_t_v = np.array(rf_t_v).T
        rf_var_row1 = rf_n_v  # np.expand_dims(arr2, 1)
        rf_var_row2 = np.add(f * rf_n_v[:t_amount], rf_t_v)
        rf_var_row3 = np.subtract(f * rf_n_v[:t_amount], rf_t_v)
        rf_variable = np.concatenate((rf_var_row1, rf_var_row2, rf_var_row3))
        e = np.identity(n_amount + t_amount * 2, dtype=float)
        p = np.ones(shape=(n_amount + t_amount * 2, 1), dtype=float)
        if not np.any(self.rf_const):  # if there is no constant load (gaps or forces) than change table:
            self.table = np.concatenate((e, r, -p, -rf_variable, p), axis=1)
        else:  # general case
            self.table = np.concatenate((e, r, -p, -rf_variable, self.rf_const), axis=1)








