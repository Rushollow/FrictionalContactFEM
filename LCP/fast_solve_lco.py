import numpy as np
from FEM.scheme import LoadVector
from input_data import FRICTION_COEFFICIENT


class LCPsolver:
    """
    solver of the Linear Complementary Problem
    """

    def __init__(self, element_null, stiffness_matrix, u_linear):
        """
        Initialize contact stiffness matrix and vector
        """
        self.unit_u = np.array([], dtype=float)
        self.u_linear = u_linear  # vector of global displacements in system with bounded contact conditions
        self.stiffness_matrix = stiffness_matrix

        element_null.order_elements_null()  # order the elements with only correct order
        self.element_null = element_null
        self.solution = np.array([], dtype=float)

        self._form_unit_u()  # displacements from unit dislocation
        self._form_csm()  # contact stiffness matrix
        self._form_clv()  # contact load vector

        self.n_amount, self.t_amount = 0, 0  # amount of normal and tangent unilateral connections in scheme
        for null_el in self.element_null:
            if null_el.orientation == 'n':
                self.n_amount += 1  # we need to remember how much normal
            else:
                self.t_amount += 1  # and how much tangential joints in contact pairs
        self.M = np.array([], dtype=float)
        self.q = np.array([], dtype=float)
        self._form_matrix_and_vector()  # form M and q: matrix and load vector for LCP

        self.xn = np.zeros(self.n_amount, dtype=float)  # interaction forces along the normal to the contact zone
        self.xt = np.zeros(self.t_amount, dtype=float)  # interaction forces tangential to the contact zone
        self.zn = np.zeros(self.n_amount, dtype=float)  # mutual displacements along the normal to the contact zone
        self.zt = np.zeros(self.t_amount, dtype=float)  # mutual displacements tangential to the contact zone

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
        self.clv = np.zeros(len(self.element_null), dtype=float)
        k = 0
        for element in self.element_null:
            self.clv[k] = element.get_strain_effort(self.u_linear)
            k += 1

    def _form_eta(self, n_amount):
        """
        Forming a vector of gaps in contact pairs which have null-elements
        :param n_amount: amount of normal unilateral contact pairs in scheme
        :return: None
        """
        self.eta = np.zeros([n_amount, 1], dtype=float)
        for i in range(n_amount):
            self.eta[i] = self.element_null.elements_list[i].gap_length

    def _form_matrix_and_vector(self):
        """
        Forming an initial table for LCP (linear complementary problem)
        General view of the table in matrices:
        e|r|p|rf : where e - identity matrix, r - modified CSM, p - matrix with ones (1), rf - modified CLV
        p - coefficients of tightening weight in each contact pair on directions where null_elements were added
        f - coefficient of friction
        How its forming:
        --------------------------------------------------------------
        | -r_nn                 -r_nt               r_nt            |
        | -r_tn - f * r_nn      -r_tt - f * r_nt    r_tt + f * r_nt |
        |  r_tn - f * r_nn       r_tt - f * r_nt   -r_tt + f * r_nt |
        --------------------------------------------------------------
        |  rf_n            |
        |  rf_t + f * rf_n |
        |  rf_t - f * rf_n |
        If we take into account gaps we change load vector:
        ----------------------------------------------
         rf_n - r_nn * eta
         rf_t + f * (rf_n - r_nn * eta) - r_tn * eta
         -rf_t + f * (rf_n - r_nn * eta) + r_tn * eta
        ------------------ ---------------------------
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
        rf_n = np.array([self.clv[:n_amount]]).T  # CLV in normal (n) directions
        rf_t = np.array([self.clv[n_amount:]]).T  # CLV in tangential (t) directions
        # r_row1 = np.concatenate((-r_nn, -r_nt, r_nt), axis=1)
        # r_row2 = np.concatenate((np.subtract(-r_tn, f * r_nn[:t_amount, :n_amount]),
        #                          np.subtract(-r_tt, f * r_nt[:t_amount, :t_amount]),
        #                          np.add(r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        # r_row3 = np.concatenate((np.subtract(r_tn, f * r_nn[:t_amount, :n_amount]),
        #                          np.subtract(r_tt, f * r_nt[:t_amount, :t_amount]),
        #                          np.add(-r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        # rf_row1 = np.subtract(rf_n, r_nn.dot(self.eta))
        # rf_row2 = np.add(rf_t,
        #                  np.subtract(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
        #                              r_tn.dot(self.eta)))
        # rf_row3 = np.add(-rf_t,
        #                  np.add(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
        #                         r_tn.dot(self.eta)))
        #
        # self.M = np.concatenate((r_row1, r_row2, r_row3), axis=0)  # get modified CSM
        # self.q = np.concatenate((rf_row1, rf_row2, rf_row3))

        r_row2 = np.concatenate((
                                 np.subtract(-r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        r_row3 = np.concatenate((
                                 np.subtract(r_tt, f * r_nt[:t_amount, :t_amount]),
                                 np.add(-r_tt, f * r_nt[:t_amount, :t_amount])), axis=1)
        rf_row2 = np.add(rf_t,
                         np.subtract(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                     r_tn.dot(self.eta)))
        rf_row3 = np.add(-rf_t,
                         np.add(f * (np.subtract(rf_n[:t_amount], r_nn.dot(self.eta)[:t_amount])),
                                r_tn.dot(self.eta)))

        self.M = np.concatenate((r_row2, r_row3), axis=0)  # get modified CSM
        self.q = np.concatenate((rf_row2, rf_row3))

    def _results(self):
        """
        Form and write (remember):
        interaction forces normal: xn
        interaction forces tangent: xt
        mutual displacements normal: zn
        mutual displacements tangent: zt
        by using react vector from Lemke algorithm table
        :return: None
        """
        # form results ????
        self._form_xn()
        self._form_xt()
        self._form_zn()
        self._form_zt()

    def matrix_minor(self, arr, i, j):
        return np.linalg.det(np.delete(np.delete(arr, i, axis=0), j, axis=1))

    def solve_lcp(self):
        """
        LCP solving algorithm
        :return:
        """
        M = self.M
        q = self.q.reshape(-1, 1)
        # M = np.array([[4,2,0,3],
        #                [-1,4,-3,-6],
        #                [1,-1,1,1],
        #                [0,1,0,5]])
        # q = np.array([-1,2,-1,-1]).reshape(-1, 1)
        # solution is [0, 1, 2, 0]
        n = q.shape[0]
        p = np.zeros(n)
        I = np.identity(n)

        z = np.sign(np.linalg.inv(I + M) @ q)  # @ - is dot product of 2 matrices
        Tz = np.diag(z.reshape(-1))
        test_m = (I + M) + (I - M) @ Tz
        print(f'eigenvalues M: {np.linalg.eigvals(M)}\n'
                f'def M: {np.linalg.det(M)}')
        for i in range(M.shape[0]):
            print(self.matrix_minor(M, i, i))
        print(f'rank: {np.linalg.matrix_rank(test_m)}\n'
              f'rows: {test_m.shape[0]}\n'
              f'cond1 {np.linalg.cond(test_m)}\n'
              f'det {np.linalg.det(test_m)}')
        H = np.linalg.inv((I + M) + (I - M) @ Tz)
        x = H @ q
        C = -H @ (I - M)

        # Проверка на наличие отрицательных чисел в векторе
        xz = x * z
        neg_xz_vals = np.where(xz < 0)[0]
        k = 0
        while neg_xz_vals.size > 0 and k < 100:
            k += 1
            j = neg_xz_vals[0]

            if 1 + (2 * z[j] * C[j, j]) <= 0:
                print('No solution! 1st check')
                self.solution = np.abs(x) - x
                return

            p[j] = p[j] + 1
            if np.log2(p[j]) > n - j:
                print('No solution! 2nd check')
                self.solution = np.abs(x) - x
                return

            if 1 + 2 * z[j] * C[j, j] > 0:
                ej = I[:, j].reshape(n, 1)
                Tz = Tz - 2 * z[j] * ej * ej.T
                z[j] = -z[j]
                alpha = (2 * z[j]) / (1 - 2 * z[j] * C[j, j])
                x = x + alpha * x[j] * C[:, j].reshape(-1, 1)
                C = C + alpha * C[:, j] @ C[j, :]

            xz = x * z
            neg_xz_vals = np.where(xz < 0)[0]
        self.solution = np.abs(x) - x
