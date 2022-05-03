import numpy as np
from input_data import LEMKE_LIMIT_STEPS, ACCURACY_OF_LCP
from input_data import INITIALIZE_GAP


class Lemke:
    """
    Class for Lemke functions and solver of LCP
    """
    def __init__(self, intl_table):
        """
        Constructor
        :param intl_table: Initial table formed for Lemke's algorithm
        if False: usual Lemke's algorithm will be held
        if True: force incrementation algorithm will be held
        """
        self._rows_table = int((intl_table.table.shape[0]))
        self.force_inc = intl_table.force_inc
        # step number of the force increment algorithm, if = 0, than solving normal Lemke
        self.force_inc_step = 0
        # used to different logic of getting the results xn,xt,zn,zt from table from last step
        self.last_force_inc_step = False
        # if force incrementation algorithm was chosen than change initial table:
        self.react_vector = intl_table.table[:, -1]  # solving normal Lemke
        self.table = intl_table.table
        self.n_amount = intl_table.n_amount
        self.t_amount = intl_table.t_amount
        # we need pairs of "variable" and "variable_next" to know values on previous steps
        self.basis = np.arange(self._rows_table)
        self._basis_next = self.basis.copy()
        self._leading_column = self._rows_table * 2 + self.force_inc_step
        self._leading_column_next = self._leading_column
        self._min_ratio = np.zeros(self._rows_table, dtype=float)
        if self.force_inc:
            self._form_min_ratio()
            self._form_leading_row()
        else:
            self._leading_row_next = np.argmin(self.react_vector)
        self._leading_row = self._leading_row_next
        # Result interaction forces and mutual displacements
        self.xn = np.zeros(self.n_amount, dtype=float)  # interaction forces along the normal to the contact zone
        self.xt = np.zeros(self.t_amount, dtype=float)  # interaction forces tangential to the contact zone
        self.zn = np.zeros(self.n_amount, dtype=float)  # mutual displacements along the normal to the contact zone
        self.zt = np.zeros(self.t_amount, dtype=float)  # mutual displacements tangential to the contact zone
        self.p_value = 0
        # Results as animation
        self.xn_anim = []
        self.xt_anim = []
        self.zn_anim = []
        self.zt_anim = []
        self.p_anim = []  # value of the parameter step by step
        # number of steps in algorithm (how many steps were held in total)
        # start from -1 case at the beginning we add the data and make += 1 step, so for program it will be step num. 0
        self.steps = -1

    def clear(self):
        """
        Clear all data in class
        :return:
        """
        self._rows_table = None
        self.table = None
        self.react_vector = None
        self.basis = None
        self._basis_next = None
        self._leading_column = None
        self._leading_column_next = None
        self._leading_row = None
        self._leading_row_next = None

    # region Different checks
    def _trivial_solution(self):
        """
        Checks if it is trivial solution. If so, we do not need to solve LCP
        :return: True or False. Bool value
        """
        if np.min(self.react_vector) > 0 and not self.force_inc:
            return True
        else:
            return False

    def _ray_solution(self):
        """
        Checks if we have ray solution
        :return: bool value
        """
        arr = np.where(self._min_ratio >= 0, self._min_ratio, np.inf)
        if np.min(arr) == np.Infinity:
            return True
        else:
            return False

    def _rough_solution(self):
        """
        Checks if we got rough solution solving LCP. We can get rough solution only if ray_solution occurred
        :return: bool value
        """
        if self._ray_solution():
            # p - tightening weight in each contact pair on directions where null_elements were added
            p_index = np.where(self.basis == self._rows_table * 2)  # get the position in basis where is p value
            if self.force_inc:
                p_value = self.table[p_index, -2]  # get the p value
            else:
                p_value = self.table[p_index, -1]  # get the p value
            if p_value < ACCURACY_OF_LCP:
                return True
        else:
            return False

    def _p_in_basis(self):
        """
        Check if p (tightening weight) in basis
        :return: bool value
        """
        p_index = np.where(self.basis == self._rows_table * 2)[0]
        return False if p_index.size == 0 else True
    # endregion

    def _get_p(self):
        """
        Get p (tightening weight) value
        :return: p
        """
        p_row_index = np.where(self.basis == self._rows_table * 2)[0]
        if not p_row_index:
            return 0
        return self.table[p_row_index[0], self._rows_table * 2]

    # region Lemke step
    def _form_table(self):
        """
        Forming new table for next Lemke step
        :return: None
        """
        leading_element = self.table[self._leading_row, self._leading_column]
        table_tmp = np.zeros(self.table.shape, dtype=float)
        for i, j in np.ndindex(self.table.shape):
            if i != self._leading_row:
                table_tmp[i, j] = self.table[i, j] - \
                                  self.table[i, self._leading_column] * self.table[self._leading_row, j] / leading_element
            else:
                table_tmp[i, j] = self.table[i, j] / leading_element
        self.table = table_tmp

    def _form_basis(self):
        """
        Forming next basis for next Lemke step
        :return: None
        """
        self._basis_next[self._leading_row] = self._leading_column

    def _form_leading_column(self):
        """
        Finding next leading column for new Lemke step
        :return: None
        """
        tmp = self.basis[self._leading_row]
        if tmp < self._rows_table:
            self._leading_column_next = tmp + self._rows_table
        else:
            self._leading_column_next = tmp - self._rows_table

    def _form_min_ratio(self):
        for i in range(self._rows_table):
            element = self.table[i, self._leading_column_next]  # element in table in leading column
            if element != 0:
                self._min_ratio[i] = self.table[i, -1] / element
            else:
                self._min_ratio[i] = np.Infinity

    def _form_leading_row(self):
        """
        Find next leading row for next step for Lemke
        :return: None
        """
        # take all elements < 0 and make them Infinity
        min_ratio_above_0 = np.where(self._min_ratio >= 0, self._min_ratio, np.inf)
        # get indices of the sorted array
        indx = np.argsort(min_ratio_above_0)
        # check if in this raw(index) is p (tightening weight)
        for i in range(1, len(indx)):  # iterate over array from min to max values
            # if next value CLOSE to minimum
            if np.isclose(min_ratio_above_0[indx[0]], min_ratio_above_0[indx[i]], atol=ACCURACY_OF_LCP):
                # check if p in this row
                if self.basis[indx[i]] == self._rows_table * 2:
                    # if so - choose row with p value to end lemke's algorithm
                    self._leading_row_next = indx[i]
                    break
            else:
                self._leading_row_next = indx[0]
                break
    # endregion

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
        # take react vector from table
        if self.force_inc:
            self._results_force_inc()
            self.p_value = -self._get_p()
        else:
            self.react_vector = self.table[:, -1]
            self.p_value = self._get_p()
        # form results
        self._form_xn()
        self._form_xt()
        self._form_zn()
        self._form_zt()

    def _results_force_inc(self):
        """
        Form results for force increment algorithm
        :return:
        """
        # if it's the zero step - just take existing values from table
        if self.last_force_inc_step or self.steps == -1:
            self.react_vector = self.table[:, -1] - self.table[:, self._rows_table * 2 + self.force_inc_step]
            self.last_force_inc_step = False  # to allow use another logic for forming the results
            return
        # on other steps make additional step to get the table which is needed and then get the results
        table_state = self.table.copy()  # remember the table condition
        row_state = self._leading_row  # remember leading row
        # make one step of algorithm to remove 'p' from basis and form results
        a = np.where(self.basis == self._rows_table * 2)
        print(a[0])
        p_row_index = np.where(self.basis == self._rows_table * 2)[0]
        self._leading_row = p_row_index
        self._form_table()
        self.react_vector = self.table[:, -1] - self.table[:, self._rows_table * 2 + self.force_inc_step]
        # take variables back
        self.table = table_state
        self._leading_row = row_state

    def _form_xn(self):
        # set previous results to zero
        self.xn = np.zeros(self.n_amount, dtype=float)
        # form xn (get all indices where numbers in basis < n_amount)
        for i in np.where(self.basis < self.n_amount):
            k = self.basis[i]
            self.xn[k] = self.react_vector[i]

    def _form_xt(self):
        # set previous results to zero
        self.xt = np.zeros(self.t_amount, dtype=float)
        xtp = np.zeros(self.t_amount, dtype=float)  # form xtp (xt plus)
        for i in np.where((self.basis >= self.n_amount) &
                          (self.basis < self.n_amount + self.t_amount)):
            k = self.basis[i] - self.n_amount
            xtp[k] = self.react_vector[i]
        xtm = np.zeros(self.t_amount, dtype=float)  # form xtm (xt minus)
        for i in np.where((self.basis >= self.n_amount + self.t_amount) &
                          (self.basis < self.n_amount + self.t_amount * 2)):
            k = self.basis[i] - (self.n_amount + self.t_amount)
            xtm[k] = self.react_vector[i]
        for i in range(self.t_amount):
            self.xt[i] = (xtp[i] - xtm[i]) / 2

    def _form_zn(self):
        # set previous results to zero
        self.zn = np.zeros(self.n_amount, dtype=float)
        # form zn
        for i in np.where((self.basis >= self.n_amount + self.t_amount * 2) &
                          (self.basis < self.n_amount * 2 + self.t_amount * 2)):
            k = self.basis[i] - (self.n_amount + self.t_amount * 2)
            self.zn[k] = self.react_vector[i]

    def _form_zt(self):
        # set previous results to zero
        self.zt = np.zeros(self.t_amount, dtype=float)
        # form zt
        ztp = np.zeros(self.t_amount, dtype=float)   # form ztp (zt plus)
        indices_ztp = np.where((self.basis >= self.n_amount * 2 + self.t_amount * 2) &
                               (self.basis < self.n_amount * 2 + self.t_amount * 3))[0]
        for i in indices_ztp:
            k = self.basis[i] - (self.n_amount * 2 + self.t_amount * 2)
            ztp[k] = self.react_vector[i]
        ztm = np.zeros(self.t_amount, dtype=float)
        indices_ztm = np.where((self.basis >= self.n_amount * 2 + self.t_amount * 3) &
                               (self.basis < self.n_amount * 2 + self.t_amount * 4))[0]
        for i in indices_ztm:
            k = self.basis[i] - (self.n_amount * 2 + self.t_amount * 3)
            ztm[k] = self.react_vector[i]
        for i in range(self.t_amount):
            self.zt[i] = ztp[i] - ztm[i]

    def _results_anim(self):
        """
        Append the results of the LCP solution on each(one) step
        :return: None
        """
        self._results()  # get results
        # write (remember) data in lists
        self.xn_anim.append(self.xn)
        self.xt_anim.append(self.xt)
        self.zn_anim.append(self.zn)
        self.zt_anim.append(self.zt)
        self.p_anim.append(self.p_value)
        self.steps += 1

    def _lemke_step(self):
        """
        One step of Lemke algorithm
        :return: None
        """
        self._form_table()
        self._form_basis()
        self._form_leading_column()
        self._form_min_ratio()
        self._form_leading_row()

    def lcp_solve(self):
        """
        Linear Complementary Problem solver
        :return: None
        """
        if self._trivial_solution():
            print('LCP trivial solution')
            self._results_anim()  # end
            return
        for step in range(LEMKE_LIMIT_STEPS):
            self._results_anim()
            self._lemke_step()  # do the step
            # remember data from previous steps
            self.basis = self._basis_next.copy()
            self._leading_column = self._leading_column_next
            self._leading_row = self._leading_row_next
            # Do checks
            if self._ray_solution():
                if not self._rough_solution():
                    print('Ray solution in LCP, you got the results for "almost broken" system,'
                          ' its {} step of Lemke, p={}'
                          .format(step, self._get_p()))
                    self._results_anim()
                else:
                    print('Rough solution is occurred, you got rough results '
                          'as there in each contact pair the p (tightening weight) still exists\n'
                          'p={} is lesser than value ACCURACY_OF_LCP={} settled by user'
                          .format(self._get_p(), ACCURACY_OF_LCP))
                    self._results_anim()
                # if incrementation force algorithm than:
                if self.force_inc:
                    # if number of leading column was chosen last time penultimate (предпоследняя) column
                    if self._rows_table * 2 + self.force_inc_step == self.table.shape[1]-1:
                        self.last_force_inc_step = True  # need this to form results for this next step properly
                        p_row_index = np.where(self.basis == (self._rows_table * 2 + self.force_inc_step))[0]
                        self._leading_row = p_row_index
                    self.modify_data_for_second_step()  # continue to solve force increment algorithm
                else:
                    return  # end
            elif not self._p_in_basis():
                print("Normal solution of LCP in {} steps".format(step+1))
                if self.force_inc:  # if it is force incrementation algorithm
                    if self._rows_table * 2 + self.force_inc_step == self.table.shape[1] - 1:
                        self.last_force_inc_step = True
                        self._results_anim()
                        return  # end
                    self.modify_data_for_second_step()  # continue to solve force increment algorithm
                else:
                    self._results_anim()
                    return  # end
        print(f"LCP solver got maximum number of steps, no objective results, p={self._get_p()}")
        self._results_anim()
        return  # end

    def modify_data_for_second_step(self):
        """
        Modifying table to prepare it for using as usual table for increment force algorithm
        :return:
        """
        self.force_inc_step += 1
        # next we set table like it's first step
        self._leading_column = self._rows_table * 2 + self.force_inc_step  # choose new leading column
        # choose new leading row (must choose minimum positive value), changed all negative values to +inf and takes min
        self._leading_row = np.where(self.react_vector > 0, self.react_vector, np.inf).argmin()
