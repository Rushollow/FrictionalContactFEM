import numpy as np
from input_data import LEMKE_LIMIT_STEPS, ACCURACY_OF_LCP, WRITE_EXCEL
import xlsxwriter



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
        self.const_load = np.any(intl_table.rf_const)  # bool, if there is const load, first solve lemke for it
        # step number of the force increment algorithm, if = 0, than solving normal Lemke
        self.force_inc_without_const = False  # if there is const load (change how results form in case of ray solution)
        self.force_inc_step = 0  # there could be multiple load vectors
        if not self.const_load:
            self.force_inc_step = 1
            self.force_inc_without_const = True  # if there is no const load (change how results form if ray solution)
        # if force incrementation algorithm was chosen than change initial table:
        self.react_vector = intl_table.table[:, -1]  # solving normal Lemke
        self.table = intl_table.table
        self.n_amount = intl_table.n_amount
        self.t_amount = intl_table.t_amount
        # we need pairs of "variable" and "variable_next" to know values on previous steps
        self._basis = np.arange(self._rows_table)
        self._basis_next = self._basis.copy()
        self._leading_column = self._rows_table * 2 + self.force_inc_step
        self._leading_column_next = self._leading_column
        self._min_ratio = np.zeros(self._rows_table, dtype=float)
        if self.force_inc and not self.const_load:  # if there is no const load and force increment
            self._leading_row_next = 0
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
        # used to form u_contact_anim in "calculation" module. form_u_contact - function
        # nested arr (2 dimentions) dict: keys are step numbers, values are tuple (which lv_var numer, p_value)
        self.p_anim_variable = {}
        # number of steps in algorithm (how many steps were held in total)
        # start from -1 case at the beginning we add the data and make += 1 step, so for program it will be step num. 0
        self.steps = -1
        # variables to roll back the step of Lemke:
        self.table_previous = np.zeros(self.table.shape, dtype=float)
        self._basis_previous = self._basis.copy()

        self.write_excel = WRITE_EXCEL
        if self.write_excel:
            print("Excel is writing all tables!")
            self.workbook = xlsxwriter.Workbook('initial_table.xlsx')
            self.worksheet = self.workbook.add_worksheet(name='Sheet1')
            self.excel_table_count = 0
            for i in range(self.n_amount):  # write range (top row)
                self.worksheet.write(0, 1 + i, 'xn'+str(i))
            for i in range(self.t_amount):
                self.worksheet.write(0, 1 + i + self.n_amount, 'xt' + str(i) + '+')
            for i in range(self.t_amount):
                self.worksheet.write(0, 1 + i + self.n_amount+self.t_amount, 'xt' + str(i) + '-')
            for i in range(self.n_amount):  # write range (top row)
                self.worksheet.write(0, 1 + i + self.n_amount+self.t_amount*2, 'zn'+str(i))
            for i in range(self.t_amount):
                self.worksheet.write(0, 1 + i + self.n_amount*2+self.t_amount*2, 'zt' + str(i) + '+')
            for i in range(self.t_amount):
                self.worksheet.write(0, 1 + i + self.n_amount*2+self.t_amount*3, 'zt' + str(i) + '-')

    def clear(self):
        """
        Clear all data in class
        :return:
        """
        self._rows_table = None
        self.table = None
        self.react_vector = None
        self._basis = None
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
        if np.min(self.react_vector) > 0:
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
            p_value = self._get_p()
            if ACCURACY_OF_LCP > p_value > 0:
                return True
        else:
            return False

    def _p_in_basis(self):
        """
        Check if p (tightening weight) in basis
        :return: bool value
        """
        p_index = np.where(self._basis == self._rows_table * 2 + self.force_inc_step)[0]
        return False if p_index.size == 0 else True
    # endregion

    def _get_p(self):
        """
        Get p (tightening weight) value
        :return: p value
        """
        p_row_index = np.where(self._basis == (self._rows_table * 2 + self.force_inc_step))[0]
        if p_row_index.size == 0:
            return 0
        return self.react_vector[p_row_index[0]]

    def _get_p_row(self):
        """
        Gets row in which 'p' is
        If p not in basis return None
        :return:
        """
        if self._p_in_basis():
            return np.where(self._basis == self._rows_table * 2 + self.force_inc_step)[0][0]
        raise Exception(f'you tried to get p row at {self.steps} step of Lemke, but its not in basis')

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
        # if there are any negative numbers in react vector - change them to zero. It's error rate (погрешность)
        if not self.force_inc:
            for i in range(table_tmp.shape[0]):
                if table_tmp[i, -1] < 0:
                    table_tmp[i, -1] = 0
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
        tmp = self._basis[self._leading_row]
        if tmp < self._rows_table:
            self._leading_column_next = tmp + self._rows_table
        else:
            self._leading_column_next = tmp - self._rows_table

    def _form_min_ratio(self):
        for i in range(self._rows_table):
            element = self.table[i, self._leading_column_next]  # element in table in leading column
            # if element > 0:
            if not np.isclose(element, 0, atol=ACCURACY_OF_LCP) and element > 0:
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
        indx: list = np.argsort(min_ratio_above_0)
        # check if in this raw(index) is p (tightening weight)
        for i in range(1, len(indx)):  # iterate over array from min to max values
            # if next value CLOSE to minimum
            if np.isclose(min_ratio_above_0[indx[0]], min_ratio_above_0[indx[i]], atol=ACCURACY_OF_LCP):
                # check if p in this row
                if self._basis_next[indx[i]] == self._rows_table * 2 + self.force_inc_step:
                    # if so - choose row with p value to end lemke's algorithm
                    self._leading_row_next = indx[i]
                    break
            else:
                self._leading_row_next = indx[0]
                break
    # endregion

    def _form_xn(self):
        # set previous results to zero
        self.xn = np.zeros(self.n_amount, dtype=float)
        # form xn (get all indices where numbers in basis < n_amount)
        for i in np.where(self._basis < self.n_amount):
            k = self._basis[i]
            self.xn[k] = self.react_vector[i]

    def _form_xt(self):
        # set previous results to zero
        self.xt = np.zeros(self.t_amount, dtype=float)
        xtp = np.zeros(self.t_amount, dtype=float)  # form xtp (xt plus)
        for i in np.where((self._basis >= self.n_amount) &
                          (self._basis < self.n_amount + self.t_amount)):
            k = self._basis[i] - self.n_amount
            xtp[k] = self.react_vector[i]
        xtm = np.zeros(self.t_amount, dtype=float)  # form xtm (xt minus)
        for i in np.where((self._basis >= self.n_amount + self.t_amount) &
                          (self._basis < self.n_amount + self.t_amount * 2)):
            k = self._basis[i] - (self.n_amount + self.t_amount)
            xtm[k] = self.react_vector[i]
        for i in range(self.t_amount):
            self.xt[i] = (xtp[i] - xtm[i]) / 2

    def _form_zn(self):
        # set previous results to zero
        self.zn = np.zeros(self.n_amount, dtype=float)
        # form zn
        for i in np.where((self._basis >= self.n_amount + self.t_amount * 2) &
                          (self._basis < self.n_amount * 2 + self.t_amount * 2)):
            k = self._basis[i] - (self.n_amount + self.t_amount * 2)
            self.zn[k] = self.react_vector[i]

    def _form_zt(self):
        # set previous results to zero
        self.zt = np.zeros(self.t_amount, dtype=float)
        # form zt
        ztp = np.zeros(self.t_amount, dtype=float)   # form ztp (zt plus)
        indices_ztp = np.where((self._basis >= self.n_amount * 2 + self.t_amount * 2) &
                               (self._basis < self.n_amount * 2 + self.t_amount * 3))[0]
        for i in indices_ztp:
            k = self._basis[i] - (self.n_amount * 2 + self.t_amount * 2)
            ztp[k] = self.react_vector[i]
        ztm = np.zeros(self.t_amount, dtype=float)
        indices_ztm = np.where((self._basis >= self.n_amount * 2 + self.t_amount * 3) &
                               (self._basis < self.n_amount * 2 + self.t_amount * 4))[0]
        for i in indices_ztm:
            k = self._basis[i] - (self.n_amount * 2 + self.t_amount * 3)
            ztm[k] = self.react_vector[i]
        for i in range(self.t_amount):
            self.zt[i] = ztp[i] - ztm[i]

    def _results(self, ray=False):
        """
        Form and write (remember):
        interaction forces normal: xn
        interaction forces tangent: xt
        mutual displacements normal: zn
        mutual displacements tangent: zt
        by using react vector from Lemke algorithm table
        :param ray: bool, if it is ray solution and force inc. alg. then use another logic for results
        :return: None
        """
        # take react vector from table
        if ray:
            self._results_force_inc()
        else:
            self.react_vector = self.table[:, -1]
            self.p_value = abs(self._get_p())
        # form results
        self._form_xn()
        self._form_xt()
        self._form_zn()
        self._form_zt()

    def _results_force_inc(self):
        """
        Form results for force increment algorithm ray solution
        :return:
        """

        self._leading_row = self._get_p_row()
        self.p_value = abs(self._get_p())
        self._lemke_step()
        if self.force_inc_without_const:
            self.react_vector = - self.table[:, self._rows_table * 2 + self.force_inc_step] * self.p_value
        else:
            self.react_vector = self.table[:, -1] - \
                                self.table[:, self._rows_table * 2 + self.force_inc_step] * self.p_value

    def _results_anim(self, ray=False):
        """
        Append the results of the LCP solution on each(one) step
        :param ray: bool, if it is ray solution and force inc. alg. then use another logic for results
        :return: None
        """
        self._results(ray=ray)  # get results
        # write (remember) data in lists
        self.xn_anim.append(self.xn)
        self.xt_anim.append(self.xt)
        self.zn_anim.append(self.zn)
        self.zt_anim.append(self.zt)
        self.p_anim.append(self.p_value)
        if not self.const_load:  # if solving variable load
            # add data to form u_contact_amin in 'calculation.py' module
            self.p_anim_variable[self.steps+1] = (self.force_inc_step, abs(self.p_value))
        self.steps += 1
        if self.write_excel:
            self.add_results_to_excel()

    def add_results_to_excel(self):
        move = self.excel_table_count * (self.table.shape[0] + 2)
        for i in range(self.table.shape[0]):                     # write table
            for j in range(self.table.shape[1]):
                self.worksheet.write(i + move + 2, j + 1, self.table[i][j])
        for i in range(self._basis.size):                        # write basis (right column) and min ratio
            self.worksheet.write(2 + move + i, 0, self._basis[i])
            if not np.isinf(self._min_ratio[i]):
                self.worksheet.write(2 + move + i, 1 + self.table.shape[1], self._min_ratio[i])
        for j in range(self.table.shape[1]-2):                     # write range (top row)
            self.worksheet.write(1 + move, 1 + j, j)
        # write leading column and row
        self.worksheet.write(move + 2, self.table.shape[1] + 2, self._leading_row)
        self.worksheet.write(move + 3, self.table.shape[1] + 2, self._leading_column)
        # write columns names
        self.worksheet.write(move + 1, self.table.shape[1] + 1, 'min ratio')
        self.worksheet.write(move + 1, self.table.shape[1] + 0, 'Rf')
        # write zn, xn, zt, xt
        self.worksheet.write(move + 2, self.table.shape[1] + 3, 'zn')
        self.worksheet.write(move + 2, self.table.shape[1] + 4, 'xn')
        self.worksheet.write(move + 2, self.table.shape[1] + 5, 'zt')
        self.worksheet.write(move + 2, self.table.shape[1] + 6, 'xt')
        for i in range(self.xn.shape[0]):
            self.worksheet.write(move + 3 + i, self.table.shape[1] + 3, self.zn[i])
            self.worksheet.write(move + 3 + i, self.table.shape[1] + 4, self.xn[i])
            self.worksheet.write(move + 3 + i, self.table.shape[1] + 5, self.zt[i])
            self.worksheet.write(move + 3 + i, self.table.shape[1] + 6, self.xt[i])

        self.excel_table_count += 1

    def _lemke_step(self):
        """
        One step of Lemke algorithm
        :return: None
        """
        # remember variables from this step
        self.table_previous = self.table.copy()
        self._basis_previous = self._basis.copy()
        self._leading_column_previous = self._leading_column
        # do step
        self._form_table()
        self._form_basis()
        self._form_leading_column()
        self._form_min_ratio()
        self._form_leading_row()
        # values from previous steps update
        self._basis = self._basis_next.copy()
        self._leading_column = self._leading_column_next
        self._leading_row = self._leading_row_next

    def lcp_solve(self):
        """
        Linear Complementary Problem solver
        :return: None
        """
        if self._trivial_solution() and self.const_load:
            print('LCP trivial solution')
            if not self.force_inc:
                self._results_anim()
                if self.write_excel:
                    self.workbook.close()
                return
            self._next_load_vector()
        for step in range(LEMKE_LIMIT_STEPS):
            self._results_anim()
            self._lemke_step()  # do the step
            # Do checks
            if not self._p_in_basis():
                if not self.force_inc:  # if solving only for const load
                    print("Normal solution of LCP const load in {} steps".format(step + 1))
                    self._results_anim()
                    break

                elif self.const_load:  # if stage was for const load
                    print("Normal solution of LCP const load in {} steps".format(step + 1))
                    self._next_load_vector()  # start force inc part (for variable load)
                    continue
                else:  # if stape was about solving variable load
                    print('Rough solution is occurred in force increment, you got rough results '
                          'as there in each contact pair the p (tightening weight) still exists\n'
                          'ultimate load parameter p={} has been found with ACCURACY_OF_LCP={} accuracy'
                          .format(self._get_p(), ACCURACY_OF_LCP))
                    if not self._rows_table * 2 + self.force_inc_step == self.table.shape[1] - 2:
                        # rollback = True because if p not in basis in force inc alg
                        # then it means there are all min_ration < 0
                        self._next_load_vector()  # continue to solve force increment
                        continue
                    break
            elif self._ray_solution():
                self._results_anim()
                if not self.force_inc:  # if we are solving only for const load
                    if not self._rough_solution():
                        print('Ray solution in LCP, you got the results for "almost broken" system,'
                              ' its {} step of Lemke, p={}'
                              .format(step, self._get_p()))
                    else:
                        print('Rough solution is occurred, you got rough results '
                              'as there in each contact pair the p (tightening weight) still exists\n'
                              'p={} is lesser than value ACCURACY_OF_LCP={} settled by user'
                              .format(self._get_p(), ACCURACY_OF_LCP))
                    break
                # if incrementation force algorithm than:
                # check if variables could be determined (leading element should not be zero or really close to it)
                self._leading_row = self._get_p_row()
                normal_ray_solution = not np.isclose(self.table[self._leading_row, self._leading_column], 0,
                                                 atol=ACCURACY_OF_LCP)
                # normal solution for force inc. is ray solution, but we can kick 'p' from basis and get zn,zn ect.
                if normal_ray_solution:
                    print(f'Ray solution of force increment algorithm'
                          f' on step: {self.steps}, OK!  with p:{self.p_value}'
                          f'\nFurther incrementation of the load will not change scheme')
                    self._results_anim(ray=True)  # find values and kick 'p' from basis
                else:
                    print(f'Ray solution in force increment algorithm, indeterminate variables / mechanism'
                          f' on step: {self.steps} with p >:{self.p_value}')
                # if there are still unsolved load vectors
                if not self._rows_table * 2 + self.force_inc_step == self.table.shape[1] - 2:
                    self._next_load_vector(normal_ray_solution)  # continue to solve force increment
                    continue
                # if number of leading column was chosen last time penultimate (предпоследняя) column
                break  # end

        if self.write_excel:
            self.workbook.close()
        if self.steps == LEMKE_LIMIT_STEPS-1:
            print(f"LCP solver got maximum number of steps, no objective results, p={self._get_p()}")
        return  # end

    def _next_load_vector(self, normal_ray_solution=True):
        """
        Modifying table to prepare it for using as usual table for increment force algorithm
        :param normal_ray_solution:
        :return:
        """
        if not normal_ray_solution:  # roll back to the previous version of the table, case it is impossible to kick 'p'
            self.table = self.table_previous.copy()
            self._basis = self._basis_previous.copy()
            self._leading_column = self._leading_column_previous
            self._results_anim(ray=True)  # find values and kick 'p' from basis
        # if there is more than 1 step with variable load, than use previous "react_vec - p column" as react vector
        if self.force_inc_step > 0:
            # new react_vec---previous react_vec-------------previous p column (active load vector)--------
            self.table[:, -1] = self.table[:, -1] - self.table[:, (self._rows_table * 2 + self.force_inc_step)] * self.p_value
        # self.table[:, -1] = self.table[:, -1] + 1  # add small artificial "pc" - "tightening weight"
        self.force_inc_step += 1  # choose next load vector
        # next we set table like it's first step
        self._leading_column = self._rows_table * 2 + self.force_inc_step  # choose new leading column
        self._leading_column_next = self._leading_column
        # choose new leading row using min ratio
        self._form_min_ratio()
        self._form_leading_row()
        # self._leading_row_next = np.argmax(self.table[:, self._leading_column])
        self._leading_row = self._leading_row_next
        self.const_load = False  # start/continue to solve for variable load
