import numpy as np
import time
import FEM.scheme
from FEM.element_4node import Element4NodeLinearContainer
from LCP.lemke import Lemke
from LCP.initial_table import InitialTable

import xlsxwriter  # need this to pass initial table to excel


class Calculate:
    """
    Class for calculation problems
    """
    def __init__(self, nodes, sm, lv_const, lv_variable=None,
                 element_null=None, element_frame=None, element_container_obj=None,
                 *args, **kwargs):
        """
        Constructor
        :param nodes: FEM.scheme.NodeContainer
        :param sm: FEM.scheme.StiffnessMatrix
        :param lv_const: FEM.scheme.LoadVector (constant part of load)
        :param lv_variable: variable load vector for force increment algorithm. FEM.scheme.LoadVector
        :param element_null: FEM.element_null.ElementNullContainer
        :param element_frame: FEM.element_frame.ElementFrameContainer
        :param element_container_obj: usually FEM.element_4node.Element4NodeContainer
        :param args: maybe other elements?
        :param kwargs:
        """
        # data set by user (scheme):
        self.nodes = nodes
        self.sm = sm
        self.lv_const = lv_const
        self.element_frame = element_frame
        if isinstance(element_container_obj, Element4NodeLinearContainer):
            self.element_4node = element_container_obj
        else: self.element_4node = None
        self.element_null = element_null
        if lv_variable is None:  # if there was no variable load vector, create and make it nested
            self.lv_variable = FEM.scheme.LoadVector()
            self.lv_variable.rf = [self.lv_variable.rf]
        else:
            self.lv_variable = lv_variable
            if len(lv_variable.rf.shape) == 1:  # check if there is only one variable load vector and make it nested
                self.lv_variable.rf = [lv_variable.rf]
            else:  # if more than 1 -> make it nested list of arrays
                temp_arr = []
                print(lv_variable.rf[0])
                for col in range(lv_variable.rf.shape[1]):
                    temp_arr.append(lv_variable.rf[:, col])
                self.lv_variable.rf = temp_arr

        # Variables to calculate
        self.u_linear_const = None  # linear global nodes displacements from constant load
        self.u_linear_variable = None  # linear global nodes displacements from variable load
        self.intl_table = None
        self.lemke = None
        self.u_contact_anim = None

    def solve_slae(self):
        """
        Solve system of linear equations. Find global displacements.
        :return: vector of global displacements
        """
        self.sm.form_r()
        self.lv_const.supports = self.sm.supports
        self.lv_const.support()
        self.u_linear_const = np.linalg.solve(self.sm.r, self.lv_const.rf)
        self.u_linear_variable = np.linalg.solve(self.sm.r, self.lv_variable.rf)

    def form_initial_table(self):
        """
        Forming initial table for Lemke's algorithm
        :return:
        """
        self.intl_table = InitialTable(self.element_null, self.sm, self.lv_const, self.u_linear_const,
                                       self.u_linear_variable, self.lv_variable)
        self.intl_table.form_initial_table()

    def solve_lcp(self):
        """
        Solving LCP (by Lemke's algorithm)
        :return:
        """
        self.lemke = Lemke(self.intl_table)
        self.lemke.lcp_solve()

    def form_u_contact(self):
        """
        Forming vector of global displacements for unilateral frictional contact
        :return:
        """
        zn = self.lemke.zn_anim
        zt = self.lemke.zt_anim
        # create object of LoadVector class with amount of lemke steps width (vectors_amount = lemke step amount)
        lv_contact = FEM.scheme.LoadVector(vectors_amount=len(zn))
        for column in range(lv_contact.rf.shape[1]):
            lv_contact.rf[:, column] = self.lv_const.rf.copy()  # copy existing values of external load
        lv_contact.supports = self.lv_const.supports.copy()  # copy supports
        for i in range(len(zn)):
            zn_temp_i = zn[i].copy()  # copy, we need initial values to visualize mutual displacements
            for k in range(len(zn_temp_i)):
                # subtract gap_length from mutual displacement to ignore it, while forming global nodes displacements
                if self.element_null[k].gap_length != 0:
                    zn_temp_i[k] -= self.element_null[k].gap_length
                z_nt_i = np.concatenate([zn_temp_i, zt[i]])
            for el_null, z in zip(self.element_null, z_nt_i):
                Qe = el_null.form_Qe(displacement=z)
                lv_contact._add_Qe_to_RF(el_null, Qe, lemke_step_num=i)  # add load from null elements to load vector
        lv_contact.support()
        solution = np.linalg.solve(self.sm.r, lv_contact.rf)
        result = []
        for column in range(solution.shape[1]):  # make it look like list of lists
            result.append(solution[:, column])
        self.u_contact_anim = result

    def _add_displacement_data_to_nodes(self):
        """
        Adding info about nodes displacements (dx, dy, dx_c....) to nodes objects
        :return: nothing
        """
        u_contact = self.u_contact_anim[-1]
        for node in self.nodes:
            node.dx = self.u_linear_const[node.indices[0]]  # linear def info
            node.dy = self.u_linear_const[node.indices[1]]
            node.dx_c = u_contact[node.indices[0]]  # contact def info
            node.dy_c = u_contact[node.indices[1]]

    def calculate_all(self):
        """
        Calculating all data step by step
        :return:
        """
        start = time.time()
        # calculate
        self.solve_slae()
        self.form_initial_table()
        self.solve_lcp()
        self.form_u_contact()
        self._add_displacement_data_to_nodes()

        end = time.time()
        print("Calculation Time: ", end - start)

    def table_to_excel(self, file_name='Initial_table'):
        '''
        this function is needed fot sending initial table to the Excel file
        :return:
        '''
        self.solve_slae()
        self.form_initial_table()
        self.solve_lcp()

        array = self.intl_table.table
        workbook = xlsxwriter.Workbook(file_name)
        worksheet = workbook.add_worksheet()
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                worksheet.write(i, j, array[i][j])
        #worksheet.write_column(0, array.shape[1]+4, ['zn']+list(self.lemke.zn_anim))
        workbook.close()
