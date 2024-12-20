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
        if lv_const is None:
            self.lv_const = FEM.scheme.LoadVector()
        self.element_frame = element_frame
        if isinstance(element_container_obj, Element4NodeLinearContainer) and element_container_obj.elements_list:
            self.element_4node = element_container_obj
        else: self.element_4node = None
        self.element_null = element_null
        self.lv_variable = lv_variable
        if lv_variable is None:  # if there was no variable load vector, create and make it nested
            self.lv_variable = FEM.scheme.LoadVector()
        # Variables to calculate
        self.u_linear_const = None  # linear global nodes displacements from constant load
        self.u_linear_variable = None  # linear global nodes displacements from variable load
        self.intl_table = None
        self.lemke = None
        self.u_contact_anim = None  # global nodes displacements as nested array for each lemke's step

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
        if len(self.u_linear_variable.shape) == 1:
            self.u_linear_variable = [self.u_linear_variable]
        else:
            tmp_arr = []
            for i in range(self.u_linear_variable.shape[1]):
                tmp_arr.append(self.u_linear_variable[:, i])
            self.u_linear_variable = tmp_arr


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
        start = time.time()
        self.lemke = Lemke(self.intl_table)
        self.lemke.lcp_solve()
        print(f"LCP solved in:{time.time() - start} sec.")

    def form_u_contact(self):
        """
        Forming vector of global displacements for unilateral (frictional) contact
        :return:
        """
        zn = self.lemke.zn_anim
        zt = self.lemke.zt_anim
        if self.lv_variable.rf.ndim == 1:
            self.lv_variable.rf = np.expand_dims(self.lv_variable.rf, 1)
        # create object of LoadVector class with amount of lemke steps width (vectors_amount = lemke step amount)
        lv_contact = FEM.scheme.LoadVector(vectors_amount=len(zn))
        rf_addition = np.zeros(shape=(lv_contact.rf.shape[0],))  # create zeros empty list
        p_value = 0
        stage_num = 0  # to remember if it is next lv_var is going
        for step_num in range(lv_contact.rf.shape[1]):
            lv_contact.rf[:, step_num] = self.lv_const.rf.copy()  # copy existing values of external
            if step_num in self.lemke.p_anim_variable:  # if this step_num is for variable load (checking keys)
                lv_num = self.lemke.p_anim_variable[step_num][0] - 1  # number of stage in multiple load vectors
                if stage_num < lv_num:  # if it is next lv_var is going
                    rf_addition += self.lv_variable.rf[:, stage_num] * p_value
                    stage_num = lv_num
                p_value = self.lemke.p_anim_variable[step_num][1]
                lv_contact.rf[:, step_num] += (self.lv_variable.rf[:, lv_num] * p_value)  # add values from lv_variable
                lv_contact.rf[:, step_num] += rf_addition  # add previous solved variable load vector
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
        for step_num in range(solution.shape[1]):  # make it look like list of lists
            result.append(solution[:, step_num])
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
        print('Linear problem solved')
        self.form_initial_table()
        print('Initial table ready')
        self.solve_lcp()
        self.form_u_contact()
        self._add_displacement_data_to_nodes()
        print('Data to GUI ready')
        print(f'Linear problem size: {self.sm.r.shape[0]}')
        print(f'Nonlinear problem size: {self.intl_table.table.shape[0]}')

        end = time.time()
        print("Calculation Time: ", end - start)
