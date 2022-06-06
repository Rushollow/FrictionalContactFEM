import pyqtgraph as pg
from FEM.element_null import ElementNull
from FEM.element_4node import Element4NodeLinear, Element4NodeLinearContainer
from SchemeForm.macro_element import ElementMacro
from FEM.scheme import list_of_all_elements, NodeContainer
from input_data import SCALE_DEF, FRICTION_COEFFICIENT
from calculation import Calculate
import numpy as np



# region Shape functions for frame element
def _frame_n1(x, L):
    return 1 - 3 * ((x * x) / (L * L)) + 2 * ((x ** 3) / (L ** 3))


def _frame_n2(x, L):
    return x - 2 / L * (x * x) + 1 / (L * L) * (x ** 3)


def _frame_n4(x, L):
    return 3 * (x * x) / (L * L) - 2 * (x ** 3) / (L ** 3)


def _frame_n5(x, L):
    return - x * x / L + x ** 3 / (L * L)


def _frame_agg(x, L, v1, fi1, v2, fi2):
    """
    Aggregation of all frame element form functions.
    :param x: axis along the element
    :param L: length of the element
    :param v1: vertical local displacement of the first node (of the element)
    :param fi1: local rotation of the first node (of the element)
    :param v2: vertical local displacement of the second node (of the element)
    :param fi2: local rotation of the second node (of the element)
    :return: deformed form of the frame element, depend on the parameters
    """
    return v1 * _frame_n1(x, L) + fi1 * _frame_n2(x, L) + \
           v2 * _frame_n4(x, L) + fi2 * _frame_n5(x, L)


# endregion


def _modify_coordinates_for_plotting_4node(x, y):
    """
    Moving lines of sides of the element to the center for good looking
    :param x: vector of horizontal x coordinates
    :param y: vector of vertical y coordinates
    :return: x, y new vectors
    """
    avg_x = np.mean(x)
    avg_y = np.mean(y)
    x = [i + (avg_x - i) / 100 for i in x]
    y = [i + (avg_y - i) / 100 for i in y]
    return x, y


class PlotScheme(Calculate):
    """
    Class for plotting. methods for schemes, results.
    """

    def __init__(self, nodes, sm, partition=10, scale_def=SCALE_DEF,
                 text=True, autorun=False,
                 *args, **kwargs):
        """
        Initialize object to work with graphics. Form data to show for user (make arrays for QT)
        :param lv_variable: load vector used for Force increment Lemke's algorithm
        :param partition: division by n parts the frame elements to get smother curve while deformed.
        :param scale_def: scale of global displacements.
        :param text: bool value if text values should be visualized on plot
        :param autorun: autoruns calculation of the given problem before loading UI.
        May be used to get some calculation data in console through PlotScheme.lemke.*some_data or some other
        """
        # initialize parent class
        super().__init__(nodes, sm, *args, **kwargs)
        # parameters (attributes)
        self.partition = partition
        self.scale = scale_def
        self.text = text
        self.size_supports = 100
        self.size_force = 200
        self.scale_lcp = None

        # add arrays for data to plot
        self.arr_nodes_pos = None
        self.arr_nodes_pos_deformed = None
        self.arr_nodes_pos_deformed_frame = None
        self.arr_nodes_pos_contact = None
        self.arr_null_el_1st_nodes_pos_angle = None
        self.arr_supp_pos_angle = None
        self.arr_frame_en = None
        self.arr_frame_en_deformed = None
        self.arr_nodes_pos_contact_frame = None
        self.arr_4node_en = None
        self.arr_4node_en_deformed = None
        self.initialize_arrays()  # form arrays for plotting (from results data)

        if autorun:
            self.calculate_all()

    def initialize_arrays(self):
        """
        Creating and filling all data arrays with zeros with specific shape.
        Creates np.array for plot data if such data is needed (if elements exist in scheme)
        :return:
        """
        # forming 2d numpy arrays:
        self.arr_nodes_pos = np.zeros((len(self.nodes), 2))  # for all nodes in scheme positions
        # for deformed scheme
        self.arr_nodes_pos_deformed = np.zeros((len(self.nodes), 2))  # for all nodes in scheme positions deformed
        # unilateral contact deformed
        self.arr_nodes_pos_contact = np.zeros((len(self.nodes), 2))  # for all nodes in scheme positions deformed
        # for all first nodes of null_elements
        self.arr_null_el_1st_nodes_pos_angle = np.zeros((len(self.element_null), 3))
        self.arr_supp_pos_angle = np.zeros((len(self.sm.supports), 3))  # for all supports


        if self.element_frame is not None:  # for frame
            self.arr_frame_en = np.zeros((len(self.element_frame), 2), dtype=int)  # for all ElementNodes (EN)
            # linear deformation
            self.arr_nodes_pos_deformed_frame = np.zeros((len(self.element_frame) * (self.partition + 1), 2))
            self.arr_frame_en_deformed = np.zeros((len(self.element_frame) * self.partition, 2), dtype=int)
            # contact deformed scheme
            self.arr_nodes_pos_contact_frame = np.zeros((len(self.element_frame) * (self.partition + 1), 2))
        # for all ElementNodes (EN) for 4node_element
        if self.element_4node is not None:  # for 4node element
            self.arr_4node_en = np.zeros((len(self.element_4node) * 4, 2), dtype=int)  # for scheme
            self.arr_4node_en_deformed = np.zeros((len(self.element_4node) * 4, 2), dtype=int)  # for deformed scheme

    def fill_arrays_scheme(self):
        """
        fill the arrays with data to plot scheme
        :return: None
        """
        # add nodes coordinates to single arr
        for i, node in enumerate(self.nodes):
            self.arr_nodes_pos[i] = np.array([node.x, node.y])
        # add pos and angle info about supports
        for i, dof in enumerate(self.sm.supports):
            for node in self.nodes:
                if dof in node.indices:
                    indx = node.indices.index(dof)
                    if indx == 0: angle = 180
                    elif indx == 1: angle = 90
                    else: angle = -45
                    self.arr_supp_pos_angle[i] = np.array([node.x, node.y, angle])
        # if there are null elements in scheme - add info about their first nodes positions and angle
        # to arr (x, y, angle)
        if self.element_null is not None:
            for i, el_null in enumerate(self.element_null):
                x, y = el_null.nodes_coordinates(self.nodes)
                self.arr_null_el_1st_nodes_pos_angle[i] = np.array([x[0], y[0], np.degrees(el_null.alpha)])
        # add frame ElementNodes info to the array
        if self.element_frame is not None:
            self.arr_frame_en = np.array(self.element_frame.EN, dtype=int)
        # add 4node element ElementNodes info to the array
        if self.element_4node is not None:
            tmp_arr = np.array(self.element_4node.EN)
            rows, cols = tmp_arr.shape
            # add first nodes to the right of array to make 4nodeElements complete
            tmp_arr = np.concatenate((tmp_arr, np.reshape(tmp_arr[:, 0], (rows, 1))), axis=1)
            for i in range(rows):
                for j in range(cols):
                    self.arr_4node_en[i * 4 + j] = tmp_arr[i, j:j + 2]

    def fill_arrays_deformed_linear(self):
        """
        Add data to arrays to display info to user.
        Filling the info for linear deformation scheme and nodes to connect (frame is divided into multiple segments)
        :return: nothing
        """
        # get coordinates of the nodes of the element that we want to draw
        nodes_def_linear = self.__create_new_deformed_nodes(self.u_linear_const)
        # add nodes coordinates to single arr deformed
        for i, node_def in enumerate(nodes_def_linear):
            self.arr_nodes_pos_deformed[i] = np.array([node_def.x, node_def.y])
            self.nodes[i].dx = node_def.x - self.nodes[i].x
            self.nodes[i].dy = node_def.y - self.nodes[i].y
        if self.element_frame is not None:
            n = 0  # node number
            for i, one_frame_element in enumerate(self.element_frame):  # iterate over all frame elements
                x, y = self.__create_points_for_deformed_frame(one_frame_element, self.u_linear_const)  # def form by points
                for j in range(len(x)):  # could choose either range(len(y))
                    # linear deformation data
                    self.arr_nodes_pos_deformed_frame[i*(self.partition+1) + j] = np.array([x[j], y[j]])  # add position
                for j in range(len(x)-1):  # could choose either range(len(y))
                    # add element nodes
                    self.arr_frame_en_deformed[i * self.partition + j] = np.array([n, n + 1], dtype=int)
                    n += 1  # increment node number
                n += 1  # new element has new node number (in logic of partition each element)

    def fill_arrays_deformed_contact(self, i_step=-1):
        """
        Add data to arrays to display info to user
        :param i_step: step of Lemke's algorithm. default -1 is the last step - final results
        :return: nothing
        """
        # get coordinates of the nodes of the element that we want to draw
        nodes_def_contact = self.__create_new_deformed_nodes(self.u_contact_anim[i_step])
        # add nodes coordinates to single arr deformed
        for i, node in enumerate(nodes_def_contact):
            self.arr_nodes_pos_contact[i] = np.array([node.x, node.y])
        if self.element_frame is not None:
            for i, one_frame_element in enumerate(self.element_frame):  # iterate over all frame elements
                xc, yc = self.__create_points_for_deformed_frame(one_frame_element, self.u_contact_anim[i_step])  # def contact form by points
                for j in range(len(xc)):  # could choose either range(len(y))
                    # contact nonlinear data
                    self.arr_nodes_pos_contact_frame[i * (self.partition + 1) + j] = np.array([xc[j], yc[j]])

    def __create_new_deformed_nodes(self, u):
        """
        Creates FEM.scheme.NodeContainer() object with new nodes positions in deformable scheme
        :param u: vector of global displacements. Result of solving SLAE, Ex: numpy.linalg.solve(sm.R, lv.RF)
        :return: NodeContainer() object with new nodes
        """
        new_nodes = NodeContainer()
        for node in self.nodes:
            move_horizontal = u[node.indices[0]] * self.scale
            move_vertical = u[node.indices[1]] * self.scale
            x = node.x + move_horizontal
            y = node.y + move_vertical
            new_nodes.add_node(x, y, rotation=node.rotation)
        return new_nodes

    def __create_points_for_deformed_frame(self, one_frame_el, u):
        """
        Plot foo from point1 to point2
        :param one_frame_el: frame element object of FEM.element_frame.ElementFrame() class
        :param u: vector of global nodes displacements (result of numpy.linalg.solve(r,lv))
        :return: tuple ([x], [y]) nodes coordinates
        """
        # get nodes coordinates of the frame element
        [x1_init, _], [y1_init, _] = one_frame_el.nodes_coordinates(self.nodes)  # x2_init, y2_init are _, _
        # displacements of the frame element nodes, global
        don = one_frame_el.global_displacements(u)
        x1def = don[0] * self.scale  # horizontal 1st node
        x2def = don[3] * self.scale  # horizontal 2nd node
        y1def = don[1] * self.scale  # vertical 1st node
        y2def = don[4] * self.scale  # vertical 2nd node
        # global node displacements, rotation
        fi1 = don[2] * self.scale  # rotation of the node in the beginning of the element (1st node)
        fi2 = don[5] * self.scale  # rotation of the node in the ending of the element (2nd node)
        # vertical displacements on the frame element, local
        v1 = y1def * one_frame_el.cosa - x1def * one_frame_el.sina
        v2 = y2def * one_frame_el.cosa - x2def * one_frame_el.sina
        # horizontal displacements on the frame element, local
        w1 = x1def * one_frame_el.cosa + y1def * one_frame_el.sina
        w2 = x2def * one_frame_el.cosa + y2def * one_frame_el.sina
        # length of the element after deformation
        L = one_frame_el.length + w2 - w1
        if L <= 0:
            print(f'Element frame number: {one_frame_el.number} has scaled deformations so big that its length <= 0,\n'
                  f'please reduce SCALE_DEF parameter or make external load lower,\n'
                  f'results may be unreadable')
        # initial coordinate system. Create array of points to rotate and move
        xi = np.linspace(0, L, self.partition+1)
        yi = _frame_agg(xi, L, v1, fi1, v2, fi2)  # using 4 form functions of frame element
        # modified coordinate system. Create array of points to plot
        x_modified = one_frame_el.cosa * xi - one_frame_el.sina * yi + x1_init + w1 * one_frame_el.cosa
        y_modified = one_frame_el.sina * xi + one_frame_el.cosa * yi + y1_init + w1 * one_frame_el.sina
        return x_modified, y_modified

