import numpy as np
from typing import List
from input_data import SCHEME_THICKNESS
from FEM.element import ElementContainer
import weakref
import gc
from typing import Tuple
from numba import njit, prange
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_frame import ElementFrameContainer


"""
Here is Nodes, Stiffness Matrix, Load Vector
"""


class Node:
    """
    Class for creating nodes
    Each node has x, y fields
    """
    def __init__(self, x, y, indices, number, rotation=False):
        """
        Creates a node
        :param x: X Axis coordinate of the adding node (horizontal Axis)
        :param y: Y Axis coordinate of the adding node (vertical Axis)
        :param indices: list[], degrees of freedom (dof) of the node [dof1, dof2].
        When node added, new indices are added in system, but we can set specific indices to the node.
        :param number: number of the node. This value is the index of this object in NodeContainer.nodes_list field
        :param rotation: boolean, if the node has 3d degree of freedom - rotation.
        In that case indices of the node are become [dof1, dof2, dof3]
        """
        self.x = x
        self.y = y
        self.rotation = rotation
        self.indices = indices
        self.number = number  # IT ISN'T WORKING IF USE FRAGMENT of macro elements!
        # list of element that have this node
        self.elements = []
        # additional info about node
        self.dx = 0  # delta x (node dislocation) in linear problem
        self.dy = 0
        self.dx_c = 0  # delta x contact (node dislocation) in contact problem
        self.dy_c = 0


class NodeContainer:
    """
    Class for nodes in scheme. Contains a list of nodes objects
    Each node have x,y coordinates
    Each node have a number, it is matches (the same as) with index in list of nodes
    """
    instances = []  # all instances of the class will be here

    def __init__(self):
        # List of all nodes objects in scheme
        self.nodes_list = []
        # setting max degree of freedom equals -1, we have no any yet.
        self.max_dof = -1
        # adding node number. This number is equal to the index in nodes_list of future Node object
        self.number = -1
        # add a weakref to each instance so the garbage collector can clean up instances
        # when they're no longer needed.
        self.__class__.instances.append(weakref.proxy(self))

    def add_node(self, x: float, y: float, rotation=False):
        """
        Adds a node to a scheme
        :param x: x coordinate of a new node. horizontal
        :param y: y coordinate of a new node. vertical
        :param rotation: Is there a rotation degree of freedom in this node? False = no
        we need rotation in node for frame element (FEM.element_frame.ElementFrame())
        :return: None, void type
        """
        # if we set a rotation in this node, add third degree of freedom in it.
        if rotation is False:
            indices = [self.max_dof + 1, self.max_dof + 2]
            self.max_dof += 2
        else:
            indices = [self.max_dof + 1, self.max_dof + 2, self.max_dof+3]
            self.max_dof += 3
        # node is added, so increase the number of the node
        self.number += 1
        # create a node
        node = Node(x, y, indices, self.number, rotation)
        self.nodes_list.append(node)

        assert isinstance(rotation, bool), 'rotation must be boolean value'

    def __getitem__(self, node_number):
        """
        Node's list access through the indices at class
        :param node_number:number of node witch we need to access
        :return:Object of Node from nodes_list[node_number]
        """
        return self.nodes_list[node_number]

    def __len__(self):
        """
        Get amount of nodes in list
        :return:
        """
        return len(self.nodes_list)

    def __str__(self):
        """
        show nodes list in scheme
        :return:
        """
        string_to_print = ''
        i = 0
        for node in self.nodes_list:
            string_to_print += f'{i}:({node.x},{node.y}) '
            i += 1
        return string_to_print

    def delete(self, node_number):
        """
        Delete node from scheme
        :param node_number: number of the node to delete
        :return:
        """
        del(self.nodes_list[node_number])

    def find_nodes_numbers_along_segment(self, point1:Tuple, point2:Tuple, sorted_by_x:bool=True,
                                         sorted_by_y=True, relative_toletance:float=1e-5):
        """
        Finds all nodes along segment that are 'close' to the segment between point1 and point2
        :param point1:
        :param point2:
        :param sorted_by_x:
        :param sorted_by_y:
        :param relative_toletance:
        :return: list of nodes' numbers sorted in chosen way
        """
        # first form list of Node objects that have needed coordinates and list of numbers
        list_of_nodes = []
        list_of_numbers = []
        part_of_cp1 = point2[0] - point1[0]
        part_of_cp2 = point2[1] - point1[1]
        for i, node in enumerate(self.nodes_list):
            # The absolute value of the cross product is twice the area of the triangle formed by the three points
            cross_product = (node.y - point1[1]) * part_of_cp1 - (node.x - point1[0]) * part_of_cp2
            if np.isclose(cross_product, 0, rtol=relative_toletance) and\
                    min(point1[0], point2[0]) <= node.x <= max(point1[0], point2[0]) and\
                    min(point1[1], point2[1]) <= node.y <= max(point1[1], point2[1]):
                list_of_nodes.append(node)
                list_of_numbers.append(i)
        # now form sorted list of found nodes, sort by attribute wished
        if sorted_by_x is True and sorted_by_y is True:
            sorted_list_of_nodes_numbers = [number for (number, node) in sorted(zip(list_of_numbers, list_of_nodes),
                                                                                key=lambda pair: (
                                                                                pair[1].x, pair[1].y))]
        elif sorted_by_x is True:
            sorted_list_of_nodes_numbers = [number for (number, node) in sorted(zip(list_of_numbers, list_of_nodes),
                                                                                key=lambda pair: pair[1].x)]
        elif sorted_by_y is True:
            sorted_list_of_nodes_numbers = [number for (number, node) in sorted(zip(list_of_numbers, list_of_nodes),
                                                                                key=lambda pair: pair[1].y)]
        else:  # if sort is not needed
            sorted_list_of_nodes_numbers = list_of_nodes
        return sorted_list_of_nodes_numbers


class StiffnessMatrix:
    """
    All field and methods about stiffness matrix
    """
    def __init__(self, nodes, el_frame, el_4node, el_null):
        self._add_inidices_to_elements(nodes, el_frame, el_4node, el_null)
        self.nodes = nodes
        self.el_frame = el_frame
        self.el_4node = el_4node
        self.el_null = el_null
        max_dof = find_max_degree_of_freedom()
        # Global stiffness matrix
        self.r = np.zeros([max_dof + 1, max_dof + 1], dtype=float)
        self.r = np.zeros(shape=(max_dof + 1, max_dof + 1), dtype=float)

        self.supports = np.array([], dtype=int)

    def _add_inidices_to_elements(self, nodes, el_frame, el_4node, el_null):
        """
        Add indices to MI (matrices of indices) in each element in scheme
        :param nodes: nodes in scheme. FEM.scheme.NodeContainer
        :param el_frame: frame elements in scheme. FEM.element_frame.ElementFrameContainer
        :param el_4node: 4 node elements in scheme. FEM.element_4node.Element4NodeContainer
        :param el_null: null-elements in scheme. FEM.element_null.ElementNullContainer
        :return:
        """
        # adding indices to 4node elements
        if el_4node is not None:
            for element in el_4node:
                for node_num in element.EN:
                    element.MI += nodes[node_num].indices[:2]
        # adding indices to frame elements
        if el_frame is not None:
            for element in el_frame:
                for node_num in element.EN:
                    if nodes[node_num].rotation == False:
                        new_dof = find_max_degree_of_freedom() + 1  # add extra dof for rotation if there is not
                        nodes[node_num].rotation = True
                        nodes[node_num].indices.append(new_dof)
                        print('Rotation for node {} was added. New {} degree of freedom'.format(node_num, new_dof))
                    element.MI += nodes[node_num].indices[:3]
        # adding indices to null-elements
        if el_null is not None:
            for element in el_null:
                for node_num in element.EN:
                    element.MI += nodes[node_num].indices[:2]
                element.MI.append(find_max_degree_of_freedom() + 1)  # add extra dof for central node in null-element

    def ensemble_ke_gsc_to_r(self, element, ke_gsc):
        """
        Ensembles stiffness matrix of the element in global system of coordinates to Global Stiffness Matrix
        :param element: object of the element that we want to ensemble
        :param ke_gsc: stiffness matrix of the element in global system of coordinates
        :return: None, void
        """
        for i in range(len(element.MI)):
            for j in range(len(element.MI)):
                self.r[element.MI[i], element.MI[j]] += ke_gsc[i, j]

    def form_r(self):
        """
        Ensemble stiffness matrix (Do the Global Stiffness Matrix (GSM))
        :return: None
        """
        # Cycle for frame elements
        if self.el_frame is not None:
            for element in self.el_frame:
                ke_gsc = element.form_ke_gsc()
                self.ensemble_ke_gsc_to_r(element, ke_gsc)

        # Cycle for 4node elements
        if self.el_4node is not None:
            for element in self.el_4node:
                ke4 = element.form_ke(nodes=self.nodes)
                self.ensemble_ke_gsc_to_r(element, ke4)

        # Cycle for null-elements
        if self.el_null is not None:
            for element in self.el_null:
                ke_gsc = element.form_ke_gsc()
                self.ensemble_ke_gsc_to_r(element, ke_gsc)
        self.supports_to_r()  # add supports to global stiffness matrix

    def supports_to_r(self):
        """
        Add supports to the scheme.
        :return: None, modifies self.R of obj of FEM.scheme.StiffnessMatrix() class
        """
        max_dof = find_max_degree_of_freedom()
        for i in self.supports:
            self.r[:, i] = np.zeros(max_dof+1)  # set i column equal to zeros
            self.r[i] = np.zeros(max_dof+1)  # set i row equal to zeros
            self.r[i, i] = 1  # set 1 in a cross


    def support_nodes(self, list_of_nodes:List[int], direction:str='hv'):
        """
        Adding supports (boundary conditions) to the system node-wise.
        numpy 2d array [n,2], n - number of supported degrees of freedom (dof).
        ex: support = [[sup_dof_1, sup_dof_2, ... , sup_dof_n], [0, 0, ... , 0].
        Second column used to model buckling of supports
        :param list_of_nodes: iterable 1D list of node numbers in scheme.
        :param direction: if not set, it use 'hv' - to add horizontal and vertical supports in node
        'h' - add only horizontal support (hinged) in node
        'v' - add only vertical support (hinged) in node
        'r' - rotation support in node
        :return: None, adding supports to stiffness matrix 'supports' array
        """
        sup_amount = len(direction)*len(list_of_nodes)
        supports = np.zeros(sup_amount, dtype=int)
        i = 0
        for node_number in list_of_nodes:
            for d in direction:
                if d == 'h':
                    # its degree of freedom in 0th place in node.MI array
                    supports[i] = self.nodes[node_number].indices[0]
                elif d == 'v':
                    # its degree of freedom in 1st place in node.MI array
                    supports[i] = self.nodes[node_number].indices[1]
                elif d == 'r':
                    # its degree of freedom in 2nd place in node.MI array
                    supports[i] = self.nodes[node_number].indices[2]
                i += 1
        self.supports = np.concatenate((self.supports, supports), axis=0)  # add information to the object

        assert isinstance(direction, str), 'direction must be string'
        assert len(direction) <= 3 and len(direction) > 0, 'direction must not be longer than 3 chars. EX: "hv" or "vh" + "r". And must be not empty'
        assert isinstance(list_of_nodes, list), 'list of nodes must be list'

    def add_spring(self, degree_of_freedom, stiffness):
        """
        Adding spring to the stiffness matrix directly. 
        Simply adding some stiffness value to the r[i, i] on main diagonal
        :param degree_of_freedom: dof where additional stiffness will be applied
        :param stiffness: 
        :return: 
        """
        self.r[degree_of_freedom, degree_of_freedom] += stiffness


class LoadVector:
    """
    Load vector of the system
    """

    def __init__(self, element_null=None, vectors_amount=None):
        """
        Initialize load vector
        :param element_null: null_elements in scheme. Use if we making load vector for forming Contact Stiffness Matrix
        :param vectors_amount: need number of steps if we forming deformed scheme on every step
        FEM.element_null.ElementNullContainer() class object
        """
        # maximum degree of freedom
        max_dof = find_max_degree_of_freedom()
        self.rf = []
        # if we create load vector(matrix) for CSM
        if element_null is not None:  # creating 2d array for creating csm
            if len(element_null) > 0:  # check if we have null elements
                self.rf = np.zeros([max_dof + 1, len(element_null)],
                                   dtype=float)
            else:
                raise ValueError("You can't create load vector for_csm without null_elements!")
        elif vectors_amount is not None:
            self.rf = np.zeros([max_dof + 1, vectors_amount], dtype=float)
        # if we creating simple load vector
        else:
            self.rf = np.zeros(max_dof + 1, dtype=float)
        self.supports = []

    def support(self):
        """
        Add edge conditions (supports) to the load vector. Make supported degree of freedoms values  = 0
        :return:
        """
        if len(self.rf.shape) > 1:  # if array have more than 1 column
            for i in self.supports:
                self.rf[i] = np.zeros(self.rf.shape[1])
        else:
            for i in self.supports:
                self.rf[i] = 0


    def add_own_weight_to_rf(self, nodes_scheme, element_container_list=None):
        """
        Adds own weight of the elements to Global Load Vector RF
        :param nodes_scheme: Object of the nodes in scheme
        :param element_container_list: list of element containers. Could set of:
        el_4node_linear: Object of FEM.element_4node.Create4NodeElement
        el_frame: Object of FEM.element_frame.CreateFrameElement
        :return: None, void
        """
        for element_container in element_container_list:
            for element in element_container:
                mass_in_nodes_of_el = element.mass_part_in_nodes(nodes_scheme)
                counter = 0
                for node_num in element.EN:
                    index = nodes_scheme[node_num].indices[1]
                    self.rf[index] += mass_in_nodes_of_el[counter] * element.own_weight
                    counter += 1

        assert isinstance(element_container_list, list), 'element_container_list must be a list!'

    def add_concentrated_force(self, force:[int, float], degree_of_freedom:int, vector_num:int=None):
        """
        Func for adding concentrated force to load vector
        :param force: the value of the external force F, [N]
        :param degree_of_freedom: number of degree of freedom where we should put on the force
        :param vector_num: number of the load vector if we have multiple rf in LoadVector
        :return: None
        """
        assert isinstance(degree_of_freedom, int), 'degree of freedom must be integer'
        if vector_num is not None:
            assert len(self.rf.shape) != 1, f'you need to add multiple rf in LoadVector'
            assert isinstance(vector_num, int), 'vector number variable when adding force must be an integer'
            self.rf[degree_of_freedom][vector_num] += force
        else:
            self.rf[degree_of_freedom] += force



    def _add_Qe_to_RF(self, element, Qe, for_csm=False, el_null_num=None, lemke_step_num=None):
        """
        Adding element's global load vector to Load vector of the system
        :param element: element which Qe that we want to add to RF (load vector for the scheme)
        :param Qe: global load vector of the element
        :param for_csm: if we adding null_elements load (dislocation) to form ContactStiffnessMatrix we set for_csm=True
        :param el_null_num: number of the null-element which load we want to add to the LoadVector (RF)
        :param: lemke_step_num: if we forming load vector for visualizing deformed scheme
        :return: None
        """
        k = 0
        if for_csm:
            for dof in element.MI:
                self.rf[dof, el_null_num] += Qe[k]
                k += 1
        elif lemke_step_num is not None:
            for dof in element.MI:
                self.rf[dof, lemke_step_num] += Qe[k]
                k += 1
        else:
            for dof in element.MI:
                self.rf[dof] += Qe[k]
                k += 1

def add_displacement(degree_of_freedom, stiffness_matrix, load_vector):  # TODO make this to work, end code write
    """
    Adding certain displacement
    :param degree_of_freedom:
    :param stiffness_matrix:
    :param load_vector:
    :return:
    """
    # modifying load vector
    if len(load_vector.rf.shape) > 1:  # if array have more than 1 column
        load_vector.rf[degree_of_freedom] = np.zeros(load_vector.rf.shape[1])
    else:
        load_vector.rf[degree_of_freedom] = 0
    # modifying stiffness matrix




def find_max_degree_of_freedom():
    """
    Finds maximum degree of freedom in scheme
    :return: degree of freedom, int
    """
    # this part for searching max dof if in scheme there are no elements (only nodes)
    max_dof = 0
    for instance in NodeContainer.instances:
        for node in instance:
            max_dof_on_node = max(node.indices)
            if max_dof_on_node > max_dof:
                max_dof = max_dof_on_node
    for obj in gc.get_objects():
        if isinstance(obj, ElementContainer):
            for element in obj.elements_list:
                max_dof_el = np.amax(element.MI, initial=0)
                if max_dof_el > max_dof:
                    max_dof = max_dof_el
    return max_dof


def list_of_all_elements():
    """
    Collecting all elements is scheme. Taking instances from FEM.element.ElementCollector() class
    :return: list of all elements as objects
    """
    all_elements = []
    for obj in gc.get_objects():
        if isinstance(obj, ElementContainer):
            for element in obj.elements_list:
                # collecting ALL elements in scheme
                all_elements.append(element)
    return all_elements

# also possible if use "prange" from numba it can be:
# @njit(fastmath=True, parallel=True) #parallel calculations
#@njit(nopython=True, fastmath=True, parallel=True)  # it calculates slower on numba
# def linalg_solve(a, b):
#     return np.linalg.solve(a, b)
#
# def solve_slae(sm, lv):
#     """
#     Solve system of linear equations. Find global displacements.
#     :param sm: stiffness matrix
#     :param lv: load vector
#     :return: vector of global displacements
#     """
#     sm.form_r()
#     lv.supports = sm.supports[:]
#     lv.support()
#     return linalg_solve(sm.r, lv.rf)
#
# def form_u_contact(sm, lv, zn, zt, element_null):
#     """
#     Forming vector of global displacements for unilateral frictional contact
#     :param sm: global stiffness matrix. FEM.scheme.StiffnessMatrix() class object
#     :param lv: global load vector. FEM.scheme.LoadVector() class object
#     :param zn: vector of vectors of normal to the contact zone mutual displacements. LCP.lemke.Lemke.zn attribute
#     :param zt: vector of vectors of tangent to the contact zone mutual displacements. LCP.lemke.Lemke.zt attribute
#     :param element_null:
#     :return:
#     """
#     # create object of LoadVector class with amount of lemke steps width (vectors_amount = lemke step amount)
#     lv_contact = LoadVector(vectors_amount=len(zn))
#     for column in range(lv_contact.rf.shape[1]):
#         lv_contact.rf[:, column] = lv.rf.copy()  # copy existing values of external load
#     lv_contact.supports = lv.supports.copy()  # copy supports
#     for i in range(len(zn)):
#         zn_temp_i = zn[i].copy()  # copy, we need initial values to visualize mutual displacements
#         for k in range(len(zn_temp_i)):
#         # subtract gap_length from mutual displacement to ignore it, while forming global nodes displacements
#             if element_null[k].gap_length != 0:
#                 zn_temp_i[k] -= element_null[k].gap_length
#             z_nt_i = np.concatenate([zn_temp_i, zt[i]])
#         for el_null, z in zip(element_null, z_nt_i):
#             Qe = el_null.form_Qe(displacement=z)
#             lv_contact._add_Qe_to_RF(el_null, Qe, lemke_step_num=i)  # add load from null elements to load vector
#     lv_contact.support()
#     solution = linalg_solve(sm.r, lv_contact.rf)
#     result = []
#     for column in range(solution.shape[1]):  # make it look like list of lists
#         result.append(solution[:, column])
#     return result
