import numpy as np
from typing import List
from FEM.element import ElementContainer
from FEM.element import ElementMethods
import math

class ElementFrameContainer(ElementContainer):
    """
    Some methods for Frame element
    """

    def add_element(self, EN: List[int], MI:List[int]=None, E:float=None, A:float=None, I:float=None, own_weight:float=0):
        """
        Adds an element to the system.
        :param EN: Element nodes
        :param MI: Matrix of indices
        :param E: Elastic modulus
        :param A: Cross-section area
        :param I: Moment of inertia
        :param own_weight: own weight of element
        :return: nothing, Void type
        """
        if MI is None:
            MI = []
        number = len(self)
        frame_el = ElementFrame(self.nodes_scheme, EN, MI, E, A, I, own_weight, number)
        # write element to container
        ElementContainer.add_element(self, element=frame_el)

        assert isinstance(MI, list), 'MI must be a list'
        assert isinstance(EN, list), 'EN must be a list'
        assert own_weight >= 0, "own_weight can't be negative"

    def set_stiffness(self, el_num, E=None, A=None, I=None, own_weight=None):
        """
        Sets the stiffness parameters of the element "el_num" from "ElementClass.elements_list"
        :param el_num: Number of the element, which we want to set stiffness
        :param E: Elastic modulus (Young's modulus). For all element types, except null-element
        :param A: cross-section area. For rod elements
        :param I: The moment of inertia. For rod elements
        :param own_weight: Own weight of frame element (linear meter)
        :return: None, void
        """
        if E is not None:
            self.elements_list[el_num].E = E
        if A is not None:
            self.elements_list[el_num].A = A
        if I is not None:
            self.elements_list[el_num].I = I
        if own_weight is not None:
            self.elements_list[el_num].own_weight = own_weight


class ElementFrame(ElementMethods):
    """
    Frame element class, contains list of frame elements in scheme
    """

    def __init__(self, nodes, EN:List[int], MI:List[int], E:float=None, A:float=None, I:float=None, own_weight:float=0,
                 number:int=None):
        """
        Creates an element
        :param EN: element nodes. There is 2 nodes for frame-element [node1, node2]
        :param MI: matrix of indices. [index1, index2, index3, index4, index5, index6]
        index - the number of degree of freedom
        :param E: elastic modulus (Young's ration). [N/m^2] or [Pa]
        :param A: cross sectional area of the element [m^2]
        :param I: moment of inertia [m^4]
        :param own_weight: own weight of the frame element [kg/m]
        """
        self.EN = EN
        self.MI = MI
        self.E = E
        self.A = A
        self.I = I
        self.own_weight = own_weight
        self.number = number

        l_c_s = self.add_parameters(self.EN, nodes)
        self.length, self.cosa, self.sina = l_c_s

    def add_parameters(self, EN, nodes):
        """
        Adds parameters: length of the element, cos(alpha) and sin(alpha),
        where alpha is angle between horizontal axis and positive direction of the element.
        The "positive direction" is chosen from 1st to 2nd node of the element
        :param EN: Element Nodes. [node1, node2] list
        :param nodes: existing nodes in scheme. Object of the class FEM.scheme.Node()
        :return: length, cos(a), sin(a) of the element
        """
        x1 = nodes[EN[0]].x  # horizontal coordinate of 1st node of the element
        x2 = nodes[EN[1]].x  # horizontal coordinate of 2nd node of the element
        y1 = nodes[EN[0]].y  # vertical coordinate of 1st node of the element
        y2 = nodes[EN[1]].y  # vertical coordinate of 2st node of the element
        dx = x1 - x2
        dy = y1 - y2
        length = math.sqrt(dx*dx + dy*dy)
        if length == 0:
            raise ValueError('Length of frame element #{} = 0 (equals zero).'
                             ' It wont work, it is needed to calculate cos(a) and sin(a)')

        return length, -dx/length, -dy/length

    def set_stiffness(self, E=None, A=None, I=None, own_weight=None):
        """
        Sets the stiffness parameters of the element "el_num" from "ElementClass.elements_list"
        :param E: Elastic modulus (Young's modulus). For all element types, except null-element
        :param A: cross-section area. For rod elements
        :param I: The moment of inertia. For rod elements
        :param own_weight: Own weight of frame element (linear meter)
        :return: None, void
        """
        if E is not None:
            self.E = E
        if A is not None:
            self.A = A
        if I is not None:
            self.I = I
        if own_weight is not None:
            self.own_weight = own_weight

    def form_ke_gsc(self):
        """
        Forms the stiffness matrix of the element in global system of coordinates
        :return: Stiffness matrix of the element in Global System of Coordinates
        """
        el_ke = self.__form_ke()
        el_lambda = self.__form_lambda__()
        el_lambda_t = el_lambda.transpose()
        return np.matmul(np.matmul(el_lambda_t, el_ke), el_lambda)

    def __form_lambda__(self):
        """
        Forms the matrix of direction cosines of frame element
        Frame element as object must have (cos, sin) the Methods2nodeElement.parameters_bar
        :return: matrix of direction cosines
        """
        cosa = self.cosa
        sina = self.sina
        return np.array([[cosa , sina, 0, 0    , 0   , 0],
                         [-sina, cosa, 0, 0    , 0   , 0],
                         [0    , 0   , 1, 0    , 0   , 0],
                         [0    , 0   , 0, cosa , sina, 0],
                         [0    , 0   , 0, -sina, cosa, 0],
                         [0    , 0   , 0, 0    , 0   , 1]])

    def __form_ke(self):
        """
        Forms the local stiffness matrix of the frame element
        :return: local stiffness frame matrix of the element el_num
        """
        EA = self.E * self.A
        EI = self.E * self.I
        L = self.length
        L2 = self.length**2
        L3 = self.length**3
        return np.array([[EA/L , 0        , 0       , -EA/L, 0        , 0       ],
                         [0    , 12*EI/L3 , 6*EI/L2 , 0    , -12*EI/L3, 6*EI/L2 ],
                         [0    , 6*EI/L2  , 4*EI/L  , 0    , -6*EI/L2 , 2*EI/L  ],
                         [-EA/L, 0        , 0       , EA/L , 0        , 0       ],
                         [0    , -12*EI/L3, -6*EI/L2, 0    , 12*EI/L3 , -6*EI/L2],
                         [0    , 6*EI/L2  , 2*EI/L  , 0    , -6*EI/L2 , 4*EI/L  ]])

    def mass_part_in_nodes(self):
        """
        Find how much mass in each node of the element
        :return: vector of masses parts in nodes of the frame element [node1_mass_part, node2_mass_part]
        """
        mass = (self.length/2)*self.own_weight
        return np.array([mass, mass])






