import numpy as np
import math
from FEM.element import ElementContainer, ElementMethods
from FEM.scheme import find_max_degree_of_freedom
from input_data import ACCURACY_OF_LCP


class ElementNullContainer(ElementContainer):
    """
    Class for adding null-elements
    """
    def __init__(self, nodes_scheme):
        """
        Initialize new object -> element
        :param nodes_scheme: nodes in scheme. FEM.scheme.ElementNodeContainer() instance
        """
        ElementContainer.__init__(self, nodes_scheme)
        # list of contact pairs for comfortable iteration. List have the same elements as self.elements_list
        self.contact_pairs_amount = 0
        # but group by contact pairs
        self.contact_pairs_list = []

    def add_element(self, EN, cke=None, alpha=None, cosa=None, sina=None, gap_length=0, add_t_el=True, orientation='n'):
        """
        Adds an element to the system.
        :param EN: [1st node number, 2nd node number] list of 2 nodes numbers (contact pair) for null element
        :param cke: stiffness of the null-element (it must be average over Global Stiffness Matrix elements or near)
        :param alpha: angle from the horizontal axis counter-clockwise to the direction of null-element
        Assumed that horizontal axis directed to the right
        :param cosa: cosine of alpha
        :param sina: sinus of alpha
        :param gap_length: length of the gap from 1 node to second of the contact pair
        :param add_t_el: If true, adding a null-element tangential to the contact zone in the same nodes
        Use only if you are adding null-element 'normal' to the contact zone.
        :param orientation: orientation of the null-element to the contact zone.
        :return: None
        """
        MI = []
        if gap_length == 0:  # check if gap length could be not zero
            dx = self.nodes_scheme[EN[1]].x - self.nodes_scheme[EN[0]].x
            dy = self.nodes_scheme[EN[1]].y - self.nodes_scheme[EN[0]].y
            distance = math.sqrt(dx * dx + dy * dy)  # distance between nodes
            if distance > ACCURACY_OF_LCP:
                alpha2 = math.acos(dx / distance)
                gap_length = math.fabs(math.cos(math.pi - alpha + alpha2)) * distance
                print('Added gap_length ={} to null-element {}'.format(gap_length, len(self) + 1))
        null_el_n = ElementNull(EN, MI, self.contact_pairs_amount, cke, alpha, cosa, sina, gap_length,
                                orientation=orientation)
        ElementContainer.add_element(self, element=null_el_n)
        null_el_t = None
        if add_t_el is True:
            MI_t = []
            null_el_t = ElementNull(EN, MI_t, self.contact_pairs_amount, cke, alpha - math.pi / 2, cosa, sina,
                                    gap_length=0, orientation='t')
            ElementContainer.add_element(self, element=null_el_t)
        # contact pairs
        contact_pair = ContactPair(null_el_n, null_el_t)
        self.contact_pairs_list.append(contact_pair)
        self.contact_pairs_amount += 1

        assert isinstance(EN, list), 'EN must be a list!'

    def _order_contact_pairs(self):
        """
        Sorts contact pairs: first - frictional (with normal and tangent null-elements)
                            last - frictionless (with normal only null elements)
        :return: None, changes initial list of contact pairs
        """
        pairs_frictional = []
        pairs_frictionless = []
        for contact_pair in self.contact_pairs_list:
            if contact_pair.type == 'frictional':
                pairs_frictional.append(contact_pair)
            else:
                pairs_frictionless.append(contact_pair)
        self.contact_pairs_list = pairs_frictional + pairs_frictionless

    def order_elements_null(self):
        """
        Orders all null elements: first comes 'normal' null-elements than 'tangent' to the contact zone
        If we had 2 arrays with 'normal' and 'tangent' null-elements:
        order is that some k 'normal' null-element must be in the same contact pair as k null-element for tangent.
        HOW TO understand text BELOW: nulli_jk: null - null-element, i - number of the null-element,
                        j - number of the contact pair, k - orientation (could be 'n' or 't' - normal of tangent)
        initial order = [null1_1n, null2_1t, null3_2n, null4_2t,...nulln_mnt]
        or if we have 'normal' more than 'tangent' it could be:
                      = [null1_1n, null2_1t, null3_2n, null_4_3n,...nulln_mnt]
        we need sorted = [null1_n, null3_n... nulln_mn, null2_1t, null4_2t,..., nulln_mt]

        It means that we need to sort in 2 steps:
        1) orientation 'n' - first, than 't'
        2) 'n' oriented elements with 't' in the same contact pair should be first
        :return: None
        """
        # sort by orientation
        elements_null_n = []
        elements_null_t = []
        # sort by availability(presence) of 't' null-element in the contact pair with 'n' null-element
        self._order_contact_pairs()  # first sort contact pairs by frictional and frictionless
        for contact_pair in self.contact_pairs_list:
            elements_null_n.append(contact_pair.null_element_n)
            if contact_pair.null_element_t is not None:
                elements_null_t.append(contact_pair.null_element_t)
        self.elements_list = elements_null_n + elements_null_t

    def correct_indices(self):
        """
        if operations with nodes have been done - it is needed to correct MI of the element
        :return:
        """
        # first correct horizontal and vertical indices from nodes
        for element in self.elements_list:
            MI = []
            for node_num in element.EN:
                MI += self.nodes_scheme[node_num].indices[:2]
            element.MI = MI
        max_dof = find_max_degree_of_freedom()
        for element in self.elements_list:
            max_dof += 1
            element.MI.append(max_dof)


class ContactPair:
    """
    Contact pair. Each contact pair is 2 nodes connected by 1 or 2 null-elements
    """
    def __init__(self, null_element_n, null_element_t):
        self.null_element_n = null_element_n
        self.null_element_t = null_element_t
        self.type = 'frictional'
        if null_element_t is None:
            self.type = 'frictionless'


class ElementNull(ElementMethods):
    """
    Null element class, 5 degrees of freedom
    """
    def __init__(self, EN, MI, contact_pair, cke, alpha=None, cosa=None, sina=None, gap_length=0, orientation=None):
        """
        Initialize null-element
        :param EN: element nodes.
        :param MI: matrix of indices
        :param contact_pair: number of the contact pair in which null-element is
        :param cke: Longitudinal stiffness for null-element
        We should choose this stiffness to avoid a degenerate global stiffness matrix
        :param alpha: Angle between horizontal axis and positive direction of the null-element.
        :param cosa: cosine of the angle
        :param sina: sinus of the angle
        :param gap_length: value of the length between nodes of contact bodies
        :param orientation: orientation of the null-element.
        It could be 'n' - normal to the contact zone and 't' - tangent to the contact zone
        """
        self.EN = EN
        self.MI = MI
        self.cke = cke
        self.gap_length = gap_length
        self.alpha = alpha
        self.orientation = orientation
        self.contact_pair = contact_pair
        # If we set alpha, calc cos and sin
        if alpha is not None:
            if alpha == math.pi/2:
                self.cosa = 0
                self.sina = 1
            elif alpha == 0:
                self.cosa = 1
                self.sina = 0
            else:
                self.cosa = math.cos(alpha)
                self.sina = math.sin(alpha)
        # If we set cos and sin when create null-element, use it
        if cosa is not None: self.cosa = cosa
        if sina is not None: self.sina = sina

    def form_ke_gsc(self):
        """
        Forms the stiffness matrix of the element in global system of coordinates
        :return: Stiffness matrix of the element in Global System of Coordinates
        """
        el_ke = self.form_ke()
        el_lambda = self.form_lambda()
        el_lambda_t = el_lambda.transpose()
        return np.matmul(np.matmul(el_lambda_t, el_ke), el_lambda)

    def form_ke(self):
        """
        Stiffness matrix for null-element in local system of coordinates
        :return: local stiffness matrix
        """
        cke = self.cke
        return np.array([[ cke,    0, -cke],
                         [   0, -cke,  cke],
                         [-cke,  cke,    0]])

    def form_qe(self, displacement=1):
        """
        Load vector for null-element in local system of coordinates
        :param displacement: the displacement for contact pair which is connected by this null-element [m]
        :return: local load vector
        """
        cke = self.cke * displacement
        return np.array([[-cke],
                         [   0],
                         [ cke]])

    def form_lambda(self):
        """
        Forms the matrix of direction cosines of null-element
        Null-element as object must have (cosa, sina). CreateNullElement.add_parameters_null_el(*) can add them
        :return: matrix of direction cosines
        """
        cosa = self.cosa
        sina = self.sina
        return np.array([[cosa, sina, 0   , 0   , 0],
                         [0   , 0   , cosa, sina, 0],
                         [0   , 0   , 0   , 0   , 1]])

    def form_Qe(self, displacement=1):
        """
        load vector for null-element in global system of coordinates
        :param displacement: the displacement for contact pair which is connected by this null-element [m]
        :return: global element load vector
        """
        lambda_el = self.form_lambda()
        qe = self.form_qe(displacement)
        return np.matmul(np.transpose(lambda_el), qe)

    def get_strain_effort(self, u_vector):
        """
        Get the strain effort in null-element
        :param u_vector: global displacements vector of the system
        :return: float number, strain effort
        """
        ke2nd_row = self.form_ke()[1]
        lambda_null_el = self.form_lambda()
        u_global_on_element = [u_vector[index] for index in self.MI]    # Find global displacements on element
        return -np.matmul(ke2nd_row, np.matmul(lambda_null_el, u_global_on_element))
