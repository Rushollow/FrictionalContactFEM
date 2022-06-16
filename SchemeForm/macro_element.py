# do the ME
from typing import List
from FEM.scheme import NodeContainer
from FEM.element import ElementContainer, ElementMethods
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_null import ElementNull
import weakref
import numpy as np
from input_data import ACCURACY_OF_STITCHING


def node_adjacent(sigma, eta):
    """
    Function with local coordinated on element
    :param sigma: local horizontal axis
    :param eta: local vertical axis
    :return: numpy.ndarray.__len__ = 4
    """
    return 1 / 4 * np.array([(1 - sigma) * (1 - eta),
                             (1 + sigma) * (1 - eta),
                             (1 + sigma) * (1 + eta),
                             (1 - sigma) * (1 + eta)], dtype=float)


class ElementMacroContainer(ElementContainer):
    """
    Container for macro elements. Used for adding ME to the scheme
    """
    instances = []  # instances of the class

    def __init__(self, nodes_scheme):
        """
        initializing container
        :param nodes_scheme: nodes in scheme. FEM.scheme.Node() instance
        element 4node container. FEM.element_4node.Element4NodeContainer() instance
        """
        ElementContainer.__init__(self, nodes_scheme=nodes_scheme)
        # 4 node elements in scheme. FEM.element_4node.Element4NodeContainer() instance
        self.element_4node = Element4NodeLinearContainer(nodes_scheme=nodes_scheme)
        # later all nodes used for ME will be added
        self.initial_nodes_numbers_of_macro_elements = []
        # used to calculate renumbering nodes for elements
        self.nodes_amount_in_scheme_before_fragmentation = len(nodes_scheme)
        # used to calculate future nodes' numbers for elements that are not fragmented
        self.future_nodes_in_scheme = 0
        # need this to add 4node elements fragment to scheme 4node elements
        self.element_4node_scheme = None

    def __getitem__(self, element_number):
        """
        What object of this class will return if we want to get access through the index
        Element's list access through the indices at class object
        :param element_number: number of element witch we need to get
        :return: object from elements_list[element_number]
        """
        return self.elements_list[element_number]

    def add_element(self, EN:List[int], frag_size:float=1, E:float=None, mu:float=None, t:float=None, own_weight:float=0,
                    frag_amount_v:int=None, frag_amount_h:int=None, frag_size_v:float=None, frag_size_h:float=None,
                    stitch:bool=False, stitch_list:List[int]=None):
        """
        Adds ME to the scheme.
        :param EN: element nodes. List [] of indices that refer to FEM.scheme.NodeContainer().elements object field
        First 2 nodes are bottom nodes, last 2 are top nodes.
        Nodes must be in order counter-clockwise
        :param frag_size: preferred size of one finite element after fragmentation
        assumed that vertical and horizontal sizes and equal
        :param E: Elastic modulus (Yong's modulus)
        :param mu: poisson's ratio
        :param t: thickness of the element
        :param own_weight: own_weight of the element (further calculated as nodes load, added to LoadVector)
        :param frag_amount_v: amount of 4node elements in vertical direction after fragmentation
        :param frag_amount_h: amount of 4node elements in horizontal direction after fragmentation
        :param frag_size_v: preferred vertical size
        :param frag_size_h: preferred horizontal size
        :param stitch: stitch this ME with others or not
        :param stitch_list: list of numbers of ME that should be stitched to this element.
        YOU MUST add only ME numbers that are already exists!
        :return:
        """
        if stitch_list is None:
            stitch_list = []
        new_me = ElementMacro(self.nodes_scheme, self.element_4node, EN, frag_size, E, mu, t, own_weight, frag_amount_v,
                              frag_amount_h, frag_size_v, frag_size_h, stitch, stitch_list=stitch_list)
        ElementContainer.add_element(self, new_me)

        assert frag_size > 0, 'frag_size must be positive value'

    def fragment_all(self, element_4node, element_frame, element_null):
        """
        Fragment all macro elements in scheme
        :return:
        """
        self._renumber_existing_elements_nodes_numbers(element_4node, element_frame, element_null)
        for me in self.elements_list:
            me.fragment()
        self._delete_coincident_nodes(self.elements_list, element_4node, element_frame, element_null)
        self._correct_nodes_indices()
        # adding fragmented 4node elements to scheme
        self.element_4node_scheme = element_4node
        self.element_4node_scheme.elements_list += self.element_4node.elements_list
        self.element_4node_scheme.EN += self.element_4node.EN

    def _renumber_existing_elements_nodes_numbers(self, element_4node, element_frame, element_null):
        """
        Renumbers existing elements' nodes to future numbers after fragmentation
        :param element_4node:
        :param element_frame:
        :param element_null:
        :return:
        """
        list_of_element_containers = [i for i in [element_4node, element_frame, element_null] if i is not None]
        self.future_nodes_in_scheme = len(self.nodes_scheme) - 1
        for me in self:
            n1 = self.future_nodes_in_scheme + 1
            n2 = n1 + me.elements_amount_h
            n4 = n1 + (me.elements_amount_h + 1) * me.elements_amount_v
            n3 = n4 + me.elements_amount_h
            future_me_element_nodes_numbers = [n1, n2, n3, n4]
            for element_container in list_of_element_containers:
                for element in element_container:
                    for i, node_number in enumerate(element.EN):
                        if node_number in me.EN:
                            index = me.EN.index(node_number)
                            element.EN[i] = future_me_element_nodes_numbers[index]
            self.future_nodes_in_scheme = n3

    def _delete_coincident_nodes(self, me_list, element_4node, element_frame, element_null):
        """
        Delete coincident nodes (nodes with same coodinates / overlapped)
        :param me_list: list of macro elements
        :param element_4node: 4 node elements in scheme
        :param element_frame: frame elements in scheme
        :param element_null: null-elements in scheme
        :return:
        """
        # remember which nodes were deleted (after that change EN in elements)
        deleted = np.zeros(len(self.nodes_scheme), dtype=int)
        # delete coincident nodes in me.EN and nodes added during fragmentation
        self._del_coincident_nodes_between_me_and_fragmented(me_list, deleted)
        # delete coincident nodes comparing nodes from 2 different ME
        self._del_coincident_nodes_between_two_me(me_list, deleted)
        #
        self._renumber_element_nodes_numbers(element_4node, element_frame, element_null, deleted)

    def _del_coincident_nodes_between_me_and_fragmented(self, me_list, deleted):
        """
        Delete coincident nodes between macro element and nodes that were added during fragmentation
        :param me_list: list of macro elements
        :param deleted: list with deleted nodes numbers,
        where deleted[i] equals the number of nodes that were deleted before i node
        :return:
        """
        # form list of initial_nodes_numbers_of_macro_elements
        for me in me_list:
            for node_number in me.EN:
                if node_number not in self.initial_nodes_numbers_of_macro_elements:
                    self.initial_nodes_numbers_of_macro_elements.append(node_number)
        # delete initial nodes of ME elements
        for node_number in self.initial_nodes_numbers_of_macro_elements:
            self.nodes_scheme.delete(node_number - deleted[node_number])  # delete node
            for k in range(node_number, len(deleted)):
                deleted[k] += 1

    def _del_coincident_nodes_between_two_me(self, me_list, deleted):
        """
        delete coincident nodes between two macro element nodes that were added during fragmentation
        :param me_list:
        :param deleted:list with deleted nodes numbers,
        where deleted[i] equals the number of nodes that were deleted before i node
        :return:
        """
        me_checked = 0
        checked_nodes = []
        for me1_num in range(len(me_list)):  # iterate over one ME
            me_checked += 1
            for me2_num in range(me_checked, len(me_list)):  # iterate over another ME
                # next value is int 0 or 1 or 2. It is 2 if both True to stitch and 0 if both is False to stitch
                # and 1 if at least 1 is True to stitch
                to_stitch = me_list[me1_num].stitch + me_list[me2_num].stitch
                if to_stitch == 0:  # if both False
                    continue  # do not do the stitch
                if to_stitch == 1:  # if only 1 of them is True to stitch
                    stitch_list = me_list[me2_num].stitch_list + me_list[me1_num].stitch_list
                    if me1_num not in stitch_list:  # if macro element number is not in stitch list
                        continue  # do not do the stitch
                for i in me_list[me1_num].nodes_numbers_me:  # iterate over nodes of one ME  # TODO: HERE! STICHING
                    node_me1 = self.nodes_scheme[i - deleted[i]]
                    for j in me_list[me2_num].nodes_numbers_me:  # iterate over nodes of another ME
                        node_me2 = self.nodes_scheme[j - deleted[j]]
                        if node_me1 is not node_me2:
                            # if nodes coincident
                            if np.absolute(node_me1.x - node_me2.x) < ACCURACY_OF_STITCHING and \
                                    np.absolute(node_me1.y - node_me2.y) < ACCURACY_OF_STITCHING:
                                self.nodes_scheme.delete(j - deleted[j])  # delete node in 2nd ME
                                for k in range(j, len(deleted)):
                                    if k not in checked_nodes:
                                        deleted[k] += 1
                                deleted[j] = j - (i - deleted[i])  # take the node number from first ME
                                checked_nodes.append(j)
                                break


    # def _node_have_null_element(self, node):
    #     """
    #     Check if node have null element
    #     In this case coinsident nodes could exist in scheme in case of contact boundary. (gap = 0)
    #     :param node: node in scheme. FEM.scheme.Node
    #     :return:
    #     """
    #     for element in node.elements:
    #         if isinstance(element, ElementNull):
    #             return True
    #     return False

    def _renumber_element_nodes_numbers(self, element_4node, element_frame, element_null, deleted):
        """
        Renumber element's nodes numbers in element.EN vectors. Because of the deletion of coincident nodes
        :param element_4node: 4node elements in scheme
        :param element_frame: frame elements in scheme
        :param element_null: null-elements in scheme
        :param deleted: list of numbers of deleted nodes from scheme
        :return:
        """
        # amount_of_nodes_used_for_me = len(self.initial_nodes_numbers_of_macro_elements)
        # amount_of_nodes_without_me = self.nodes_amount_in_scheme_before_fragmentation - amount_of_nodes_used_for_me
        # for i in self.initial_nodes_numbers_of_macro_elements:
        #     deleted[i] -= amount_of_nodes_without_me
        list_of_element_containers = [i for i in [self.element_4node, element_4node, element_frame, element_null]
                                      if i is not None]
        for element_container in list_of_element_containers:
            for element in element_container:
                for i, node_number in enumerate(element.EN):
                    element.EN[i] -= deleted[node_number]

    def _correct_nodes_indices(self):
        """
        correct indices and numbers of the nodes (case some nodes where deleted and numeration is different)
        :return:
        """
        for i, node in enumerate(self.nodes_scheme):
            node.indices = [2 * i, 2 * i + 1]
            node.number = i


class ElementMacro(ElementMethods):
    """
    Element for fragmentation. It is 4 node element.
    First 2 nodes are bottom nodes, last 2 are top nodes
    'bottom' and 'top' are need to understand variables: vertical and horizontal fragmentation

    MacroElement:
                      'top'
                3 ------------- 2 - node number             |---------------|
                |               |                        |\ |       |       |
         'left' |               | 'right'          size_v|| |       | size_h|
                |               |                        \| |-------|<----->|
                |               |                           |       |       |
                0 ------------- 1                           |       |       |
                    'bottom'                                -----------------
    """
    def __init__(self, nodes, element_4node, EN, frag_size, E, mu, t, own_weight,
                 frag_amount_v, frag_amount_h, frag_size_v, frag_size_h, stitch, stitch_list):
        """
        Initialize macro element (ME)
        :param nodes: Nodes in scheme
        :param element_4node: 4 node elements in scheme
        :param EN: element nodes numbers of the ME. List [] of numbers.
        :param frag_size: preferred size of one finite element after fragmentation
        :param E: Elastic modulus
        :param mu: Poisson's ratio
        :param t: thickness of the element
        :param own_weight: own weight of the element
        :param frag_amount_v: amount of 4node elements in vertical direction after fragmentation
        :param frag_amount_h: amount of 4node elements in horizontal direction after fragmentation
        :param frag_size_v: preferred vertical size
        :param frag_size_h: preferred horizontal size
        :param stitch: stitch this ME with others or not
        :param stitch_list: list of numbers of the ME that will be stitched to this ME
        """
        self.nodes_scheme = nodes
        self.element_4node = element_4node
        self.EN = EN
        self.E = E
        self.mu = mu
        self.t = t
        self.own_weight = own_weight
        self.frag_size = frag_size
        self.frag_size_v = frag_size_v
        self.frag_size_h = frag_size_h
        self.stitch = stitch
        self.stitch_list = stitch_list
        # list of elements into which ME was fragmented.
        self.elements_list = []
        # nodes numbers im macro element
        self.nodes_numbers_me = []
        # set a value if vertical or horizontal fragmentation size is None
        if self.frag_size_v is None: self.frag_size_v = self.frag_size
        if self.frag_size_h is None: self.frag_size_h = self.frag_size
        # calculate sizes of element_4node
        self.elements_amount_v = None
        self.elements_amount_h = None
        self._find_number_of_elements()
        if frag_amount_v is not None:
            self.elements_amount_v = frag_amount_v
        if frag_amount_h is not None:
            self.elements_amount_h = frag_amount_h

    def _find_number_of_elements(self):
        """
        Find amount of elements that should be in horizontal and vertical direction
        :return:
        """
        # average horizontal size of the ME
        h_bot = self._distance_between_2nodes(self.nodes_scheme[self.EN[0]], self.nodes_scheme[self.EN[1]])
        h_top = self._distance_between_2nodes(self.nodes_scheme[self.EN[2]], self.nodes_scheme[self.EN[3]])
        self.h_avg_size = (h_bot + h_top) / 2
        # average vertical size of the ME
        v_left = self._distance_between_2nodes(self.nodes_scheme[self.EN[0]], self.nodes_scheme[self.EN[3]])
        v_right = self._distance_between_2nodes(self.nodes_scheme[self.EN[1]], self.nodes_scheme[self.EN[2]])
        self.v_avg_size = (v_left + v_right) / 2
        if self.elements_amount_h is None:
            self.elements_amount_h = int(np.around(self.h_avg_size / self.frag_size_h))
        if self.elements_amount_v is None:
            self.elements_amount_v = int(np.around(self.v_avg_size / self.frag_size_v))
        if self.elements_amount_h < 1: self.elements_amount_h = 1
        if self.elements_amount_v < 1: self.elements_amount_v = 1

    def _distance_between_2nodes(self, node1, node2):
        """
        Calculates the distance between two nodes in scheme.
        :param node1: FEM.scheme.Node() object
        :param node2: FEM.scheme.Node() object
        :return: euclidean distance
        """
        x1 = node1.x  # horizontal coordinate of 1st node of the element
        x2 = node2.x  # horizontal coordinate of 2nd node of the element
        y1 = node1.y  # vertical coordinate of 1st node of the element
        y2 = node2.y  # vertical coordinate of 2st node of the element
        a = np.array([x1 , y1])
        b = np.array([x2, y2])
        return np.linalg.norm(a-b)

    def fragment(self):
        self._add_nodes()
        self._add_elements()

    def _add_nodes(self):
        """
        Creating net of nodes on a square 2x2 meters
        :return:
        """
        node_number = len(self.nodes_scheme)
        # add nodes for fragmentation
        for i in range(self.elements_amount_v+1):
            for j in range(self.elements_amount_h+1):
                x = j * 2 / self.elements_amount_h
                y = i * 2 / self.elements_amount_v
                self.nodes_scheme.add_node(x, y)
                # add nodes to array to know which to consider when fragment ME
                self.nodes_numbers_me.append(node_number)
                node_number += 1
        # initial nodes coordinates of ME
        nodes_coordinates_me = self.nodes_coordinates(self.nodes_scheme)
        self._adjust_nodes_to_positions(nodes_coordinates_me)

    def _adjust_nodes_to_positions(self, nodes_coordinates_me):
        x, y = nodes_coordinates_me
        x, y = np.array(x), np.array(y)
        for i in self.nodes_numbers_me:
            node = self.nodes_scheme[i]
            sigma = node.x - 1
            eta = node.y - 1
            adjusted_locally = node_adjacent(sigma, eta)  # adjust on local coordinates on 2x2 square
            # adjust to real size
            node.x = np.sum(x.T * adjusted_locally)
            node.y = np.sum(y.T * adjusted_locally)

    def _add_elements(self):
        """
        Add 4 node elements on created nodes
        :return:
        """
        # cycle through future elements numbers. i - number of adding element
        for i in range(self.elements_amount_v * self.elements_amount_h):
            # calculate nodes numbers (4 is amount of nodes in ME)
            nodes_added_fragment = (self.elements_amount_v + 1) * (self.elements_amount_h + 1)
            k1 = i + i//self.elements_amount_h + len(self.nodes_scheme) - nodes_added_fragment
            k2 = k1 + 1
            k3 = self.elements_amount_h + k1 + 2
            k4 = k3 - 1
            self.element_4node.add_element(EN=[k1, k2, k3, k4],
                                           E=self.E,
                                           mu=self.mu,
                                           t=self.t,
                                           own_weight=self.own_weight)


