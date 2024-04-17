from FEM.element import ElementContainer, ElementMethods
from input_data import PLANE_STRAIN
import numpy as np
from typing import List


class Element4NodeLinearContainer(ElementContainer):
    """
    Class for storing and adding 4node elements
    """

    # Some methods for 4node elements
    def add_element(self, EN: List[int], MI: List[int] = None,
                    E: float = None, mu: float = None, t: float = None, own_weight: float = 0):
        """
        Adding element 4NodeLinear to the scheme
        :param EN: element nodes: [node1, node2, node3, node4]. Numbers of the element nodes
        :param MI: matrix of indices [index1, index2, index3, index4, index5, index6, index7, index8]
        :param E: Elastic modulus (Young's modulus)
        :param mu: Poisson's ratio (for plate elements)
        :param t: thickness of the element (for plate elements)
        :param own_weight: Own weight of the element. (kilograms/meter^2) [kg/m^2] - Si
        :return: None, just add element
        """
        assert isinstance(EN, list), 'EN must be list'
        assert E > 0, 'E (elastic modulus must be positive'
        assert (mu > 0 and mu <= 0.5), "mu (Poisson's ratio) must be 0<mu<=0.5"
        assert t > 0, 't (thickness) must be positive'
        assert own_weight >= 0, 'own weight cannot be negative'

        if MI is None:
            MI = []
            # for node_num in EN:
            #     MI += self.nodes_scheme[node_num].indices[:2]
        element = Element4NodeLinear(EN, MI, E, mu, t, own_weight)
        # write instance to the list and remember other parameters
        ElementContainer.add_element(self, element)

    def set_stiffness(self, el_num, E=None, mu=None, t=None, own_weight=0):
        """
        Sets the stiffness parameters of the element "el_num" from "ElementClass.elements_list"
        :param el_num: Number of the element, which we want to set stiffness
        :param E: Elastic modulus (Young's modulus). For all element types, except null-element
        :param mu: Poisson's ratio. For plate elements
        :param t: Thickness of the element. For plate elements
        :param own_weight: Own weight of 4node element (square meter)
        :return: None, void
        """
        if E is not None:
            self.elements_list[el_num].E = E
        if mu is not None:
            self.elements_list[el_num].mu = mu
        if t is not None:
            self.elements_list[el_num].t = t
        if own_weight is not None:
            self.elements_list[el_num].own_weight = own_weight

    def correct_indices(self):
        """
        if operations with nodes have been done - it is needed to correct MI of the element
        :return:
        """
        for element in self.elements_list:
            MI = []
            for node_num in element.EN:
                MI += self.nodes_scheme[node_num].indices[:2]
            element.MI = MI

    def find_elements_along_segment(self, nodes_container, point1, point2, rtol=1e-1) -> (list, list[int]):
        """
        Finds all 4node elements between 2 given point
        calculates center on the element, if it is near the segment -> OK!
        :param nodes_container:
        :param point1: first point (x1, y1)
        :param point2: second point (x2, y2)
        :param rtol: relative tolerance
        :return:
        """
        central_points = []  # central points of elements 4node (average coordinates x and y)
        elements = []
        list_of_numbers = []
        for element in self.elements_list:
            avg_x, avg_y = 0, 0  # average x coordinate
            for nn in element.EN:
                avg_x += nodes_container[nn].x
                avg_y += nodes_container[nn].y
            avg_x /= 4  # case its 4 nodes in this element
            avg_y /= 4
            central_points.append((avg_x, avg_y))

        part_of_cp1 = point2[0] - point1[0]  # part of cross product
        part_of_cp2 = point2[1] - point1[1]
        for i, cp in enumerate(central_points):
            # The absolute value of the cross product is twice the area of the triangle formed by the three points
            cross_product = (cp[1] - point1[1]) * part_of_cp1 - (cp[0] - point1[0]) * part_of_cp2
            if np.isclose(cross_product, 0, rtol=rtol) and\
                    min(point1[0]-rtol, point2[0]-rtol) <= cp[0] <= max(point1[0]+rtol, point2[0]+rtol) and\
                    min(point1[1]-rtol, point2[1]-rtol) <= cp[1] <= max(point1[1]+rtol, point2[1]+rtol):
                elements.append(self.elements_list[i])
                list_of_numbers.append(i)
        return elements, list_of_numbers


class Element4NodeLinear(ElementMethods):
    # Class of 4node linear element general form with 2 degrees of freedom in each node.
    # It can't be degenerate (have intersection of opposite sides or have inner corner above 180 degree)

    def __init__(self, EN, MI, E=None, mu=None, t=None, own_weight=0.0):
        """
        Creates an element
        :param EN: element nodes: [node1, node2, node3, node4]. Numbers of the element nodes
        :param MI: matrix of indices [index1, index2, index3, index4, index5, index6, index7, index8]
        :param E: Elastic modulus (Young's modulus)
        :param mu: Poisson's ratio (for plate elements)
        :param t: thickness of the element (for plate elements)
        :param own_weight: Own weight of the element. (kilograms/meter^2) [kg/m^2] - Si
        """
        self.MI = MI
        self.EN = EN
        self.E = E
        self.mu = mu
        self.t = t
        if own_weight > 0:
            own_weight = -own_weight
        self.own_weight = own_weight * t
        # stresses calculated by def stress
        self.sigma_x: float = 0
        self.sigma_y: float = 0

    def elastic_matrix(self, plane_deformation=PLANE_STRAIN):
        """
        Elastic matrix D4 of the element.
        :param self: Object of 4node element
        :param plane_deformation: Checks if it is "plane deformation state" or "plane stress state"
        :return: elastic matrix
        """
        E = self.E
        mu = self.mu
        if plane_deformation:  # modify Elastic modulus (E) and Poisson's ratio (mu)
            E = E / (1 - mu ** 2)
            mu = mu / (1 - mu)
        # else it's plane stress
        return E / (1 - mu ** 2) * np.array([[1, mu, 0],
                                             [mu, 1, 0],
                                             [0, 0, (1 - mu) / 2]])

    def __det_jacobian_matrix__(self, nodes, k, n):
        """
        Calculates the determinant of Jacobian matrix for 4node element
        :param self: Object of 4node element
        :param nodes: Object of Scheme.Node class with nodes in nodes_list
        :param k: ksi variable to know where we should calculate
        :param n: nu variable to know where we should calculate
        :return: determinant of Jacobian matrix for 4node element
        """
        x, y = self.nodes_coordinates_positive(nodes)
        res = (x[0] * y[1] - x[1] * y[0] - x[0] * y[3] + x[1] * y[2] - x[2] * y[1] + x[3] * y[0] + x[2] * y[3] - x[3] *
               y[2] - n * x[0] * y[1] + n * x[1] * y[0] + n * x[0] * y[2] - n * x[2] * y[0] - n * x[1] * y[3] + n *
               x[3] * y[1] + n * x[2] * y[3] - n * x[3] * y[2] - k * x[0] * y[2] + k * x[2] * y[0] + k * x[0] * y[3] +
               k * x[1] * y[2] - k * x[2] * y[1] - k * x[3] * y[0] - k * x[1] * y[3] + k * x[3] * y[1]) / 8
        return np.abs(res)

    @staticmethod
    def __ksi_nu__(point_num):
        """
        Finds the local coordinates ksi and nu of the 4node element
        :param point_num: number of the point where we need to do numerical integration
        :return:
        """
        loc = 0.577350269
        res = []
        if point_num == 0:
            res.append(-loc);
            res.append(-loc)
        elif point_num == 1:
            res.append(-loc);
            res.append(loc)
        elif point_num == 2:
            res.append(loc);
            res.append(-loc)
        elif point_num == 3:
            res.append(loc);
            res.append(loc)
        return res

    def form_ke(self, nodes):
        """
        Forms the local stiffness matrix of the 4node element
        :param nodes: Object of the Scheme.Node class. Expect to have nodes in "list_nodes"
        :return: local stiffness matrix of the element
        """
        summ = 0
        # cycle for numerical integration of Bt*D*B*detJ equation
        for point in range(4):  # it's 4 point of numerical integration
            [ksi, nu] = self.__ksi_nu__(point)
            det_j = self.__det_jacobian_matrix__(nodes, ksi, nu)
            b = self.__gradient_matrix__(nodes, ksi, nu)
            b_t = b.transpose()
            d = self.elastic_matrix()
            # sum += b_t*d*b*det_j
            summ += (np.matmul(np.matmul(b_t, d), b)) * det_j
        # t is thickness of the element
        return self.t * summ

    # Gradient matrix for 4-node linear element (general)
    def __gradient_matrix__(self, nodes, k, n):
        """
        Gradient matrix B4. It is for local coordinates ksi = k , and nu = n
        :param nodes: The object of Scheme.Node class with nodes in nodes_list
        :param k: Local horizontal coordinates on element
        :param n: Local vertical coordinates on element
        :return: Gradient matrix. float numpy 2D array
        """
        x, y = self.nodes_coordinates_positive(nodes)

        # denominator

        den = self.__det_jacobian_matrix__(nodes, k, n) * 8

        n11 = -(n * (y[1] - y[2]) + k * (y[2] - y[3]) - y[1] + y[3]) / (den)
        n12 = (n * (x[1] - x[2]) + k * (x[2] - x[3]) - x[1] + x[3]) / (den)
        n21 = (n * (y[0] - y[3]) + k * (y[2] - y[3]) - y[0] + y[2]) / (den)
        n22 = -(n * (x[0] - x[3]) + k * (x[2] - x[3]) - x[0] + x[2]) / (den)
        n31 = -(n * (y[0] - y[3]) + k * (-y[0] + y[1]) + y[1] - y[3]) / (den)
        n32 = (n * (x[0] - x[3]) + k * (-x[0] + x[1]) + x[1] - x[3]) / (den)
        n41 = (n * (y[1] - y[2]) + k * (-y[0] + y[1]) + y[0] - y[2]) / (den)
        n42 = -(n * (x[1] - x[2]) + k * (-x[0] + x[1]) + x[0] - x[2]) / (den)

        return np.array([[n11, 0, n21, 0, n31, 0, n41, 0],
                         [0, n12, 0, n22, 0, n32, 0, n42],
                         [n12, n11, n22, n21, n32, n31, n42, n41]])

    def mass_part_in_nodes(self, nodes):
        """
        Calculate masses share in 4node element as concentrated load in nodes
        :param nodes: Object of Scheme.Node class. Must have nodes in nodes_list
        Use nodes_coordinates method on element to find them
        :return: Masses as concentrated load [Fn1, Fn2, Fn3, Fn4], where Fni - Load value of i node of the element
        """

        x, y = self.nodes_coordinates(nodes)

        a1 = (2 * x[0] * (y[1] - y[3]) + x[1] * (y[2] - 2 * y[0] + y[3]) + x[2] * (y[3] - y[1]) + x[3] * (
                    2 * y[0] - y[1] - y[2])) / 12
        a2 = (x[0] * (2 * y[1] - y[2] - y[3]) + 2 * x[1] * (y[2] - y[0]) + x[2] * (y[0] - 2 * y[1] + y[3]) + x[3] * (
                    y[0] - y[2])) / 12
        a3 = (x[0] * (y[1] - y[3]) + x[1] * (2 * y[2] - y[0] - y[3]) + 2 * x[2] * (y[3] - y[1]) + x[3] * (
                    y[0] + y[1] - 2 * y[2])) / 12
        a4 = (x[0] * (y[1] + y[2] - 2 * y[3]) + x[1] * (-y[0] + y[2]) + x[2] * (-y[0] - y[1] + 2 * y[3]) + x[3] * (
                    +2 * y[0] - 2 * y[2])) / 12
        return np.abs(np.array([a1, a2, a3, a4]))

    def set_stiffness(self, E=None, mu=None, t=None, own_weight=0):
        """
        Sets the stiffness parameters of the element "el_num" from "ElementClass.elements_list"
        :param E: Elastic modulus (Young's modulus). For all element types, except null-element
        :param mu: Poisson's ratio. For plate elements
        :param t: Thickness of the element. For plate elements
        :param own_weight: Own weight of 4node element (square meter)
        :return: None, void
        """
        if E is not None:
            self.E = E
        if mu is not None:
            self.mu = mu
        if t is not None:
            self.t = t
        if own_weight is not None:
            self.own_weight = own_weight

    def __str__(self):
        return f'4 node element: {self.EN}'

    def _element_displacements(self, U: [float]):
        element_dis = []
        for indx in self.MI:
            element_dis.append(U[indx])
        return element_dis

    def _bilinear_shape_function(self, x, y):
        """
        Compute the bilinear shape function derivative
        at a given natural coordinate (xi, eta).
        """
        dN = np.array([[-(1 - y), (1 - y), (1 + y), -(1 + y)],
                      [-(1 - x), -(1 + x), (1 + x), (1 - x)]]) / 4
        return dN

    def stress(self, U:[float], nodes):
        """
        Calculate stresses (sigma_x, sigma_y)
        :param U: vector of global displacements
        :param nodes: scheme nodes container
        :return: None
        """
        el_dis = np.array(self._element_displacements(U))
        x_coord, y_coord = self.nodes_coordinates_positive(nodes=nodes)

        sigma_x, sigma_y = 0, 0

        for i in range(4):
            x, y = x_coord[i], y_coord[i]
            # Compute bilinear shape functions and their derivatives
            dN = self._bilinear_shape_function(x, y)
            # Compute Jacobian matrix
            J = np.dot(dN, el_dis.reshape(-1, 2))
            # Compute strains
            epsilon_xi = np.dot(np.linalg.inv(J), dN) @ el_dis.reshape(-1, 2)
            epsilon_yi = np.dot(np.linalg.inv(J), dN) @ el_dis.reshape(-1, 2)
            # Compute stresses
            sigma_x += self.E / (1 - self.mu ** 2) * (epsilon_xi + self.mu * epsilon_yi)
            sigma_y += self.E / (1 - self.mu ** 2) * (epsilon_yi + self.mu * epsilon_xi)

        return np.mean(sigma_x), np.mean(sigma_y)
