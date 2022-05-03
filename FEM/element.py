

class ElementContainer(object):
    """
    Class for adding any kind of element and count them, listing
    """
    def __init__(self, nodes_scheme):
        """
        add new object -> element
        """
        # nodes in scheme
        self.nodes_scheme = nodes_scheme
        # List of added elements
        self.elements_list = []
        # array of all element nodes for elements in container
        self.EN = []

    def __len__(self) -> int:
        return len(self.elements_list)

    def __str__(self) -> str:
        """
        What we are returning if an access to our object of class Element was as to string
        :return: Type of element and Amount of elements
        """
        return "Element type: {}, Amount {}".format(type(self), len(self))

    def __getitem__(self, element_number) -> int:
        """
        Element's list access through the indices at class object
        :param element_number: number of element witch we need to access
        :return: Object from elements_list[element_number]
        """
        return self.elements_list[element_number]

    def add_element(self, element):
        """
        Adds an element to the system. Here is general operation for any element type.
        This method is realized fully in specific element type container
        """
        self.elements_list.append(element)  # remember the element in list
        # adding element to the node. Node will have the info about which elements it has
        for node_number in element.EN:
            self.nodes_scheme[node_number].elements.append(element)
        # add EN of the element to greater container
        self.EN.append(element.EN)

    def clear(self):
        """
        Clears the data in the object
        :return: None
        """
        self.elements_amount = 0
        self.elements_list = []

    def nodes_coordinates(self, el_num):
        """
        Finds the coordinates of the nodes of the element
        :param el_num: number of the element (index of element it Container)
        :return: [[x1,x2...,xn],[y1,y2...,yn]] n - nodes amount of the element
        """
        x, y = [], []
        for node_number in range(len(self.elements_list[el_num].EN)):
            x.append(self.nodes_scheme[self.elements_list[el_num].EN[node_number]].x)
            y.append(self.nodes_scheme[self.elements_list[el_num].EN[node_number]].y)
        return x, y


class ElementMethods:
    """
    Class for creating new elements and do operations with them
    """
    EN = []
    MI = []

    def nodes_coordinates(self, nodes):
        """
        Finds the coordinates of the nodes of the element
        :param nodes: Object of class Node from FEM/Scheme. Object 4nodeElement must have added nodes
        :return: [[x1,x2...,xn],[y1,y2...,yn]] n - nodes amount of the element
        """
        x, y = [], []
        for node_number in self.EN:
            x.append(nodes[node_number].x)
            y.append(nodes[node_number].y)
        return x, y

    def nodes_coordinates_positive(self, nodes):
        """
        Finds the coordinates of the nodes of the element
        And then change them all to positive values with keeping ratios between them:
        |x2 - x1|, |y2 - y1| ect. are the same.
        It's like mooving all element to first quarter of Cartesian Coordinates system where x>0 and y>0
        :param nodes: Object of class Node from FEM/Scheme. Object 4nodeElement must have added nodes
        :return: [[x1,x2...,xn],[y1,y2...,yn]] n - nodes amount of the element
        """
        x, y = [], []
        for node_number in self.EN:
            x.append(nodes[node_number].x)
            y.append(nodes[node_number].y)
        min_x = min(x)
        min_y = min(y)
        if min_x < 0:
            for i in range(len(x)):
                x[i] -= min_x
        if min_y < 0:
            for j in range(len(y)):
                y[j] -= min_y
        return x, y

    def global_displacements(self, U):
        """
        Find global displacements on the element
        :param U: Global displacements of all scheme (result of np.linalg.solve())
        :return: list [] of global displacements on the element.
        Order like in MI of the element
        """
        displacements_on_element = [U[i] for i in self.MI]
        return displacements_on_element



