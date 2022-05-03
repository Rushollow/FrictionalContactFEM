import matplotlib as mpl
from matplotlib import markers
from FEM.element_null import ElementNull
from FEM.element_4node import Element4NodeLinear
from SchemeForm.macro_element import ElementMacro
from FEM.scheme import list_of_all_elements, NodeContainer, form_u_contact
from input_data import SCALE_DEF, FRICTION_COEFFICIENT
import math
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
    return  v1 * _frame_n1(x, L) + fi1 * _frame_n2(x, L) + \
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


class PlotScheme:
    """
    Class for plotting. methods for schemes, results.
    """
    def __init__(self, nodes, element_null, sm, lv, u_linear, lemke,
                 element_container_obj, element_frame,
                 partition=10, scale_def=SCALE_DEF, text=True):
        """
        Working with solution information, visualizing it. Uses matplotlib
        :param nodes: nodes in scheme. FEM.scheme.NodeContainer() class obj
        :param element_null: null-elements in scheme. FEM.element_null.ElementNullContainer() obj
        :param sm: filled stiffness matrix. FEM.scheme.StiffnessMatrix() obj
        :param lv: filled load vector. FEM.scheme.LoadVector() obj
        :param u_linear: vector of global displacements for linear formulation.
        :param lemke: solved contact problem. LCP.lemke.Lemke() obj
        :param element_container_obj: any type of the elements in scheme which need to plot.
        Ex.: FEM.element_4node.Element4NodeContainer() obj
        :param element_frame: frame elements in scheme. FEM.element_frame.ElementsFrameContainer() obj
        :param partition: division by n parts the frame elements to get smother curve while deformed.
        :param scale_def: scale of global displacements.
        :param text: bool value if text values should be visualized on plot
        """
        # data to use
        self.nodes = nodes
        self.element_null = element_null
        self.sm = sm
        self.lv = lv
        self.u_linear = u_linear
        self.lemke = lemke
        self.element_container_obj = element_container_obj  # we can choose plot another elements or not
        self.element_frame = element_frame  # we can choose plot frame elements or not
        self.partition = partition
        self.scale = scale_def
        self.u_contact_anim = []  # all vectors with global displacements during Lemke algorithm
        if lemke is not None:
            self.u_contact_anim = form_u_contact(self.sm, self.lv,
                                                 self.lemke.zn_anim, self.lemke.zt_anim, self.element_null)
        # parameters
        self.size_nodes = 10
        self.size_supports = 100
        self.size_force = 200
        self.scale_lcp = None
        self.text = text

    def plot_element(self, plt, element, u=None):
        """
        Function for plotting one element
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param element: element that we want to plot
        :param u: vector of global displacements
        :return: None
        """
        # get coordinates of the nodes of the element that we want to draw
        if u is not None:  # if we have displacements get new coordinates
            nodes = self.__create_new_deformed_nodes(u)
        else:
            nodes = self.nodes
        x, y = element.nodes_coordinates(nodes)
        # if it is 4 node element, add another element in list to draw 4th line with plt.plot()
        # 4 node elements plot
        if isinstance(element, Element4NodeLinear):
            # lines will be moved a little to the center of the element
            x, y = _modify_coordinates_for_plotting_4node(x, y)
            x.append(x[0])  # to encircle the element
            y.append(y[0])  # to encircle the element
            plt.plot(x, y, color='grey', zorder=1)
            # if el_num is not None:  # this fot plotting number of the element 4node
            #     plt.text(np.average(x), np.average(y), str("%.0f" % el_num), color='blue')
        elif isinstance(element, ElementMacro):  # do not plot ME
            pass
        else:
            # null elements. (numsides, style, angle) -  markers
            size = 100
            if isinstance(element, ElementNull):
                angle_deg = element.alpha * 180 / np.pi  # get angle in degrees
                t = mpl.markers.MarkerStyle(marker=1)  # make marker
                t._transform = t.get_transform().rotate_deg(angle_deg)  # rotate it
                plt.scatter(x[0], y[0], s=size, color='gold', marker=t, zorder=3, alpha=0.5)  # plot
            # frame elements
            else:
                plt.plot(x, y, color='blue', linewidth=1.5, zorder=2)

    def plot_list_of_elements(self, plt, deformed=False, u=None):
        """
        Function for plotting list of elements
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param deformed: if True then it uses self.nodes_def to plot def scheme, else it uses self.nodes
        :param u: vector of global displacements
        :return:
        """
        if deformed is False:
            list_of_elements = list_of_all_elements()  # get all elements in scheme
            for element in list_of_elements:
                self.plot_element(plt, element, u)
        else:
            if self.element_container_obj is not None:
                for element in self.element_container_obj.elements_list:
                    self.plot_element(plt, element, u)

    def plot_nodes(self, plt, deformed=False, linear=True):
        """
        Function for plotting all nodes in scheme
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param deformed: if True then use self.nodes_def to plot def scheme? else use self.nodes
        :param linear: if True than use self.nodes_def_linear, else use self.nodes_contact
        """
        if deformed:
            if linear:
                nodes = self.__create_new_deformed_nodes(self.u_linear)
            else:
                nodes = self.__create_new_deformed_nodes(self.u_contact_anim[-1])
        else:
            nodes = self.nodes

        for node in nodes:
            plt.scatter(node.x, node.y, s=self.size_nodes, color='lime', marker='s')

    def plot_supports(self, plt, supports):
        """
        Function for plotting supports (boundary conditions)
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param supports: support matrix
        :return:
        """
        for node in self.nodes:
            arr_index = np.intersect1d(node.indices, supports)  # find all intersection of 2 arrays
            for i in arr_index:
                order_of_index = node.indices.index(i)  # find 0-horizontal, 1-vertical or 2-rotation degree of freedom
                if order_of_index == 0:
                    angle_deg = 0
                    t = mpl.markers.MarkerStyle(marker=0)  # make marker
                    t._transform = t.get_transform().rotate_deg(angle_deg)  # rotate it
                    plt.scatter(node.x, node.y, s=self.size_supports, color='red', marker=t, zorder=3)  # plot
                elif order_of_index == 1:
                    angle_deg = 90
                    t = mpl.markers.MarkerStyle(marker=0)  # make marker
                    t._transform = t.get_transform().rotate_deg(angle_deg)  # rotate it
                    plt.scatter(node.x, node.y, s=self.size_supports, color='red', marker=t, zorder=3)  # plot
                elif order_of_index == 2:
                    plt.scatter(node.x, node.y, s=self.size_supports, facecolors='none', edgecolors='red', zorder=3)

    def plot_external_forces(self, plt):
        """
        Plotting external forces on matplotlib.pyplot of axis
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :return: None
        """
        for force, dof in zip(self.lv.rf, range(len(self.lv.rf))):  # dof - degree of freedom
            if force != 0:  # if we have load on some i dof
                for node in self.nodes:  # iterate over all nodes (nodes have coordinates)
                    if dof in node.indices:  # check which node have this dof
                        if node.indices[0] == dof:  # check the direction in which to put the marker
                            if force > 0:
                                angle_deg = 0
                            else:
                                angle_deg = 180
                        else:
                            if force > 0:
                                angle_deg = 90
                            else:
                                angle_deg = 270
                        t = mpl.markers.MarkerStyle(marker=1)  # make marker
                        t._transform = t.get_transform().rotate_deg(angle_deg)  # rotate it
                        plt.scatter(node.x, node.y, s=self.size_force, color='green', marker=t, zorder=3, alpha=0.7)

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
        :return:
        """
        # get nodes coordinates of the frame element
        [x1_init, x2_init], [y1_init, y2_init] = one_frame_el.nodes_coordinates(self.nodes)
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
                  f'please reduce SCALE_DEF parameter of make external load lower,\n'
                  f'results may be unreadable')
        # initial coordinate system. Create array of points to rotate and move
        xi = np.linspace(0, L, self.partition)
        yi = _frame_agg(xi, L, v1, fi1, v2, fi2)  # using 4 form functions of frame element
        # modified coordinate system. Create array of points to plot
        x_modified = one_frame_el.cosa * xi - one_frame_el.sina * yi + x1_init + w1 * one_frame_el.cosa
        y_modified = one_frame_el.sina * xi + one_frame_el.cosa * yi + y1_init + w1 * one_frame_el.sina
        return x_modified, y_modified

    def plot_frame_deformed(self, plt, u):
        """
        Plot all frame deformed elements
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param u: vector of global displacements
        :return: None, void type
        """
        for one_frame_element in self.element_frame:
            x, y = self.__create_points_for_deformed_frame(one_frame_element, u)
            plt.plot(x, y, color='blue', linewidth=1.5, zorder=2)

    def _get_frame_deformed(self, plt, u):
        """
        Get all frame deformed elements as image
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param u: vector of global displacements
        :return: image of all deformed frame elements
        """
        el_imgs = []
        for one_frame_element in self.element_frame:
            x, y = self.__create_points_for_deformed_frame(one_frame_element, u)
            el_imgs.append(plt.plot(x, y, color='blue', linewidth=1.5, zorder=2))
        res_img = el_imgs[0]
        for i in range(1, len(el_imgs)):
            res_img += el_imgs[i]
        return res_img

    def plot_deformed_scheme(self, plt, u=None):
        """
        Plot deformed scheme with 'element_container_obj' element type and element_frame.
        :param plt: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param u: vector of global displacements
        :return:
        """
        if self.element_frame is not None:
            self.plot_frame_deformed(plt, u)
        if self.element_container_obj is not None:
            self.plot_list_of_elements(plt, deformed=True, u=u)

    def plot_lcp_nt_i_step_lemke(self, plt_n, plt_t, i, text=True):
        """
        Plots interaction forces and mutual displacements along the normal and tangent to the contact zone
        :param plt_n: matplotlib.pyplot object or axis for plotting contact along the normal to the contact zone
        :param plt_t: matplotlib.pyplot object or axis for plotting contact along the tangent to the contact zone
        :param i: step number of Lemke algorithm
        :return: None
        """
        zn, zt, xn, xt = self.lemke.zn_anim[i], self.lemke.zt_anim[i], self.lemke.xn_anim[i], self.lemke.xt_anim[i]
        null_el_n = [element for element in self.element_null if element.orientation == 'n']
        null_el_t = [element for element in self.element_null if element.orientation == 't']
        x_range = range(len(null_el_n))
        # operations below are needed for that case if we have tangential null elements not in all contact pairs
        # so we make zt = 0, xt = 0, for frictionless contact pairs (without tangent null element)
        zt_tmp, xt_tmp = np.zeros(len(null_el_n)), np.zeros(len(null_el_n))
        for i, element in zip(range(len(null_el_t)), null_el_t):
            zt_tmp[element.contact_pair] = zt[i]
            xt_tmp[element.contact_pair] = xt[i]
        zt, xt = zt_tmp, xt_tmp

        if self.scale_lcp is None:
            if np.max(zn) > 1e-8:  # some smaller value will make plots unreadable
                self.scale_lcp = int(np.max(xn) / np.max(zn))
            elif np.max(zt) > 1e-8:
                self.scale_lcp = int(np.max(xt) / np.max(zt))
            else:
                self.scale_lcp = 1
        # plot for normal
        plt_n.bar(x_range, xn, color='red', width=0.2)  # plot bars
        plt_n.bar(x_range, zn * self.scale_lcp, color='blue', width=0.2)
        # plot for tangent
        plt_t.bar(x_range, xt, color='red', width=0.2)
        plt_t.bar(x_range, zt * self.scale_lcp, color='blue', width=0.2)
        plt_t.plot(x_range, xn*FRICTION_COEFFICIENT, color='orange')
        plt_t.plot(x_range, -xn * FRICTION_COEFFICIENT, color='orange')
        if self.text is True:
            # plot for normal
            for x, y in zip(x_range, xn):  # plot values
                plt_n.text(x, y, str("%.2f" % y), color='lightcoral')
            for x, y in zip(x_range, zn):
                plt_n.text(x, y * self.scale_lcp, str("%.4f" % y), color='cornflowerblue')
            # for tangent
            for x, y in zip(x_range, xt):  # plot values
                plt_t.text(x, y, str("%.2f" % y), color='lightcoral')
            for x, y in zip(x_range, zt):
                plt_t.text(x, y * self.scale_lcp, str("%.4f" % y), color='cornflowerblue')


    def plot_def_i_step_lemke(self, plt_def_contact, i):
        """
        Plot deformed scheme of Lemke algorithm
        :param plt_def_contact: matplotlib.pyplot object or axis for plotting def scheme
        :param i: step number
        :return: None
        """
        if self.element_frame is not None:
            self._get_frame_deformed(plt_def_contact, self.u_contact_anim[i])
        if self.element_container_obj is not None:
            self.plot_list_of_elements(plt_def_contact, deformed=True, u=self.u_contact_anim[i])

    def plot_nodes_i_step_lemke(self, plt_def_contact, i):
        """
        Function for plotting all nodes in scheme
        :param plt_def_contact: where to plot element. Must be Matplotlib.pyplot obj, or matplotlib.axes obj
        :param i: number of Lemke's algorithm step
        """
        nodes = self.__create_new_deformed_nodes(self.u_contact_anim[i])
        for node in nodes:
            plt_def_contact.scatter(node.x, node.y, s=self.size_nodes, color='lime', marker='s')

    def plot_scheme_i_step_lemke(self, plt_def_contact, plt_n, plt_t, i):
        """
        Plot all that we need for each step of Lemke's algorithm
        All plt_ variables must me instances of matplotlib.axes object (or maybe pyplot)
        :param plt_def_contact: axes for plotting deformed scheme
        :param plt_n: axes for plotting normal contact interaction forces and mutual displacements
        :param plt_t: axes for plotting tangent contact interaction forces and mutual displacements
        :param i: step of Lemke's algorithm to plot
        :return:
        """
        self.plot_def_i_step_lemke(plt_def_contact, i)
        self.plot_lcp_nt_i_step_lemke(plt_n, plt_t, i)
        self.plot_nodes_i_step_lemke(plt_def_contact, i)





