#!/urs/bin/python3
# -*- coding: utf-8 -*-
#  pyuic5 name.ui -o name.py - запускаем из папки с файлом ui в conda prompt (if anaconda is not installed then use cmd)
# pyuic5 ShowPlotWindow.ui -o plot_window.py

import copy
import math
import time
from PyQt5 import QtCore
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtCore import QRect, QPoint
import sys
import os
import numpy as np
import pyqtgraph as pg
from Visualize.plot_data_qt import PlotScheme
from input_data import FRICTION_COEFFICIENT
from input_data import BACKGROUND_COLOR


class MainWindow(QMainWindow):
    """
    Main window class
    """
    def __init__(self, graph=None):
        """
        Constructor of main window
        :param graph: info to plot (from Visualize.plot_data)
        """

        super(MainWindow, self).__init__()
        root = os.path.dirname(os.path.realpath(__file__))
        self.ui = uic.loadUi(os.path.join(root, 'ContactFEM.ui'), baseinstance=self)  # user interface from QT designer
        self.ui.setStatusBar(self.ui.statusbar)  # Activate statusBar for additional info
        self.lemke_step_shown = None  # number of lemke step shown on

        # widget and item that double click event is working with
        self.active_widget = None
        self.active_plot_item = None
        # double click and new window
        self.new_window = None  # plot new window on full screen
        pg.PlotWidget.mouseDoubleClickEvent = self.mouseDoubleClickEvent

        if graph is not None:
            self.graph = graph
        else:
            try:
                self.graph = PlotScheme(None, None, None)  # Just empty thing (helps to see PlotScheme methods in IDE)
            except AttributeError: print('graph data was not passed to UI (nodes, sm, elements and etc.')

        if BACKGROUND_COLOR is not None:
            self.bg_color = BACKGROUND_COLOR  # choose background color
        self._connect_actions_menu_bar()  # connect actions so they could work
        self._show_scheme()

    def _connect_actions_menu_bar(self):
        """
        Connect menu bar actions to functions
        """
        self.ui.actionRun.triggered.connect(self.run)

    def _show_scheme(self):
        """
        Method to show the initial scheme user set by user
        Scheme showed in graphicsView (QT)
        :return: nothing
        """
        # adding data and connect for scheme
        # nodes
        item_nodes = pg.ScatterPlotItem(pos=self.graph.arr_nodes_pos, size=7, symbol='s', brush='g',
                                        data=np.arange(self.graph.arr_nodes_pos.shape[0]),
                                        tip='x: {x:.3g}\ny: {y:.3g}\nnum={data}'.format,
                                        hoverable=True,
                                        #pxMode=False,  # Set pxMode=False to allow spots to transform with the view
                                        # hoverPen=pg.mkPen('g'),
                                        hoverSize=1e-1
                                        )
        # frame
        if self.graph.arr_frame_en is not None:
            frame_pen = pg.mkPen(color='b', width=5)
            item_frame_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos, adj=self.graph.arr_frame_en,
                                               symbolBrush=None,
                                               pen=frame_pen, symbolPen=None)
            self.ui.graphicsView_Scheme.addItem(item_frame_elements)
        # 4node elements
        if self.graph.arr_4node_en is not None:
            item_4node_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos, adj=self.graph.arr_4node_en,
                                               symbolBrush=None, symbolPen=None)
            self.ui.graphicsView_Scheme.addItem(item_4node_elements)
        # null-elements
        for pos_angle in self.graph.arr_null_el_1st_nodes_pos_angle:
            arrow = pg.ArrowItem(pos=pos_angle[:2], angle=-pos_angle[2], brush='white', headWidth=1, headLen=15)
            self.ui.graphicsView_Scheme.addItem(arrow)
        # supports
        for pos_angle in self.graph.arr_supp_pos_angle:
            arrow = pg.ArrowItem(pos=pos_angle[:2], angle=pos_angle[2], brush='red', headWidth=5, headLen=10)
            self.ui.graphicsView_Scheme.addItem(arrow)
        # add nodes in the end
        self.ui.graphicsView_Scheme.addItem(item_nodes)
        # set parameters
        self.ui.graphicsView_Scheme.setAspectLocked()
        # set background color
        self.ui.graphicsView_Scheme.setBackground(self.bg_color)
        self.ui.graphicsView_Scheme_deformed.setBackground(self.bg_color)
        self.ui.graphicsView_Scheme_unilateral_contact.setBackground(self.bg_color)
        self.ui.graphicsView_contact_info_tangent.setBackground(self.bg_color)
        self.ui.graphicsView_contact_info_normal.setBackground(self.bg_color)


    def run(self):
        """
        Run the calculation of the problem.
        Solving the problem and adding data to 'graphicsView' (QT)
        :return: nothing
        """
        # check if calculation is needed:
        if self.graph.lemke is None:  # if there is a data, then problem is solved (probably)
            self.graph.calculate_all()  # calculate the problem
        self.graph.fill_arrays_deformed_linear()  # form data to plot linear
        self._show_scheme_deformed()  # add data to UI linear deformation
        self.graph.fill_arrays_deformed_contact()  # form data to plot contact
        self._show_scheme_deformed_contact()  # add data to UI contact deformation
        self._show_contact_forces_and_displacements()  # add data to UI contact mutual forces and displacements
        # set parameters
        self.ui.graphicsView_Scheme_deformed.setAspectLocked()
        self.ui.graphicsView_Scheme_unilateral_contact.setAspectLocked()
        # remember shown Lemke's algorithm step
        self.lemke_step_shown = len(self.graph.u_contact_anim) - 1

    def _show_scheme_deformed(self):
        """
        Show deformed scheme.
        Scheme showed in graphicsView (QT)
        :return: nothing
        """
        self.ui.graphicsView_Scheme_deformed.clear()  # clear previous data
        # creating plot items and connecting them
        # frame
        if self.graph.arr_frame_en_deformed is not None:
            frame_pen = pg.mkPen(color='b', width=5)
            item_frame_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos_deformed_frame,
                                               adj=self.graph.arr_frame_en_deformed,
                                               symbolBrush=None, pen=frame_pen, symbolPen=None)
            self.ui.graphicsView_Scheme_deformed.addItem(item_frame_elements)
        # 4node
        if self.graph.arr_4node_en is not None:
            item_4node_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos_deformed, adj=self.graph.arr_4node_en,
                                           symbolBrush=None, symbolPen=None)
            self.ui.graphicsView_Scheme_deformed.addItem(item_4node_elements)
        # nodes
        item_nodes = pg.ScatterPlotItem(pos=self.graph.arr_nodes_pos_deformed, size=7, symbol='s', brush='g',
                                        data=self.graph.nodes,  # list of nodes objects
                                        tip='dx: {data.dx:.3g}\ndy: {data.dy:.3g}\n num={data.number}'.format,
                                        hoverable=True,
                                        # pxMode=False,  # Set pxMode=False to allow spots to transform with the view
                                        hoverPen=pg.mkPen('r'),
                                        hoverSize=1e-2)
        self.ui.graphicsView_Scheme_deformed.addItem(item_nodes)

    def _show_scheme_deformed_contact(self, i_step=-1):
        """
        Showing chosen i step of Lemke's algorithm
        :param i_step: number of the step. -1 means it will show last step
        :return: nothing
        """
        # clear data
        self.ui.graphicsView_Scheme_unilateral_contact.clear()
        # form data to plot contact
        self.graph.fill_arrays_deformed_contact(i_step=i_step)
        # displaying data
        # 4node elements
        if self.graph.arr_4node_en is not None:
            item_4node_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos_contact, adj=self.graph.arr_4node_en,
                                               symbolBrush=None, symbolPen=None)
            self.ui.graphicsView_Scheme_unilateral_contact.addItem(item_4node_elements)
        if self.graph.arr_nodes_pos_contact_frame is not None:
            frame_pen = pg.mkPen(color='b', width=5)
            item_frame_elements = pg.GraphItem(pos=self.graph.arr_nodes_pos_contact_frame,
                                               adj=self.graph.arr_frame_en_deformed,
                                               symbolBrush=None, pen=frame_pen, symbolPen=None)
            self.ui.graphicsView_Scheme_unilateral_contact.addItem(item_frame_elements)
        # nodes
        item_nodes_def = pg.GraphItem(pos=self.graph.arr_nodes_pos_contact, symbol='s', brush='g',
                                      data=self.graph.nodes,  # list of nodes objects
                                      tip='dx: {data.dx_c:.3g}\n dy: {data.dy_c:.3g}\n num={data.number}'.format,
                                      hoverable=True,
                                      # pxMode=False,  # Set pxMode=False to allow spots to transform with the view
                                      hoverPen=pg.mkPen('r'),
                                      hoverSize=1e-2)
        self.ui.graphicsView_Scheme_unilateral_contact.addItem(item_nodes_def)

    def _show_contact_forces_and_displacements(self, i_step:[int]=-1, text=True):
        """
        Show contact forces xn, xt and displacements zn, zt
        Data showed in graphicsView (QT)
        :param i_step: step if Lemke's algorithm to show. -1 is the last one
        :return: nothing
        """
        # clear all previous data
        self.ui.graphicsView_contact_info_normal.clear()
        self.ui.graphicsView_contact_info_tangent.clear()
        n_range = np.arange(len(self.graph.lemke.zn))  # amount of contact pairs (nodes that could be in contact)
        bar_width = 0.2  # width of bars...
        zn, zt = self.graph.lemke.zn_anim[i_step], self.graph.lemke.zt_anim[i_step]  # mutual displacement
        xn, xt = self.graph.lemke.xn_anim[i_step], self.graph.lemke.xt_anim[i_step]  # contact forces
        # scale up mutual displacements, so they could be seen on the same chart with contact forces
        scale_lcp_n = 1
        scale_lcp_t = 1
        if zn.any() and np.max(zn) > 1e-8:  # some smaller value will make plots unreadable
            small_range = np.max(zn) - np.min(xn)
            large_range = np.max(xn) - np.min(zn)
            small_max_abs = np.max(np.abs(zn))
            scale_lcp_n = large_range / (1 * small_max_abs)
        if zt.any() and (np.max(zt) > 1e-8 or np.min(zt) < -1e-8) :
            small_range = np.max(zt) - np.min(xt)
            large_range = np.max(xt) - np.min(zt)
            small_max_abs = np.max(np.abs(zt))
            scale_lcp_t = large_range / (1 * small_max_abs)
        self.ui.statusbar.showMessage(f'Scale of LCP is {scale_lcp_n, scale_lcp_t}')
        # normal
        if zn.any():
            zn_item = pg.BarGraphItem(x=n_range, height=zn * scale_lcp_n, width=bar_width, brush='b',
                                      data=zn,
                                      tip='x: {x:.3g}\ny: {y:.3g}\nvalue={data}'.format,
                                      hoverable=True)
            zn_item.setToolTip(str(zn))  # TODO here, adding tooltip to bar graph
            self.ui.graphicsView_contact_info_normal.addItem(zn_item)
        if xn.any():
            xn_item = pg.BarGraphItem(x=n_range, height=xn, width=bar_width, brush='r')
            self.ui.graphicsView_contact_info_normal.addItem(xn_item)
        # tangent
        if zt.any():
            zt_item = pg.BarGraphItem(x=n_range, height=zt * scale_lcp_t, width=bar_width, brush='b')
            self.ui.graphicsView_contact_info_tangent.addItem(zt_item)
        if xt.any():
            xt_item = pg.BarGraphItem(x=n_range, height=xt, width=bar_width, brush='r')
            self.ui.graphicsView_contact_info_tangent.addItem(xt_item)
        # adding ultimate forces line (to tangent contact forces, friction)
        # ultimate_forces_item1 = pg.PlotCurveItem(n_range, xn*FRICTION_COEFFICIENT,  pen=pg.mkPen(color=(255,170,100), width=2))
        # ultimate_forces_item2 = pg.PlotCurveItem(n_range, -xn*FRICTION_COEFFICIENT, pen=pg.mkPen(color=(255,170,100), width=2))
        # self.ui.graphicsView_contact_info_tangent.addItem(ultimate_forces_item1)
        # self.ui.graphicsView_contact_info_tangent.addItem(ultimate_forces_item2)

        bg1 = pg.BarGraphItem(x=n_range, y=xn*FRICTION_COEFFICIENT, height=0.01, width=0.9, color='orange')
        bg2 = pg.BarGraphItem(x=n_range, y=-xn*FRICTION_COEFFICIENT, height=0.01, width=0.9, color='orange')
        self.ui.graphicsView_contact_info_tangent.addItem(bg1)
        self.ui.graphicsView_contact_info_tangent.addItem(bg2)

        if text is True:  # TODO this
            return
            n = 1
            x, y = self.graph.arr_nodes_pos[n]
            text_item = pg.TextItem(
                html=f'<span style="color: #FFF;font-size: 11pt">{n}</span><br></div>',
                anchor=(0.0, 0.6), angle=0, border=None, fill=None)
            text_item.setPos(x, y)
            self.ui.graphicsView_contact_info_normal.addItem(text_item)

    def mouseDoubleClickEvent(self, event):
        """
        Event on doubleclick in App
        On doubleclick in graphicsView(QT) widget opens new window with items from this widget
        to look closer / full-screen
        :param event:
        :return: nothing
        """
        print("Show pg.PlotWidget in another window (mouseDoubleClickEvent)")
        plot_item = None  # thing to plot in new window
        # get position where click happened in global coordinates
        pos = event.globalPos()
        for i in range(self.ui.gridLayout.count()):  # iterate over all widgets in layout
            widget = self.ui.gridLayout.itemAt(i).widget()
            if isinstance(widget, pg.PlotWidget):  # find plotWidgets
                size = widget.size()  # get size of plot
                global_point = widget.mapToGlobal(QPoint(0, 0))  # get global pos of top left point
                global_rect = QRect(global_point, size)  # create rectangle of the plot
                if global_rect.contains(pos):  # if clicked in rectangle:
                    plot_item = widget.getPlotItem()  # get item that widget has
                    self.active_widget = widget
                    self.active_plot_item = plot_item
                    break
        if plot_item is not None:  # if found clicked widget
            self.new_window = ShowPlotWindow(self.active_widget, self.active_plot_item)  # create new window
            # create plot widget with plot item in it
            self.new_window.ui.graphicsView = pg.PlotWidget(parent=self.new_window.ui.verticalLayoutWidget, plotItem=plot_item)
            # add widget to layout
            self.new_window.ui.verticalLayout.addWidget(self.new_window.ui.graphicsView)
            self.new_window.ui.graphicsView.setBackground(self.bg_color)
            self.new_window.show()  # show window

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Left:
            if self.lemke_step_shown - 1 >= 0:
                self.lemke_step_shown -= 1
                self.redraw_graphs(self.lemke_step_shown)
        elif event.key() == QtCore.Qt.Key_Right:
            if self.lemke_step_shown + 1 < len(self.graph.u_contact_anim):
                self.lemke_step_shown += 1
                self.redraw_graphs(self.lemke_step_shown)

    def redraw_graphs(self, step_to_show):
        self._show_scheme_deformed_contact(step_to_show)
        self._show_contact_forces_and_displacements(step_to_show)
        p = self.graph.lemke.p_anim[self.lemke_step_shown]  # parameter 'p' value
        if step_to_show not in self.graph.lemke.p_anim_variable:
            stage = 0
        else:
            stage = self.graph.lemke.p_anim_variable[step_to_show][0]
        self.ui.statusbar.showMessage(f'Shown {self.lemke_step_shown} step. p={p}. stage={stage}')


class ShowPlotWindow(QMainWindow):
    """
    Show window class
    This window used in doubleclick event (to look closer at data displayed on widget)
    """
    def __init__(self, active_widget, active_plot_item):
        """
        Constructor of show-window
        """
        super(ShowPlotWindow, self).__init__()
        root = os.path.dirname(os.path.realpath(__file__))
        self.ui = uic.loadUi(os.path.join(root, 'ShowPlotWindow.ui'), baseinstance=self)

        self.active_widget = active_widget
        self.active_plot_item = active_plot_item
        #self.showMaximized()

###Make it work! pass plot item back to initial PlotWidget of pyqtgraph
    # def closeEvent(self, event):
    #     dummy = pg.PlotItem(name='Copy of plot2')
    #     self.ui.graphicsView.plotItem = dummy
    #     print(self.ui.graphicsView.plotItem)
    #     print(self.active_plot_item)
    #     plot_widget = pg.PlotWidget(plotItem=self.active_plot_item)
    #     view_box = plot_widget.plotItem.getViewBox()
    #     # print(plot_widget.plotItem)
    #     # self.active_widget.plotItem.vb = view_box





def application(graph=None):
    """
    Execution of the application (show window to user)
    :param graph: Information to plot - instance of Visualize.plot_data_scheme_qt.PlotScheme()
    :return: nothing
    """
    if graph is None:
        raise AttributeError('application need graph object')
    app = QApplication(sys.argv)  # create object of our app and pass system data to it (info about computer)
    window = MainWindow(graph)

    window.show()
    sys.exit(app.exec_())

def min_max_normalize(np_vec):
    min_val = min(np_vec)
    max_val = max(np_vec)
    normalized_vector = (np_vec - min_val) / (max_val - min_val)
    return normalized_vector
