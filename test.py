import os
import datetime
import socket
import math
import time
import numpy as np
import pyqtgraph as pg
from scipy.optimize import minimize

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

from input_data import FRICTION_COEFFICIENT
assert FRICTION_COEFFICIENT == 0.5, 'Friction coef need to be 0.5'


start_time = time.time()

start = 0
stop = 100
char_limit = 3

y = []
def true_len(i, j, k):
    s = ''
    if i != 0: s += str(i)
    if j != 0: s += str(j)
    if k != 0: s += str(k)
    return len(s)

for i in range(start, stop):
    if i == 1: continue
    for j in range(start, stop):
        if j ==1: continue
        for k in range(start, stop):
            if k == 1: continue
            val = i*i + j*j + k*k
            if val not in y and true_len(i, j, k) < char_limit:
                y.append(val)
y = sorted(y)
x = list(range(len(y)))

# plot = pg.plot()
# line = plot.plot(arr[0], arr[1], pen='g', symbol='x', symbolPen='g',
#                   symbolBrush=0.2, name='green')




print(f'TIME: {time.time()-start_time}')

# importing Qt widgets
from PyQt5.QtWidgets import *
import sys

# importing pyqtgraph as pg
import pyqtgraph as pg
from PyQt5.QtGui import *


# Bar Graph class
class BarGraphItem(pg.BarGraphItem):

    # constructor which inherit original
    # BarGraphItem
    def __init__(self, *args, **kwargs):
        pg.BarGraphItem.__init__(self, *args, **kwargs)

    # creating a mouse double click event
    def mouseDoubleClickEvent(self, e):
        # setting scale
        self.setScale(0.2)


class Window(QMainWindow):

    def __init__(self):
        super().__init__()

        # setting title
        self.setWindowTitle("PyQtGraph")

        # setting geometry
        self.setGeometry(100, 100, 600, 500)

        # calling method
        self.UiComponents()

        # showing all the widgets
        self.show()

    # method for components
    def UiComponents(self):
        # creating a widget object
        widget = QWidget()

        # creating a new label
        label = QLabel("Line Plot")


        # create plot window object
        plt = pg.plot()

        # adding legend
        plt.addLegend()

        # set properties of the label for y axis
        plt.setLabel('left', 'number', units='i*i + j*j + k*k')

        # set properties of the label for x axis
        plt.setLabel('bottom', 'order', units='None')

        plt.plot(x, y, pen='b', symbol='x',
                    symbolPen='g', symbolBrush=0.3, name='blue')


        # setting this widget as central widget of the main window
        self.setCentralWidget(widget)


# create pyqt5 app
App = QApplication(sys.argv)

# create the instance of our Window
window = Window()

# start the app
sys.exit(App.exec())