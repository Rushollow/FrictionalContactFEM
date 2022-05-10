import numpy as np

from numba import jit

import pyqtgraph as pg
# import pyqtgraph.examples
# pyqtgraph.examples.run()


# import sys
# from PyQt5 import QtWidgets, QtCore
#
# app = QtWidgets.QApplication(sys.argv)
# widget = QtWidgets.QWidget()
# widget.resize(400, 200)
# widget.setWindowTitle("This is PyQt Widget example")
# widget.show()
# exit(app.exec_())

vector = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18]]
a = []
print([i for vec in vector for i in vec])