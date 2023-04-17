import numpy as np

from numba import jit

import pyqtgraph as pg
# import pyqtgraph.examples
# pyqtgraph.examples.run()
import copy

class Foo():
    def __init__(self):
        self.x=1
        self.y=1

    def copy(self):
        return copy.copy(self)


f1 = Foo()
f1.x = 2
print(f1.x)
f2 = f1.copy()
f1.x = 3
print(f2.x)


