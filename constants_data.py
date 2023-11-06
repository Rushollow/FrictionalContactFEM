import input_data
import math
import numpy as np
from scipy.misc import derivative

# On/Off debug
DEBUG = True

# Soil types
class Soil:

    def __init__(self):
        if input_data.SOIL_TYPE == 'rock':
            self.E = 1
            self.gamma = 1.1

        if input_data.SOIL_TYPE == 'sand':
            self.E = 1
            self.gamma = 1


# Arch types
class Arch:

    def __init__(self):
        self.arch_span = input_data.ARCH_SPAN
        self.arch_height = input_data.ARCH_HEIGHT

    def _initial_func_parameters(self, x):
        f = self.arch_height
        l = self.arch_span
        x += l/2
        return f, l, x

    def outline_sinus(self, x):
        f, l, x = self._initial_func_parameters(x)
        return f * np.sin(np.pi * x / l)

    def outline_circle(self, x):
        f, l, x = self._initial_func_parameters(x)
        if 2 * f > l:
            f = l / 2
            print(f'Arch height f:{self.arch_height} is higher than its span divided by two l/2:{l/2}\n'
                  f'-> Taken value of arch height is {f}')
        r = f / 2 + l*l / (8 * f)
        return (r*r - (l / 2 - x)*(l / 2 - x))**0.5 - r + f

    def outline_parabola(self, x):
        f, l, x = self._initial_func_parameters(x)
        return 4 * f * x / l * (1 - x / l)

    def outline_ellipse(self, x):
        f, l, x = self._initial_func_parameters(x)
        return 2 * f * np.sqrt(x/l * (1 - x/l))

