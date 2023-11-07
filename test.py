import numpy as np

from numba import jit

import pyqtgraph as pg
# import pyqtgraph.examples
# pyqtgraph.examples.run()

import math
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

import numpy as np


def simplex_method(c, A, b):
    # m, n = A.shape
    # # Добавляем фиктивные переменные
    # c = np.concatenate((c,np.zeros(m)))
    # A = np.concatenate((A, np.eye(m)), axis=1)
    # # Создаем таблицу симплекс-метода
    # table = np.concatenate((A, b.reshape(-1, 1)), axis=1)
    # table = np.concatenate((table, np.array(c)), axis=0)

    m, n = A.shape
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)
    # Add slack variables to convert inequality constraints to equality constraints
    A = np.hstack((A,np.eye(m)))
    c = np.hstack((c, np.zeros(m)))

    # Create initial tableau
    table = np.vstack((np.hstack((A,np.expand_dims(b, axis=1)))),np.hstack((c, np.array([0]))))

    while True:
        # Находим разрешающий столбец
        entering_col = np.argmin(table[-1, :-1])
        if table[-1, entering_col] >= 0:
            break

        # Находим разрешающую строку
        ratios = table[:-1, -1] / table[:-1, entering_col]
        leaving_row = np.argmin(ratios)

        # Обновляем таблицу симплекс-метода
        pivot = table[leaving_row, entering_col]
        table[leaving_row, :] /= pivot
        for i in range(m + 1):
            if i == leaving_row:
                ratio = table[i, entering_col]
                table[i, :] -= ratio * table[leaving_row, :]

        # Извлекаем результаты
    solution = table[:-1, -1]
    objective_value = -table[-1, -1]

    return solution, objective_value

c = np.array([-2, -3])
A = np.array([[1, 1], [2, 1], [1, 0]])
b = np.array([4, 5, 3])

solution, objective_value = simplex_method(c, A, b)
print("Solution:", solution)
print("Objective value:", objective_value)



