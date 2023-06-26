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

# set inputs
Ar = 1  # Площадь поперечного сечения рамы
Er = 1  # Модуль упругости рамы
Ix = 1  # Момент инерции сечения рамы
F = 1   # Величина внешней нагрузки

nodes = NodeContainer()
# узлы для рамы
nodes.add_node(0, 0)  # 0 Узел
nodes.add_node(1, 0)  # 1 Узел
nodes.add_node(2, 0)  # 2 Узел
nodes.add_node(3, 0)  # 3 Узел

# Узлы на которые ставятся опоры
nodes.add_node(0, 0)  # 4 Узел
nodes.add_node(2, 0)  # 5 Узел
nodes.add_node(3, 0)  # 6 Узел

# Добавляются элементы
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(3):
    element_frame.add_element(EN=[i, i+1], E=Er, A=Ar, I=Ix)
# Добавляются нуль-элементы для односторонних связей
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[4, 0], cke=1, alpha=math.pi/2, add_t_el=True)
element_null.add_element(EN=[5, 2], cke=1, alpha=math.pi/2, add_t_el=True)
element_null.add_element(EN=[6, 3], cke=1, alpha=math.pi/2, add_t_el=True)

# Матрица жёсткости
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=None, el_null=element_null)
sm.support_nodes(list_of_nodes=[4, 5, 6], direction='hv')  # опирание узлов

SITUATION = 2
lv_const = LoadVector()
lv_const.add_concentrated_force(force=-F, degree_of_freedom=3)

# передать все данные для расчёта
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=None,
                   element_frame=element_frame, element_container_obj=None, element_null=element_null,
                   partition=10, scale_def=1, autorun=True)

# запуск приложения
if __name__ == "__main__":
    graph.fill_arrays_scheme()  # сформировать данные для графики
    application(graph)



