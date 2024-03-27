import math
import numpy as np
import time
from prettytable import PrettyTable

from input_data import SCALE_DEF
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_frame import ElementFrameContainer
from LCP.initial_table import InitialTable
from LCP.lemke import Lemke
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

start = time.time()

# set inputs
Ar = 0.2  # Cross-section area
Er = 2e11  # Yong's modulus
Ix = 2.3 * 10 ** (-4)  # Moment of inertia
F = 1e7  # Force external
# add nodes ---------------------------------------------------------------------
nodes = NodeContainer()
for i in range(6):  # adding nodes every 1 m (along x axis)
    nodes.add_node(i, 0)
    nodes.add_node(i, 0)
nodes.add_node(2.5, 0)  # add node in the center
# add frame elements -------------------------------------------------------------
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[0, 2], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[2, 4], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[4, 12], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[12, 6], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[6, 8], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[8, 10], E=Er, A=Ar, I=Ix)
# add 4 node elements
element_4node = None  # there are no such elements
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
for i in range(0, 12, 2):
    element_null.add_element(EN=[i + 1, i], cke=123, alpha=math.pi / 2)


# form R and solve SLAE --------------------------------------------------------
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=[], el_null=element_null)
nodes_to_support = [i for i in range(1, 12, 2)]
sm.support_nodes(nodes_to_support, direction='hv')
print('rank and shape0:', np.linalg.matrix_rank(sm.r), sm.r.shape)
print(f'cond: {np.linalg.cond(sm.r)}')
#assert np.linalg.matrix_rank(sm.r) == sm.r.shape[0]  # check if system in construction

lv_const = LoadVector()
lv_const.add_concentrated_force(force=-F, degree_of_freedom=25)
#lv.add_concentrated_force(force=-F/100, degree_of_freedom=21)
lv_variable = None

autorun = True
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=50, autorun=autorun)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    print(mytable)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)



