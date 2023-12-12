import math
import numpy as np
import time
from prettytable import PrettyTable

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

from input_data import FRICTION_COEFFICIENT
assert FRICTION_COEFFICIENT == 0.5, 'Friction coef need to be 0.5'  # 0.5

start = time.time()

# set inputs
Ar = 1  # 1
Er = 50000  # 50000
Ix = 1  # 1
F = 500  # 500  # 207.83048427075536 - F=1 нарастающий параметр
gap = 1  # 1
force_inc_alg = False

# add nodes # for 4 node element
nodes = NodeContainer()
# add frame nodes
nodes.add_node(0, 30)
nodes.add_node(10, 20)
nodes.add_node(20, 10)
nodes.add_node(30, 10)
nodes.add_node(40, 0)
# nodes for sup
nodes.add_node(10+gap, 20)
nodes.add_node(20, 10)
nodes.add_node(30, 10)
nodes.add_node(40, 0)

# add elements
element_4node = None

# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for nn in range(4):
    element_frame.add_element(EN=[nn, nn+1], E=Er, A=Ar, I=Ix)

# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[5, 1], cke=123, alpha=math.pi)
element_null.add_element(EN=[6, 2], cke=123, alpha=math.pi/2)
element_null.add_element(EN=[7, 3], cke=123, alpha=math.pi/2)
element_null.add_element(EN=[8, 4], cke=123, alpha=math.pi)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)

nodes_to_support = [5, 6, 7, 8]
sm.support_nodes(nodes_to_support, direction='hv')

lv = LoadVector()
lv_v = LoadVector()

if not force_inc_alg:
    lv.add_concentrated_force(force=F, degree_of_freedom=0)
    lv.add_concentrated_force(force=-F, degree_of_freedom=1)
else:
    lv.add_concentrated_force(force=1, degree_of_freedom=0)
    lv.add_concentrated_force(force=-1, degree_of_freedom=1)
    lv_v.add_concentrated_force(force=F, degree_of_freedom=0)
    lv_v.add_concentrated_force(force=-F, degree_of_freedom=1)


autorun = True
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=autorun)

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
