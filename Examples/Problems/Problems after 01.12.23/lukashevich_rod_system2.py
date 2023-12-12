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
assert FRICTION_COEFFICIENT == 0.3, 'Friction coef need to be 0.3'

start = time.time()

# set inputs
Ar = 1
Er = 50_000_000
Ix = 1
F1 = 45_000  # Н
F2 = 15_000
F3 = 5_000
P1 = 40_000
P2 = 20_000

# add nodes # for 4 node element
nodes = NodeContainer()
# add frame nodes
nodes.add_node(0, 5, rotation=True)  # 0
nodes.add_node(10, 0, rotation=True)  # 1
nodes.add_node(20, 0, rotation=True)  # 2
nodes.add_node(30, 5, rotation=True)  # 3
nodes.add_node(40, 5, rotation=True)  # 4
#add sup nodes
nodes.add_node(10, 0)  # 5
nodes.add_node(20, 0)  # 6
nodes.add_node(30, 5)  # 7

# add elements
element_4node = None

# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for nn in range(4):
    element_frame.add_element(EN=[nn, nn+1], E=Er, A=Ar, I=Ix)

# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[5, 1], cke=123, alpha=math.pi / 2)
element_null.add_element(EN=[6, 2], cke=123, alpha=math.pi / 2)
element_null.add_element(EN=[7, 3], cke=123, alpha=math.pi / 2)

# add hinges
nodes.add_hinge(1, element_frame_container=element_frame)
nodes.add_hinge(2, element_frame_container=element_frame)
nodes.add_hinge(3, element_frame_container=element_frame)

for i, el in enumerate(element_frame):
    print(f'Frame №:{i}, EN:{el.EN}')
for node in nodes:
    print(f'node {node.number} x:{node.x} y:{node.y} Indxs:{node.indices}')

springs = False
# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.support_nodes([0], direction='v')
if not springs:
    sm.support_nodes([4, 5, 6, 7], direction='hv')
else:  # add springs so problem would be most like in Lukashevich algorithm
    sm.support_nodes([4], direction='hv')
    sm.add_spring(degree_of_freedom=15, stiffness=Er*Ar/10)
    sm.add_spring(degree_of_freedom=17, stiffness=Er*Ar/12.4)
    sm.add_spring(degree_of_freedom=19, stiffness=Er*Ar/2)

    sm.add_spring(degree_of_freedom=16, stiffness=Er*Ar/10)
    sm.add_spring(degree_of_freedom=18, stiffness=Er*Ar/15)
    sm.add_spring(degree_of_freedom=20, stiffness=Er*Ar/11)
    print(f'Жёсткость пружин по нормали: \n {Er*Ar/10, Er*Ar/12.4, Er*Ar/2}')
    print(f'Жёсткость пружин по касательной: \n {Er*Ar/10, Er*Ar/15, Er*Ar/11}')


lv = LoadVector()
lv_v = LoadVector()
force_inc = True
if not force_inc:
    lv.add_concentrated_force(force=P1, degree_of_freedom=0)
    lv.add_concentrated_force(force=-P2, degree_of_freedom=6)
    lv.add_concentrated_force(force=-F1, degree_of_freedom=4)
    lv.add_concentrated_force(force=-F2, degree_of_freedom=7)
    lv.add_concentrated_force(force=-F3, degree_of_freedom=10)
else:
    lv_v.add_concentrated_force(force=P1, degree_of_freedom=0)
    lv_v.add_concentrated_force(force=-P2, degree_of_freedom=6)
    lv.add_concentrated_force(force=-F1, degree_of_freedom=4)
    lv.add_concentrated_force(force=-F2, degree_of_freedom=7)
    lv.add_concentrated_force(force=-F3, degree_of_freedom=10)

autorun = True
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=300, autorun=autorun)

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
