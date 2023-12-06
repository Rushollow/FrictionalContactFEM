import math
import numpy as np
import time

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
Ar = 0.2*0.2
Er = 2e10
Ix = 0.2**4/12
F1 = 45000  # Н
F2 = 15000
F3 = 5000
P1 = 40000 * 2
P2 = 20000

# add nodes # for 4 node element
nodes = NodeContainer()
# add frame nodes
nodes.add_node(0, 5, rotation=False)  # 0
nodes.add_node(10, 0, rotation=False)  # 1
nodes.add_node(20, 0, rotation=False)  # 2
nodes.add_node(30, 5, rotation=False)  # 3
nodes.add_node(40, 5, rotation=False)  # 4
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


for i, el in enumerate(element_frame):
    print(f'Frame №:{i}, EN:{el.EN}')
for node in nodes:
    print(f'node {node.number} x:{node.x} y:{node.y} Indxs:{node.indices}')


# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_support = [4, 5, 6, 7]
sm.support_nodes(nodes_to_support, direction='hv')
sm.support_nodes([0], direction='v')
lv = LoadVector()
lv.add_concentrated_force(force=P1, degree_of_freedom=0)
lv.add_concentrated_force(force=-P2, degree_of_freedom=4)
lv.add_concentrated_force(force=-F1, degree_of_freedom=3)
lv.add_concentrated_force(force=-F2, degree_of_freedom=5)
lv.add_concentrated_force(force=-F3, degree_of_freedom=7)


# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=True)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)


for i in range(len(graph.lemke.zn_anim)):
    print(f'{i} zn:{graph.lemke.zn_anim[i]} xn:{graph.lemke.zt_anim[i]}')


if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
