import math
import numpy as np
import time

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

start = time.time()

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# add nodes # for 4 node element
nodes = NodeContainer()
length = 1  # 1  # meter
n = 5  # amount of nodes of frame
# nodes for frame
nodes.add_node(0, 0)  # 0
nodes.add_node(length, 0)  # 1
nodes.add_node(length*2, 0)  # 2
nodes.add_node(length*3, 0)  # 3
nodes.add_node(length*4, 0)  # 4

# node for support
nodes.add_node(length*0, 0)  # 5
nodes.add_node(length*2, 0)  # 6
nodes.add_node(length*4, 0)  # 7

# set inputs
Ar = 1  # 1              4*7e-4
Er = 1  # 1         2.04e11
Ix = 0.5  # 0.5          (4e-8*(7**3))/12
F = 1  # 1  # Newtons        10e3

# add elements
element_4node = None
# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(4):
    element_frame.add_element(EN=[i, i+1], E=Er, A=Ar, I=Ix)
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[5, 0], cke=1, alpha=math.pi/2, add_t_el=True)
element_null.add_element(EN=[6, 2], cke=1, alpha=math.pi/2, add_t_el=True, gap_length=0)
element_null.add_element(EN=[7, 4], cke=1, alpha=math.pi/2, add_t_el=True)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.support_nodes(list_of_nodes=[5, 6, 7], direction='hv')  # sup for unilateral
# HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SITUATION = 5
lv_const = LoadVector()
lv_variable = None
if SITUATION == 1:  # just 2 const force
    lv_const.add_concentrated_force(force=-F, degree_of_freedom=3)
    lv_const.add_concentrated_force(force=-F/5, degree_of_freedom=8)
elif SITUATION == 2:  # one variable load
    lv_variable = LoadVector()
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=7)
elif SITUATION == 3:  # one variable load upward and 2 const forces
    lv_const.add_concentrated_force(force=-F, degree_of_freedom=3)
    lv_const.add_concentrated_force(force=-F/5, degree_of_freedom=8)
    lv_variable = LoadVector()
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=7)
elif SITUATION == 4:
    lv_const.add_concentrated_force(force=-F, degree_of_freedom=3)
    lv_const.add_concentrated_force(force=-F / 5, degree_of_freedom=8)
    lv_variable = LoadVector(vectors_amount=2)
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=7, vector_num=0)
    lv_variable.add_concentrated_force(force=F, degree_of_freedom=7, vector_num=1)
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=9, vector_num=1)
elif SITUATION == 5:  #
    lv_variable = LoadVector(vectors_amount=2)
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=3, vector_num=0)
    lv_variable.add_concentrated_force(force=-F, degree_of_freedom=7, vector_num=1)

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=True)

for i in range(len(graph.lemke.zn_anim)):
    print(f'{i}: p:{graph.lemke.p_anim[i]} zn:{graph.lemke.zn_anim[i]} xn:{graph.lemke.xn_anim[i]}'
          f'zt:{graph.lemke.zt_anim[i]} xt:{graph.lemke.xt_anim[i]}')
# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
