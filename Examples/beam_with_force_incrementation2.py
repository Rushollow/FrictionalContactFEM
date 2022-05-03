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
length = 1  # meter
# nodes for frame
nodes.add_node(0,0)  # 0
nodes.add_node(length,0) # 1
nodes.add_node(length*2,0)  # 2
nodes.add_node(length*3,0)  # 3
nodes.add_node(length*3.5,0)  # 4
# nodes for supports

nodes.add_node(length,0)  # 5
nodes.add_node(length*2,0)  # 6
nodes.add_node(length*3,0)  # 7


# set inputs
Ar = 1
Er = 1
Ix = 1
E = 1
F = 2  # Newtons 2

# add elements
element_4node = None
# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(4):
    element_frame.add_element(EN=[i, i+1], E=Er, A=Ar, I=Ix)

# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[5, 1], cke=123, alpha=math.pi/2, gap_length=0)
element_null.add_element(EN=[6, 2], cke=123, alpha=math.pi/2, gap_length=0)
element_null.add_element(EN=[7, 3], cke=123, alpha=math.pi/2, gap_length=0)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_support = [5, 6, 7]
sm.support_nodes(nodes_to_support, direction='hv')
sm.support_nodes(list_of_nodes=[0], direction='hvr')
lv = LoadVector()
lv.add_concentrated_force(force=-F, degree_of_freedom=4*2+1)


# get initial table to Excel file
from calculation import Calculate
calc = Calculate(nodes=nodes, sm=sm, lv_const=lv, element_frame=element_frame, element_null=element_null)
#calc.table_to_excel('initial_table.xlsx')


# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, element_frame=element_frame, element_container_obj=element_4node,
                   element_null=element_null, partition=10, scale_def=1, force_incrementation=False, autorun=True)

for i in range(len(graph.lemke.zn_anim)):
    print(f'p:{graph.lemke.p_anim[i]} zn:{graph.lemke.zn_anim[i]} xn:{graph.lemke.zt_anim[i]}')
# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
