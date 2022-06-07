import math
import numpy as np
import time

from FEM.element_frame import ElementFrameContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_null import ElementNullContainer
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

start = time.time()  # to calculate time

nodes = NodeContainer()
nodes.add_node(0, 0)
nodes.add_node(0, 1)
nodes.add_node(3, 0)
nodes.add_node(3, 1)
number_of_floors = 10
for i in range(number_of_floors+1):
    for j in range(3):
        nodes.add_node(0.5+j, i)
# add frame
E_frame = 1e12  # Pa
A_frame = 0.2*0.2  # m^2
I_frame = 1*10**(-4)  # m^3
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[0, 1], E=E_frame, A=A_frame, I=I_frame)
element_frame.add_element(EN=[2, 3], E=E_frame, A=A_frame, I=I_frame)

# add 4node element
E_4node = 1e9
mu = 0.27
t = 2
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
for i in range(0, 3 * number_of_floors, 3):
    element_4node.add_element(EN=[i + 4, i + 5, i + 8, i + 7], E=E_4node, mu=mu, t=t)
    element_4node.add_element(EN=[i + 5, i + 6, i + 9, i + 8], E=E_4node, mu=mu, t=t)

cke = E_4node*t
element_null = ElementNullContainer(nodes_scheme=nodes)
# el_null.add_element(nodes, EN=[4, 0], cke=cke, alpha=math.pi)  # link to the left
element_null.add_element(EN=[7, 1], cke=cke, alpha=math.pi)  # link to the left
# el_null.add_element(nodes, EN=[6, 2], cke=cke, alpha=math.pi)  # link to the right
element_null.add_element(EN=[3, 9], cke=cke, alpha=math.pi)  # link to the right

F = 1e6  # N
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
# Supports
sm.support_nodes(list_of_nodes=[0, 2], direction='hvr')
sm.support_nodes([4, 5, 6], direction='v')
# set stiffness of rubber bearings
bearing_stiffness = 1e6
sm.add_spring(degree_of_freedom=8, stiffness=bearing_stiffness)
sm.add_spring(degree_of_freedom=10, stiffness=bearing_stiffness)
sm.add_spring(degree_of_freedom=12, stiffness=bearing_stiffness)

print('R rank:{}, R rows:{}'.format(np.linalg.matrix_rank(sm.r), sm.r.shape[0]))

lv = LoadVector()
for i in range(0, 3 * number_of_floors, 3):
    lv.add_concentrated_force(force=-F, degree_of_freedom=(i+5) * 2)  # adding load in center of building at each floor

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=True, force_incrementation=False)


print(f'xn sum: {sum(graph.lemke.xn)}, force sum = {2*F}'
      f'xt sum: {sum(graph.lemke.xt)}')



# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
