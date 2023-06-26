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

# Add nodes
nodes = NodeContainer()
nodes.add_node(0, 0.5, rotation=True)
nodes.add_node(0, 1.5, rotation=True)
nodes.add_node(1, 1.5, rotation=True)
nodes.add_node(0.5, 0)
nodes.add_node(1.5, 0)
nodes.add_node(1.5, 1)
nodes.add_node(0.5, 1)
nodes.add_node(0, 0, rotation=True)
nodes.add_node(1.5, 1.5, rotation=True)

# Add frame
E_frame = 1e12  # Pa
A_frame = 0.1  # m^2
I_frame = 1e-4  # m^3
F = 1e9  # N
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[0, 1], E=E_frame, A=A_frame, I=I_frame)
element_frame.add_element(EN=[1, 2], E=E_frame, A=A_frame, I=I_frame)
element_frame.add_element(EN=[7, 0], E=E_frame, A=A_frame, I=I_frame)
element_frame.add_element(EN=[2, 8], E=E_frame, A=A_frame, I=I_frame)

# Add 4node elements
E_4node = 1e10
mu = 0.27
t = 2
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_4node.add_element(EN=[3, 4, 5, 6], E=E_4node, mu=mu, t=t)

# Add null-elements
cke = E_frame*A_frame
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[3, 0], cke=cke, alpha=math.pi)
element_null.add_element(EN=[5, 2], cke=cke, alpha=math.pi / 2)
# works as with 2 null elements normal in the corner
# element_null.add_element(EN=[6, 1], cke=cke, alpha=math.pi, add_t_el=False)
# element_null.add_element(EN=[6, 1], cke=cke, alpha=math.pi / 2, add_t_el=False)
# as below with 2 null elements normal and tangential in the corner
element_null.add_element(EN=[6, 1], cke=cke, alpha=math.pi * 3 / 4)

# Form stiffness matrix and load vector
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
# set stiffness springs
spring_stiffness = 1e9
sm.add_spring(degree_of_freedom=15, stiffness=spring_stiffness)
sm.add_spring(degree_of_freedom=16, stiffness=spring_stiffness)
sm.support_nodes([7, 8], direction='hv')
print('R rank:{}, R rows:{}'.format(np.linalg.matrix_rank(sm.r), sm.r.shape[0]))  # check if system is construction

lv_const = LoadVector()
lv_const.add_concentrated_force(force=-F, degree_of_freedom=11)
lv_const.add_concentrated_force(force=F, degree_of_freedom=12)
lv_variable = None


# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=0.1, autorun=True)

for i in range(len(graph.lemke.zn_anim)):
    print(f'{i}: p:{graph.lemke.p_anim[i]} zn:{graph.lemke.zn_anim[i]} xn:{graph.lemke.xn_anim[i]}'
          f'zt:{graph.lemke.zt_anim[i]} xt:{graph.lemke.xt_anim[i]}')
    print(f'sums: zn:{sum(graph.lemke.zn_anim[i])} xn:{sum(graph.lemke.xn_anim[i])}'
          f'zt:{sum(graph.lemke.zt_anim[i])} xt:{sum(graph.lemke.xt_anim[i])}')
# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)

