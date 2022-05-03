import math
import numpy as np
import time

from input_data import SCALE_DEF
from FEM.element_frame import ElementFrameContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_null import ElementNullContainer
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector, solve_slae
from LCP.initial_table import InitialTable
from LCP.lemke import Lemke
from Visualize.plot_data_scheme import PlotScheme
from GUI.tkinter_gui import ContactFEM

start = time.time()  # to calculate time

nodes = NodeContainer()
for i in range(4):
    nodes.add_node(i, 3)
for i in range(3):
    nodes.add_node(3, 2-i)
for i in range(3):
    nodes.add_node(2-i, 0)
nodes.add_node(0, 1)
nodes.add_node(0, 2)
for i in range(3):
    for j in range(3):
        nodes.add_node(0.5+j, 2.5-i)

E_frame = 1e9  # Pa
A_frame = 0.2*0.2  # m^2
I_frame = 1*10**(-4)  # m^3
F = 1e6  # N
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(11):
    element_frame.add_element(EN=[i, i + 1], E=E_frame, A=A_frame, I=I_frame)
element_frame.add_element(EN=[11, 0], E=E_frame, A=A_frame, I=I_frame)
E_4node = 1e9
mu = 0.27
t = 2
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_4node.add_element(EN=[15, 16, 13, 12], E=E_4node, mu=mu, t=t)
element_4node.add_element(EN=[16, 17, 14, 13], E=E_4node, mu=mu, t=t)
element_4node.add_element(EN=[18, 19, 16, 15], E=E_4node, mu=mu, t=t)
element_4node.add_element(EN=[19, 20, 17, 16], E=E_4node, mu=mu, t=t)

gap = 0.5
cke = E_4node*t
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[18, 10], cke=cke, alpha=math.pi)  # link to the left
element_null.add_element(EN=[15, 11], cke=cke, alpha=math.pi)  # link to the left
element_null.add_element(EN=[12, 0], cke=cke, alpha=math.pi * 3 / 4)  # to the top left corner
element_null.add_element(EN=[13, 1], cke=cke, alpha=math.pi / 2)  # to the top
element_null.add_element(EN=[14, 2], cke=cke, alpha=math.pi / 2)  # to the top
print(nodes[18].x, nodes[18])
for i in [18, 10, 15, 11]:
    print('Node {}: x:{}, y:{}'.format(i, nodes[i].x, nodes[i].y))


sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
# No need nodes to support, springs will be added
# set stiffness of deflection wall, add springs
wall_spring_stiffness = 1e8
sm.add_spring(degree_of_freedom=0, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=1, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=6, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=7, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=12, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=13, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=18, stiffness=wall_spring_stiffness)
sm.add_spring(degree_of_freedom=19, stiffness=wall_spring_stiffness)
# set stiffness of rubber bearings
bearing_stiffness = 100
sm.add_spring(degree_of_freedom=32, stiffness=bearing_stiffness)
sm.add_spring(degree_of_freedom=33, stiffness=bearing_stiffness)

lv = LoadVector()
lv.add_concentrated_force(force=-F, degree_of_freedom=32)
lv.add_concentrated_force(force=F*0.95, degree_of_freedom=33)

u_linear = solve_slae(sm, lv)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})
for i, k in zip(u_linear, range(len(u_linear))):
    print(k, i)

intl_table = InitialTable(element_null, sm, lv, u_linear)
intl_table.form_initial_table()
# do lemke or solve LCP
lemke = Lemke(intl_table)
lemke.lcp_solve()

print('zn={}\nxn={}\nF={}'.format(lemke.zn_anim[-1], lemke.xn_anim[-1]/F, F/F))

# plot --------------------------------------------------------------------------
# calculate data to plot
graph = PlotScheme(nodes, element_null, sm, lv, u_linear, lemke,
                   element_container_obj=element_4node, element_frame=element_frame, partition=10, scale_def=1)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

app = ContactFEM(graph=graph)
#app.geometry('1280x720')
app.mainloop()
