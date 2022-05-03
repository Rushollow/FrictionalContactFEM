import math
import numpy as np
import time

from FEM.element_frame import ElementFrameContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_null import ElementNullContainer
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector, solve_slae
from LCP.initial_table import InitialTable
from LCP.lemke import Lemke
from Visualize.plot_data_scheme import PlotScheme
from GUI.tkinter_gui import ContactFEM

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
#assert np.linalg.matrix_rank(sm.r) == sm.r.shape[0]

lv = LoadVector()
lv.add_concentrated_force(force=-F, degree_of_freedom=11)
lv.add_concentrated_force(force=F, degree_of_freedom=12)

u_linear = solve_slae(sm, lv)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})
for i, k in zip(u_linear, range(len(u_linear))):
    print(k, i)

intl_table = InitialTable(element_null, sm, lv, u_linear)
intl_table.form_initial_table()
# do lemke or solve LCP
lemke = Lemke(intl_table)
lemke.lcp_solve()


print('zn={}\nxn={} - normalized xn/F\nrF={} - normalized'.format(lemke.zn, lemke.xn/F, F/F))

# plot --------------------------------------------------------------------------
# calculate data to plot
graph = PlotScheme(nodes, element_null, sm, lv, u_linear, lemke,
                   element_container_obj=element_4node, element_frame=element_frame, scale_def=1)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

app = ContactFEM(graph=graph)
#app.geometry('1280x720')
app.mainloop()
