import math
import numpy as np
import time
from prettytable import PrettyTable

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

from input_data import FRICTION_COEFFICIENT
assert FRICTION_COEFFICIENT == 0.3, 'Friction coef need to be 0.3'

start = time.time()

# set inputs

Eg = 2e5
mu_g = 0.3
tg = 1
Erw = 2e8
mu_rw = 0.2
trw = 1
gamma_rw = 100  # Н/m^3
E_g_bot = 2000 # Spring stiffness 1!!!!
qn = 263
qt = 213
qx1 = 68
qx2 = 89
qx3 = 150
qx4 = 142
qy = 45
qgr1 = 471
qgr2 = 462
h0, h1, h2, h3 = 5, 1, 1, 5.5
L1, L2, L3 = 2, 2, 4.5
L2_1 = 1.5
mesh_size = 0.2
force_inc = False
autorun = False

# add nodes # for 4 node element
nodes = NodeContainer()
# add ME nodes
nodes.add_node(0, h0)  # 0
nodes.add_node(0, h0+h1)  # 1
nodes.add_node(L1, h0+h1+h2)  # 2
nodes.add_node(L1, h0)  # 3
nodes.add_node(L1+L2, h0+h1+h2)  # 4
nodes.add_node(L1+L2, h0)  # 5
nodes.add_node(L1+L2+L3, h0+h1)  # 6
nodes.add_node(L1+L2+L3, h0)  # 7
nodes.add_node(L1, h0+h1+h2+h3)  # 8
nodes.add_node(L1 + L2_1, h0+h1+h2+h3)  # 9
nodes.add_node(0, h0-0.1)  # 10 SUP NODE

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = None
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# add MEs
# ME retaining wall
element_macro.add_element(EN=[0, 1, 2, 3], frag_amount_h=int((h1+h2)/mesh_size), frag_amount_v=int(L1/mesh_size),
                          E=Erw, mu=mu_rw, t=trw, own_weight=gamma_rw, stitch=False)  # 0
element_macro.add_element(EN=[3, 2, 4, 5], frag_amount_h=int((h1+h2)/mesh_size), frag_amount_v=int(L2/mesh_size),
                          E=Erw, mu=mu_rw, t=trw, own_weight=gamma_rw, stitch=False, stitch_list=[0])  # 1
element_macro.add_element(EN=[5, 4, 6, 7], frag_amount_h=int((h1+h2)/mesh_size), frag_amount_v=int(L3/mesh_size),
                          E=Erw, mu=mu_rw, t=trw, own_weight=gamma_rw, stitch=False, stitch_list=[1])  # 2
element_macro.add_element(EN=[2, 8, 9, 4], frag_amount_h=int(h3/mesh_size), frag_amount_v=int(L2/mesh_size),
                          E=Erw, mu=mu_rw, t=trw, own_weight=gamma_rw, stitch=False, stitch_list=[0, 1, 2])  # 3

element_macro.fragment_all(element_4node, element_frame, element_null)
# find all nodes alongside bottom of the retaining wall
nodes_bot = nodes.find_nodes_numbers_along_segment(point1=(0, h0), point2=(L1+L2+L3, h0))
# n null elements and adding t null elements silently
for i in nodes_bot:
    element_null.add_element(EN=[0, i], cke=123, alpha=math.pi / 2, gap_length=0)

side1 = nodes.find_nodes_numbers_along_segment((L1, h0+h1+h2+h3), (L1+L2_1, h0+h1+h2+h3), sorted_by_y=False)
side2 = nodes.find_nodes_numbers_along_segment((L1+L2_1, h0+h1+h2+h3), (L1+L2, h0+h1+h2), sorted_by_y=False)
side3 = nodes.find_nodes_numbers_along_segment((L1+L2, h0+h1+h2), (L1+L2+L3, h0+h1), sorted_by_y=False)
side4 = nodes.find_nodes_numbers_along_segment((L1+L2+L3, h0+h1), (L1+L2+L3, h0), sorted_by_y=True)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.support_nodes([0], direction='hv')

for spring_node_num in nodes_bot:
    sm.add_spring(degree_of_freedom=spring_node_num*2, stiffness=E_g_bot)

lv = LoadVector()
lv_v = LoadVector()

for i, nn in enumerate(side1):
    length = L2_1 / len(side1)
    force = (2 * qn / L2_1 * length * i - qn) * length
    if i == 0 or i == len(side1)+1:
        force /= 2
    lv.add_concentrated_force(force, degree_of_freedom=nn*2)
    # print(f'{force=}, {length=} {nn=}')


if not force_inc:
    pass
    # lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])

else:
    pass
    # lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])


# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=autorun)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

# if autorun:
#     mytable = PrettyTable()
#     mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
#     for i in range(len(graph.lemke.zn_anim)):
#         mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
#                          graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
#     print(mytable)


if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
