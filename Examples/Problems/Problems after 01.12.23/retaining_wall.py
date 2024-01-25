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
gamma_g = 10000 # !!!!!!!!!!!!!!!!
Eg_bot = Eg*2  # !!!!!!!!!!!!
mu_g_bot = mu_g # !!!!!!!!!!!!
gamma_g_bot = gamma_g # !!!!!!!!!!!!!
Erw = 2e8
mu_rw = 0.2
trw = 1
gamma_rw = gamma_g # !!!!!!!!!!!!!!
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
L0 = L3
L4 = L0*2
h4 = 2
mesh_size = 0.3  # 0.2
force_inc = False
autorun = True

# add nodes # for 4 node element
nodes = NodeContainer()
# add ME nodes
nodes.add_node(L0+0, h0)  # 0
nodes.add_node(L0+0, h0+h1)  # 1
nodes.add_node(L0+L1, h0+h1+h2)  # 2
nodes.add_node(L0+L1, h0)  # 3
nodes.add_node(L0+L1+L2, h0+h1+h2)  # 4
nodes.add_node(L0+L1+L2, h0)  # 5
nodes.add_node(L0+L1+L2+L3, h0+h1)  # 6
nodes.add_node(L0+L1+L2+L3, h0)  # 7
nodes.add_node(L0+L1, h0+h1+h2+h3)  # 8
nodes.add_node(L0+L1 + L2_1, h0+h1+h2+h3)  # 9
nodes.add_node(0, 0)  # 10
nodes.add_node(0, h0)  # 11
nodes.add_node(L0+L1, 0)  # 12
nodes.add_node(L0+L1+L2, 0)  # 13
nodes.add_node(L0+L1+L2+L3, 0)  # 14
nodes.add_node(L0+L1+L2+L3+L4, h0)  # 15
nodes.add_node(L0+L1+L2+L3+L4, 0)  # 16
nodes.add_node(L0+L1+L2+L3, h0+h1+h2+h3+h4)  # 17
nodes.add_node(L0+L1+L2+L3+L4, h0+h1+h2+h3+2*h4)  # 18
nodes.add_node(L0+L1+L2+L3+L4, h0+h1)  # 19
nodes.add_node(L0, 0)  # 20


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
# МЕ soil
element_macro.add_element(EN=[20, 0, 3, 12], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L1/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=tg, own_weight=gamma_g_bot, stitch=False)  # 4
element_macro.add_element(EN=[12, 3, 5, 13], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L2/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=tg, own_weight=gamma_g_bot, stitch=False, stitch_list=[4])  # 5
element_macro.add_element(EN=[13, 5, 7, 14], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L3/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=tg, own_weight=gamma_g_bot, stitch=False, stitch_list=[5])  # 6
element_macro.add_element(EN=[14, 7, 15, 16], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L4/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=tg, own_weight=gamma_g_bot, stitch=False, stitch_list=[6])  # 7
element_macro.add_element(EN=[4, 9, 17, 6], frag_amount_h=int(h3/mesh_size), frag_amount_v=int(L3/mesh_size),
                          E=Eg, mu=mu_g, t=tg, own_weight=gamma_g, stitch=False)  # 8
element_macro.add_element(EN=[6, 17, 18, 19], frag_amount_h=int(h3/mesh_size), frag_amount_v=int(L4/mesh_size),
                          E=Eg, mu=mu_g, t=tg, own_weight=gamma_g, stitch=False, stitch_list=[8])  # 9
element_macro.add_element(EN=[7, 6, 19, 15], frag_amount_h=int((h1+h2)/mesh_size), frag_amount_v=int(L4/mesh_size),
                          E=Eg, mu=mu_g, t=tg, own_weight=gamma_g, stitch=False, stitch_list=[8, 7, 9, 6])  # 10
element_macro.add_element(EN=[10, 11, 0, 20], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L0/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=tg, own_weight=gamma_g_bot, stitch=False, stitch_list=[4])  # 11


element_macro.fragment_all(element_4node, element_frame, element_null)

# find nodes for contact pairs
n_contact1 = nodes.find_nodes_numbers_along_segment((L0, h0), (L0+L1+L2+L3, h0))
n_contact2 = nodes.find_nodes_numbers_along_segment((L0+L1+L2+L3, h0), (L0+L1+L2+L3, h0+h1))
n_contact3 = nodes.find_nodes_numbers_along_segment((L0+L1+L2+L3, h0+h1), (L0+L1+L2, h0+h1+h2),
                                                    sorted_by_x=False)
n_contact4 = nodes.find_nodes_numbers_along_segment((L0+L1+L2, h0+h1+h2), (L0+L1+L2_1, h0+h1+h2+h3),
                                                    sorted_by_x=False, sorted_by_y=True)
# add null elements to contact pairs
for i in range(0, len(n_contact1), 2):
    element_null.add_element(EN=[n_contact1[i+1], n_contact1[i]], cke=123, alpha=math.pi/2, gap_length=0)
    print(f'contact1 at {n_contact1[i+1], n_contact1[i]}')
angle4 = math.pi / 2 - math.atan(h3 / (L2-L2_1))
for i in range(2, len(n_contact4), 2):
    element_null.add_element(EN=[n_contact4[i], n_contact4[i+1]], cke=123, alpha=angle4, gap_length=0)
    print(f'contact4 at {n_contact4[i], n_contact4[i+1]}')
angle3 = math.pi/2 - math.atan(h2/L3)
for i in range(2, len(n_contact3), 2):
    element_null.add_element(EN=[n_contact3[i], n_contact3[i+1]], cke=123, alpha=angle3, gap_length=0)
    print(f'contact3 at {n_contact3[i], n_contact3[i+1]}')
for i in range(2, len(n_contact2), 2):
    element_null.add_element(EN=[n_contact2[i], n_contact2[i+1]], cke=123, alpha=0, gap_length=0)
    print(f'contact2 at {n_contact2[i], n_contact2[i+1]}')



# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L0+L1+L2+L3+L4, 0))
nodes_to_sup_right = nodes.find_nodes_numbers_along_segment(point1=(L0+L1+L2+L3+L4, 0),
                                                            point2=(L0+L1+L2+L3+L4, h0+h1+h2+h3+3*h4),
                                                            relative_toletance=0.5)
nodes_to_sup_left = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(0, h0))

print(nodes_to_sup_right)
sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(nodes_to_sup_right, direction='h')
sm.support_nodes(nodes_to_sup_left, direction='h')

lv = LoadVector()
lv_v = LoadVector()

if not force_inc:
    # lv.add_concentrated_force(force=qn, degree_of_freedom=0)
    lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])
else:
    lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])
    # lv_v.add_concentrated_force(force=qn, degree_of_freedom=0)


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
