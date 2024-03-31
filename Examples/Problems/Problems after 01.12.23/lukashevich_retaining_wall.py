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

from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP
assert FRICTION_COEFFICIENT == 1, 'Friction coef need to be 1'
assert PLANE_STRAIN is True, 'PLANE STRAIN need to be true!'
assert ACCURACY_OF_LCP >= 1e-6

start = time.time()

# set inputs
t = 1  # Thickness
Erw = 26500e6  # бетон stiffness retaining wall
mu_rw = 0.15  # mu for retaining wall
trw = t   # thickness
gamma_rw = 24e3  # Н/m^3  own weight for retaining wall
Eg_bot = 2e10  # Spring stiffness 1!!!!
mu_g_bot = 0.2
gamma_g_bot = 0
qn = 263e3  #263e3
qt = 213e3  # 213e3
qx1 = 68e3
qx2 = 89e3
qx3 = 142e3
qx4 = 150e3
qy = 45e3
qgr1 = 471e3  # 471e3
qgr2 = 462e3
h0, h1, h2, h3 = 5, 1, 1, 5.5
L1, L2, L3 = 2, 2, 4.5
L2_1 = 1.5
L0 = L3
L4 = L0
mesh_size = 0.2  # 0.2
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
nodes.add_node(L0, 0)  # 12
nodes.add_node(L0+L1, 0)  # 13
nodes.add_node(L0+L1+L2, 0)  # 14
nodes.add_node(L0+L1+L2+L3, 0)  # 15
nodes.add_node(L0+L1+L2+L3+L4, h0)  # 16
nodes.add_node(L0+L1+L2+L3+L4, 0)  # 17

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
element_macro.add_element(EN=[10, 11, 0, 12], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L0/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=t, own_weight=gamma_g_bot, stitch=False)  # 4
element_macro.add_element(EN=[12, 0, 3, 13], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L1/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=t, own_weight=gamma_g_bot, stitch=False, stitch_list=[4])  # 5
element_macro.add_element(EN=[13, 3, 5, 14], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L2/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=t, own_weight=gamma_g_bot, stitch=False, stitch_list=[5])  # 6
element_macro.add_element(EN=[14, 5, 7, 15], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L3/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=t, own_weight=gamma_g_bot, stitch=False, stitch_list=[6])  # 7
element_macro.add_element(EN=[15, 7, 16, 17], frag_amount_h=int(h0/mesh_size/2), frag_amount_v=int(L4/mesh_size),
                          E=Eg_bot, mu=mu_g_bot, t=t, own_weight=gamma_g_bot, stitch=False, stitch_list=[7])  # 8


element_macro.fragment_all(element_4node, element_frame, element_null)
# find nodes for contact pairs
n_contact1 = nodes.find_nodes_numbers_along_segment((L0, h0), (L0+L1+L2+L3, h0))

# n null elements and adding t null elements silently
for i in range(0, len(n_contact1), 2):
    element_null.add_element(EN=[n_contact1[i+1], n_contact1[i]], cke=123, alpha=math.pi/2, gap_length=0)
    # print(f'contact1 at {n_contact1[i+1], n_contact1[i]}')

side1 = nodes.find_nodes_numbers_along_segment((L0+L1, h0+h1+h2+h3), (L0+L1+L2_1, h0+h1+h2+h3), sorted_by_y=False)
side2 = nodes.find_nodes_numbers_along_segment((L0+L1+L2_1, h0+h1+h2+h3), (L0+L1+L2, h0+h1+h2), sorted_by_y=False)
side3 = nodes.find_nodes_numbers_along_segment((L0+L1+L2, h0+h1+h2), (L0+L1+L2+L3, h0+h1), sorted_by_y=False)
side4 = nodes.find_nodes_numbers_along_segment((L0+L1+L2+L3, h0+h1), (L0+L1+L2+L3, h0), sorted_by_y=True)
# first 2 nodes are with contact pair
side5 = nodes.find_nodes_numbers_along_segment((L0+L1+L2+L3, h0), (L0+L1+L2+L3+L4, h0), sorted_by_y=False)[2:]
print(1, side1)
print(2, side2)
print(3, side3)
print(4, side4)
print(5, side5)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L0+L1+L2+L3+L4, 0))
nodes_to_sup_right = nodes.find_nodes_numbers_along_segment(point1=(L0+L1+L2+L3+L4, 0),
                                                            point2=(L0+L1+L2+L3+L4, h0),
                                                            relative_tolerance=0.5)
nodes_to_sup_left = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(0, h0))

sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(nodes_to_sup_right, direction='h')
sm.support_nodes(nodes_to_sup_left, direction='h')

lv = LoadVector()
lv_v = LoadVector()

lv.add_own_weight_to_rf(nodes, [element_4node])
print('SIDE 1:')
# add load at the top node SIDE1
for i, nn in enumerate(side1):
    parts = len(side1)-1
    length = L2_1 / parts
    force_v = -(L2_1*qn*(parts - 2*i))/(parts*parts)
    force_h = -qt * length
    if i == 0 or i == parts:
        force_v /= 2
        force_h /= 2
    lv.add_concentrated_force(force_v, degree_of_freedom=nn*2+1)
    lv.add_concentrated_force(-qt * length, degree_of_freedom=nn * 2)
    # print(f'{i=}|{force_v=},{length=}, {nn=}, dof={nn*2+1}')
    # print(f'{force_h=}, {length=} {nn=}, dof={nn * 2}')

print('SIDE 2:')
# SIDE 2 top right nodes
for i, nn in enumerate(side2):
    parts = len(side2) - 1
    length = h3 / parts
    force_h = -(h3*(i*qx2-i*qx1+parts*qx1))/(parts*parts)
    force_v = -qy * length
    if i == 0 or i == parts:
        force_h /= 2
        force_v /= 2
    lv.add_concentrated_force(force_v, degree_of_freedom=nn * 2 + 1)
    lv.add_concentrated_force(force_h, degree_of_freedom=nn * 2)
    # print(f'vertical   force={force_v}, {length=} {nn=}, dof={nn*2+1}')
    # print(f'horizontal force={force_h}, {length=} {nn=}, dof={nn * 2}')

print('SIDE 3:')
# SIDE 3 right top
for i, nn in enumerate(side3):
    parts = len(side3)-1
    length = L3 / parts
    force_h = -qx3 * h2/parts
    force_v = -qgr1 * length
    if i == 0 or i == parts:
        force_h /= 2
        force_v /= 2
    lv.add_concentrated_force(force_v, degree_of_freedom=nn * 2+1)
    lv.add_concentrated_force(force_h, degree_of_freedom=nn * 2)
    # print(f'{force_v=}, {length=} {nn=}, dof={nn*2+1}')
    # print(f'{force_h=}, length={h2/parts} {nn=}, dof={nn * 2}')
#
print('SIDE 4:')
# SIDE 4 right top
for i, nn in enumerate(side4):
    parts = len(side4) -1
    length = h1 / parts
    force_h = -qx4 * length
    if i == 0 or i == parts:
        force_h /= 2
    lv.add_concentrated_force(force_h, degree_of_freedom=nn * 2)
    # print(f'horizontal force={-qx4 * length}, {length=} {nn=}, dof={nn * 2-1}')
print('SIDE 5:')
for i, nn in enumerate(side5):
    parts = len(side5) - 1
    length = L4 / parts
    force_v = -qgr2 * length
    if i == parts:
        force_v /= 2
    lv.add_concentrated_force(force_v, degree_of_freedom=nn*2+1)
    # print(f'{force_v=}, {length=} {nn=}, dof={nn*2}')


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
                   partition=10, scale_def=1000, autorun=autorun)

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
