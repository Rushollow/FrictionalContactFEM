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

from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN
assert FRICTION_COEFFICIENT == 1, 'Friction coef need to be 1'
assert PLANE_STRAIN is True, 'PLANE STRAIN need to be true!'

start = time.time()

# set inputs

Eg = 2e5  # soil stiffness
mu_g = 0.3  # mu for soil
tg = 1  # thickness
Erw = 26500e6  # бетон stiffness retaining wall
mu_rw = 0.2  # mu for retaining wall
trw = 1   # thickness
gamma_rw = 24e3  # Н/m^3  own weight for retaining wall
E_g_bot = 20000e6  # Spring stiffness 1!!!!
qn = 263
qt = 213
qx1 = 68
qx2 = 89
qx3 = 142
qx4 = 150
qy = 45
qgr1 = 471
qgr2 = 462
h0, h1, h2, h3 = 5, 1, 1, 5.5
L1, L2, L3 = 2, 2, 4.5
L2_1 = 1.5
mesh_size = 0.2
force_inc = False
autorun = True

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

print('SIDE 1:')
# add load at the top node SIDE1
for i, nn in enumerate(side1):
    length = L2_1 / len(side1)
    force = (2 * qn / L2_1 * length * i - qn) * length
    if i == 0 or i == len(side1)+1:
        force /= 2
    lv.add_concentrated_force(force, degree_of_freedom=nn*2)
    lv.add_concentrated_force(-qt * length, degree_of_freedom=nn * 2 - 1)
    print(f'vertical {force=}, {length=} {nn=}, dof={nn*2}')
    print(f'horizontal force={-qt * length}, {length=} {nn=}, dof={nn * 2-1}')

print('SIDE 2:')
# SIDE 2 top right nodes
for i, nn in enumerate(side2):
    length = math.sqrt((L2-L2_1)*(L2-L2_1) + h3*h3) / len(side2)
    force_h = (2 * (qx1 + qx2) / math.sqrt((L2-L2_1)*(L2-L2_1) + h3*h3) * length * i - qx1) * length
    lv.add_concentrated_force(-qy * length, degree_of_freedom=nn * 2)
    lv.add_concentrated_force(-force_h, degree_of_freedom=nn * 2 - 1)
    print(f'vertical force={-qy * length}, {length=} {nn=}, dof={nn*2}')
    print(f'horizontal force={force_h}, {length=} {nn=}, dof={nn * 2-1}')

print('SIDE 3:')
# SIDE 3 right top
for i, nn in enumerate(side3):
    length = math.sqrt((h2*h2) + (L3*L3)) / len(side3)
    lv.add_concentrated_force(-qgr1 * length, degree_of_freedom=nn * 2)
    lv.add_concentrated_force(-qx3 * length, degree_of_freedom=nn * 2 - 1)
    print(f'vertical force={-qgr1 * length}, {length=} {nn=}, dof={nn*2}')
    print(f'horizontal force={-qx3 * length}, {length=} {nn=}, dof={nn * 2-1}')

print('SIDE 4:')
# SIDE 4 right top
for i, nn in enumerate(side4):
    length = h1 / len(side3)
    lv.add_concentrated_force(-qgr2 * length, degree_of_freedom=nn * 2)
    lv.add_concentrated_force(-qx4 * length, degree_of_freedom=nn * 2 - 1)
    print(f'vertical force={-qgr2 * length}, {length=} {nn=}, dof={nn*2}')
    print(f'horizontal force={-qx4 * length}, {length=} {nn=}, dof={nn * 2-1}')


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
                   partition=10, scale_def=200, autorun=autorun)

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
