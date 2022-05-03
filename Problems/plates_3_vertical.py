import time
import numpy as np
from matplotlib import pyplot as plt
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector, solve_slae
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from LCP.initial_table import InitialTable  # to form initial table for LCP
from LCP.lemke import Lemke  # to solve LCP
import math

from Visualize.plot_data_scheme import PlotScheme  # for visualizing
from GUI.tkinter_gui import ContactFEM
from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

# elements variables
E_plate = 3.5e10  # Pa
mu_plate = 0.27  #
t_plate = 1  # m

plate_height = 8
plate_length = 2
mesh_size = 0.5 #0.25
gap_left = 0.3 #0.01
gap_right = 0.3 #0.01
F = 1e9 # 1e8


start = time.time()

nodes = NodeContainer()
# add nodes for ME 1 left
nodes.add_node(0, 0)
nodes.add_node(plate_length, 0)
nodes.add_node(plate_length, plate_height)
nodes.add_node(0, plate_height)
# add nodes for ME 2 central
nodes.add_node(0 + plate_length + gap_left, 0)
nodes.add_node(plate_length + plate_length + gap_left, 0)
nodes.add_node(plate_length + plate_length + gap_left, plate_height)
nodes.add_node(0 + plate_length + gap_left, plate_height)
# add nodes for ME 3 right
nodes.add_node(plate_length * 2 + gap_left + gap_right, 0)
nodes.add_node(plate_length + plate_length * 2 + gap_left + gap_right, 0)
nodes.add_node(plate_length + plate_length * 2 + gap_left + gap_right, plate_height)
nodes.add_node(plate_length * 2 + gap_left + gap_right, plate_height)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = None
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=False)
element_macro.add_element(EN=[2*4, 2*4 + 1, 2*4 + 2, 2*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=False)

element_macro.fragment_all(element_4node, element_frame, element_null)

# add null elements
nodes_left_contact_1 = nodes.find_nodes_numbers_along_segment((plate_length, plate_height), (plate_length, mesh_size/2))
nodes_right_contact_1 = nodes.find_nodes_numbers_along_segment((plate_length + gap_left, plate_height), (plate_length + gap_left, mesh_size/2))
nodes_left_contact_2 = nodes.find_nodes_numbers_along_segment((plate_length*2 + gap_left, plate_height), (plate_length*2 + gap_left, mesh_size/2))
nodes_right_contact_2 = nodes.find_nodes_numbers_along_segment((plate_length*2 + gap_left + gap_right, plate_height), (plate_length*2 + gap_left + gap_right, mesh_size/2))
# contact zone 1 between central plate and left plate
if gap_left == 0: # if there is no gap and nodes in contact zone all in 1 list
    for i in range(0, len(nodes_left_contact_1), 2):
        node_num1 = nodes_left_contact_1[i]
        node_num2 = nodes_left_contact_1[i+1]
        element_null.add_element(EN=[node_num1, node_num2], cke=E_plate, alpha=0)
else:
    for i, j in zip(nodes_left_contact_1, nodes_right_contact_1):
        element_null.add_element(EN=[i, j], cke=E_plate, alpha=0)
# contact zone 2 between central plate and right plate
if gap_right == 0: # if there is no gap and nodes in contact zone all in 1 list
    for i in range(0, len(nodes_left_contact_2), 2):
        node_num1 = nodes_left_contact_2[i]
        node_num2 = nodes_left_contact_2[i+1]
        element_null.add_element(EN=[node_num1, node_num2], cke=E_plate, alpha=0)
else:
    for i, j in zip(nodes_left_contact_2, nodes_right_contact_2):
        element_null.add_element(EN=[i, j], cke=E_plate, alpha=0)




# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_support = nodes.find_nodes_numbers_along_segment((0, 0), (plate_length * 3 + gap_left + gap_right, 0))
sm.support_nodes(nodes_to_support, direction='hv')
nodes_left_side = nodes.find_nodes_numbers_along_segment((0, 0), (0, plate_height))
lv = LoadVector()
lv.add_concentrated_force(F, nodes_left_side[-1] * 2)
# lv.add_concentrated_force(F, nodes_left_side[-2] * 2)
# lv.add_concentrated_force(F, nodes_left_side[-3] * 2)

u_linear = solve_slae(sm, lv)

# set to show only first 5 numbers when printing numpy values
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# do stuff about contact SM and LV
intl_table = InitialTable(element_null, sm, lv, u_linear)
intl_table.form_initial_table()
# do lemke ot solve LCP
lemke = Lemke(intl_table)
lemke.lcp_solve()

# # plot info elements
# for i, el in enumerate(element_frame):
#     print(f'frame {i} EN:{el.EN}, MI:{el.MI}')
# for i, el in enumerate(element_4node):
#     print(f'4node {i} EN:{el.EN}, MI:{el.MI}')
# for i, el in enumerate(element_null):
#     print(f'null {i} EN:{el.EN}, MI:{el.MI}, alpha:{el.alpha}')
# Stiffness matrix info
# smt_sub_sm = sm.r.T - sm.r
# print(f'SM rank:{np.linalg.matrix_rank(sm.r)}, rows:{sm.r.shape[1]}, cond:{np.linalg.cond(sm.r)} smT-sm max:{np.max(smt_sub_sm)}, min:{np.min(smt_sub_sm)}')



# calculate time
end = time.time()
last = end - start
print("Time: ", last)

# plot -------------------------------------------------------------------------- calculate data to plot
graph = PlotScheme(nodes, element_null, sm, lv, u_linear, lemke,
                   element_container_obj=element_4node, element_frame=element_frame, partition=10, scale_def=1, text=False)
u_contact = graph.u_contact_anim[-1]

# Contact info
print(f'zn: {lemke.zn}\nzt:{lemke.zt}')
print(f'xn: {lemke.xn}\nxt:{lemke.xt}')
print(f'min zn:{min(lemke.zn)}, max zn: {max(lemke.zn)}, min zt:{min(lemke.zt)}, max zt:{max(lemke.zt)}')
print(f'min xn:{min(lemke.xn)}, max xn: {max(lemke.xn)}, min xt:{min(lemke.xt)}, max xt:{max(lemke.xt)}')
print(f'sum xn:{sum(lemke.xn)}, sum xt{sum(lemke.xt)}')
print(f'max_u_linear: {max(u_linear)}, min_u_linear: {min(u_linear)}')
print(f'max_u_contact: {min(u_contact)}, min_u_contact: {min(u_contact)}')
print(f'top left node displacement horizontal: {u_contact[nodes_left_side[-1] * 2]}')
print(f'FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP:{FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP}')
if gap_left != 0 or gap_right != 0:
    print('Scale def must be = 1!')
if mesh_size == 0.25:
    print(f'Displacements for top right node x:{u_contact[890 * 2]}, y:{u_contact[890 * 2 + 1]}')

app = ContactFEM(graph=graph)
app.mainloop()
#
# ps = PlotScheme(nodes, element_null, sm=None, lv=None, u_linear=None, lemke=None, element_container_obj=element_4node,
#                 element_frame=element_frame, partition=10, scale_def=1)

# # Scheme with nodes numbers
# fig1, ax1 = plt.subplots(1, 1)
# plt.gca().set_aspect('equal', adjustable='box')
# fig1.patch.set_facecolor('darkgray')
# ax1.patch.set_facecolor('darkgray')
# for i, node in enumerate(nodes.nodes_list):
#     ax1.scatter(node.x, node.y, color='lightgreen', zorder=3, marker='s')
# ps.plot_list_of_elements(ax1)
# for i, node in enumerate(nodes):
#     dx = 0.2
#     dy = 0.2
#     if i < (plate_length/mesh_size + 1) * (plate_height/mesh_size + 1):
#         ax1.text(node.x + dx / 5, node.y - dy / 1, str("%.0f" % i), color='brown')
#     elif i < (plate_length/mesh_size + 1) * (plate_height/mesh_size + 1)*2:
#         ax1.text(node.x + dx / 5, node.y + dy / 5, str("%.0f" % i), color='black')
#     else:
#         ax1.text(node.x + dx / 5, node.y - dy / 1, str("%.0f" % i), color='brown')

# # Scheme deformed
# fig2, ax2 = plt.subplots(1, 1)
# plt.gca().set_aspect('equal', adjustable='box')
# fig2.patch.set_facecolor('darkgray')
# ax2.patch.set_facecolor('darkgray')
# ps.plot_nodes(ax2, deformed=True, linear=False)
# ps.plot_def_i_step_lemke(plt_def_contact=ax2, i=lemke.steps)
#
# # NORMAL interaction forces and displacements


plt.show()
