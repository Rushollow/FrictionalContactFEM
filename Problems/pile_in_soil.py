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
from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

from Visualize.plot_data_scheme import PlotScheme  # for visualizing
from GUI.tkinter_gui import ContactFEM

# elements variables
E_soil = 4.5e7  # Pa
mu_soil = 0.27  #
t_soil = 1  # m
own_weight_soil = 22661.1  # N/m^3
# own_weight_soil = 0  # N/m^3
E_pile = 2e10  # Pa
A_pile = 0.09  # m^2 a = 0.3м, b = 0.3м
I_pile = 0.3*(0.3**3)/12  # m^4
# scheme data
area_len = 10
pile_len = 7
mesh_size = 0.5
# force
F = 10000  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! проверить
# F = 1000000
# PLANE_STRAIN = True
# FRICTION_COEFFICIENT = 0.19
# ACCURACY_OF_LCP = 10**(-14)


start = time.time()

nodes = NodeContainer()

# add nodes for soil ME 1
nodes.add_node(0, 0)
nodes.add_node(area_len, 0)
nodes.add_node(area_len, area_len)
nodes.add_node(0, area_len)
# nodes for soil ME 2
for node_num in range(4):
    nodes.add_node(nodes[node_num].x + area_len, nodes[node_num].y)
first_node_num_pile = len(nodes)
# add nodes for frame (pile)
for i in range(int(pile_len / mesh_size) + 1):  # number of nodes for frame
    nodes.add_node(area_len, area_len - i * mesh_size)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)
# add frame elements
for i in range(first_node_num_pile, first_node_num_pile + int(pile_len / mesh_size)):
    element_frame.add_element(EN=[i, i + 1], E=E_pile, A=A_pile, I=I_pile)
# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=0.5, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=0.5, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)

element_macro.fragment_all(element_4node, element_frame, element_null)

# adding null elements
first_node_null_el_left = int(pile_len / mesh_size + (area_len / mesh_size + 1) ** 2)
first_node_null_el_right = int(pile_len / mesh_size + (area_len / mesh_size + 1) ** 2 * 2 - (area_len / mesh_size))
for i in range(int(pile_len / mesh_size) + 1):
    node1 = first_node_null_el_left - i * int(area_len / mesh_size + 1)
    element_null.add_element(EN=[node1, i], cke=E_pile*A_pile, alpha=0)  # adding left soil to pile
    node2 = first_node_null_el_right - i * int(area_len / mesh_size + 1)
    element_null.add_element(EN=[i, node2], cke=E_pile * A_pile, alpha=0)  # adding pile to right soil

start_number_left = int(pile_len / mesh_size + area_len / mesh_size + 1)
start_number_right = int(first_node_null_el_left + 1)
for i in range(int(pile_len / mesh_size) - 1):
    increment = i * (area_len / mesh_size + 1)
    end_increment = (area_len / mesh_size + 1) * 9
    node1 = int(start_number_left + end_increment - increment)
    node2 = int(start_number_right + end_increment - increment)
    element_null.add_element(EN=[node1, node2], cke=E_pile*A_pile, alpha=0)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = []
nodes_to_sup_sides = []
for i in range(int(area_len / mesh_size + 1)):
    nodes_to_sup_bot.append(i + int(pile_len / mesh_size + 1))  # nodes numbers bottom left ME
    nodes_to_sup_bot.append(first_node_null_el_left + 1 + i)  # nodes numbers bottom right ME
    nodes_to_sup_sides.append(int(pile_len / mesh_size + 1 + i * (area_len / mesh_size + 1)))
    nodes_to_sup_sides.append(int(first_node_null_el_left + (area_len / mesh_size + 1) + i * (area_len / mesh_size + 1)))

sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(nodes_to_sup_sides, direction='h')
lv = LoadVector()
lv.add_concentrated_force(F, 0)
lv.add_concentrated_force(-F, 1)
lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])

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
# Contact info
# print(f'zn: {lemke.zn}\nzt:{lemke.zt}')
# print(f'xn: {lemke.xn}\nxt:{lemke.xt}')
# print(f'sum xn:{sum(lemke.xn)}, sum xt{sum(lemke.xt)}')



# calculate time
end = time.time()
last = end - start
print("Time: ", last)

# plot --------------------------------------------------------------------------
# calculate data to plot
# graph = PlotScheme(nodes, element_null, sm, lv, u_linear, lemke,
#                    element_container_obj=element_4node, element_frame=element_frame, partition=10, scale_def=1)
#
# app = ContactFEM(graph=graph)
# app.mainloop()

ps = PlotScheme(nodes, element_null, sm, lv, u_linear=u_linear, lemke=lemke, element_container_obj=element_4node,
                element_frame=element_frame, partition=10, scale_def=1)

print(f'0:{ps.u_contact_anim[-1][0]}m\n'
      f'1:{ps.u_contact_anim[-1][1]}m\n'
      f'19:{ps.u_contact_anim[-1][19]}m\n'
      f'20:{ps.u_contact_anim[-1][20]}m\n')

# Contact info
print(f'zn: {lemke.zn}\nzt:{lemke.zt}')
print(f'xn: {lemke.xn}\nxt:{lemke.xt}')
print(f'min zn:{min(lemke.zn)}, max zn: {max(lemke.zn)}, min zt:{min(lemke.zt)}, max zt:{max(lemke.zt)}')
print(f'min xn:{min(lemke.xn)}, max xn: {max(lemke.xn)}, min xt:{min(lemke.xt)}, max xt:{max(lemke.xt)}')
print(f'sum xn:{sum(lemke.xn)}, sum xt{sum(lemke.xt)}')
print(f'max_u_linear: {max(u_linear)}, min_u_linear: {min(u_linear)}')
print(f'max_u_contact: {min(ps.u_contact_anim[-1])}, min_u_contact: {min(ps.u_contact_anim[-1])}')
print(f'FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP:{FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP}')


# Scheme with nodes numbers
# fig1, ax1 = plt.subplots(1, 1)
# plt.gca().set_aspect('equal', adjustable='box')
# for i, node in enumerate(nodes.nodes_list):
#     ax1.scatter(node.x, node.y, color='lightgreen', zorder=3, marker='s')
# ps.plot_list_of_elements(ax1)
# for i, node in enumerate(nodes):
#     dx = 0.2
#     dy = 0.2
#     if i <= pile_len/mesh_size:
#         ax1.text(node.x - dx / 1, node.y + dy / 5, str("%.0f" % i), color='blue')
#     elif i > pile_len/mesh_size and i < 452:
#         ax1.text(node.x - dx / 1, node.y - dy / 1, str("%.0f" % i), color='red')
#     else:
#         ax1.text(node.x + dx / 5, node.y + dy / 5, str("%.0f" % i), color='black')

# Scheme deformed
fig2, ax2 = plt.subplots(1, 1)
plt.gca().set_aspect('equal', adjustable='box')
ps.plot_nodes(ax2, deformed=True, linear=False)
ps.plot_def_i_step_lemke(plt_def_contact=ax2, i=lemke.steps)
# NORMAL interaction forces and displacements
fig3, ax3 = plt.subplots(2, 4)
ax3[0, 0].set_title('xn and zn')
ax3[1, 0].set_title('xt and zt')
i = 0
for xn, zn in zip(lemke.xn, lemke.zn):
    if i % 2 == 0:
        ax3[0, 1].barh(i, xn, color='red', height=0.8)
        ax3[0, 0].barh(i, zn, color='blue', height=0.8)
    else:
        ax3[0, 2].barh(i, xn, color='red', height=0.8)
        ax3[0, 3].barh(i, zn, color='blue', height=0.8)
    i += 1
ax3[0, 0].invert_yaxis()
ax3[0, 0].invert_xaxis()
ax3[0, 1].invert_yaxis()
ax3[0, 1].invert_xaxis()
ax3[0, 2].invert_yaxis()
ax3[0, 3].invert_yaxis()

# Tangent interaction forces and displacements
i = 0
for xt, zt in zip(lemke.xt, lemke.zt):
    if i % 2 == 0:
        ax3[1, 1].barh(i, xt, color='red', height=0.8)
        ax3[1, 0].barh(i, zt, color='blue', height=0.8)
    else:
        ax3[1, 2].barh(i, xt, color='red', height=0.8)
        ax3[1, 3].barh(i, zt, color='blue', height=0.8)
    i += 1
ax3[1, 0].invert_yaxis()
ax3[1, 0].invert_xaxis()
ax3[1, 1].invert_yaxis()
ax3[1, 1].invert_xaxis()
ax3[1, 2].invert_yaxis()
ax3[1, 3].invert_yaxis()

plt.show()
