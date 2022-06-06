import time
import numpy as np
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application
from prettytable import PrettyTable

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
mesh_size = 1
# force
F = 100000  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! проверить
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
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)

element_macro.fragment_all(element_4node, element_frame, element_null)

# adding null elements\
amount_of_nodes_in_frame = pile_len / mesh_size + 1
amount_of_nodes_in_left = (area_len / mesh_size + 1) ** 2
f_l = amount_of_nodes_in_frame + amount_of_nodes_in_left
first_node_null_el_left = int(pile_len / mesh_size + (area_len / mesh_size + 1) ** 2)
first_node_null_el_right = int(pile_len / mesh_size + (area_len / mesh_size + 1) ** 2 * 2 - (area_len / mesh_size))
for i in range(int(pile_len / mesh_size) + 1):  # -1 to make without last node
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
sm.support_nodes([7], direction='hvr')
lv_const = LoadVector()
# lv_const.add_concentrated_force(F/1000, 0)
lv_const.add_concentrated_force(-F, 1)
lv_const.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])

# set to show only first 5 numbers when printing numpy values
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# calculate time
end = time.time()
last = end - start
print("Time: ", last)
print(len(sm.r), np.linalg.matrix_rank(sm.r))

# plot --------------------------------------------------------------------------
# Calculation and plotting object
autorun = True
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=40, autorun=autorun)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'var'] + [i for i in range(len(graph.lemke.zn_anim[0]))]
    for i in range(len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], 'zn'] + graph.lemke.zn_anim[i].tolist())
        mytable.add_row([i, graph.lemke.p_anim[i], 'zt'] + graph.lemke.zt_anim[i].tolist())
        mytable.add_row([i, graph.lemke.p_anim[i], 'xn'] + graph.lemke.xn_anim[i].tolist())
        mytable.add_row([i, graph.lemke.p_anim[i], 'xt'] + graph.lemke.xt_anim[i].tolist())
        print(i)
    print(mytable)


if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
