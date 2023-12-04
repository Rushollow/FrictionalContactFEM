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
assert FRICTION_COEFFICIENT == 0.6, 'Friction coef need to be 0.6'
# elements variables
E_plate_top = 2e11  # Pa
E_plate_bot = 2e11  # Pa
mu_plate = 0.3  #
t_plate = 0.01  # m

plate_height = 0.1
plate_length = 0.2
mesh_size = 0.01

F = 1e5  # сила в Н
dead_weight_bot = 75500  # удельный вес в Н/м3
dead_weight_top = 75500

start = time.time()

nodes = NodeContainer()
# add nodes for ME 1 lower
nodes.add_node(0, 0)
nodes.add_node(plate_length, 0)
nodes.add_node(plate_length, plate_height)
nodes.add_node(0, plate_height)
# add nodes for ME 2 upper
nodes.add_node(0, 0 + plate_height)
nodes.add_node(plate_length, 0 + plate_height)
nodes.add_node(plate_length, plate_height + plate_height)
nodes.add_node(0, plate_height + plate_height)


# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = None
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_plate_bot, mu=mu_plate, t=t_plate,
                          own_weight=dead_weight_bot, stitch=False)  # ME bot
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_plate_top, mu=mu_plate, t=t_plate,
                          own_weight=dead_weight_top, stitch=False)  # Me top


element_macro.fragment_all(element_4node, element_frame, element_null)

# add null elements
nodes_bot_contact = []
nodes_top_contact = []
for i in range(int(plate_length / mesh_size) + 1):
    number = i + (plate_length / mesh_size + 1) * (plate_height / mesh_size + 1) - (plate_length / mesh_size + 1)
    nodes_bot_contact.append(int(number))
    nodes_top_contact.append(int(number) + int((plate_length / mesh_size + 1)))

print('top c_nodes:', nodes_bot_contact)
print('bot c_nodes:', nodes_top_contact)

# ТУТ ДОБАВЛЯЕМ НУЛЬ ЭЛЕМЕНТЫ ПО СЕТКЕ ЛИБО С ЗАДАННЫМ ЗАЗОРОМ
for node1, node2 in zip(nodes_bot_contact, nodes_top_contact):
    element_null.add_element(EN=[node1, node2], cke=E_plate_bot*t_plate, alpha=math.pi/2)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_side_bot = []

for i in range(int(plate_height / mesh_size) + 1):
    nodes_to_sup_side_bot.append(i*(int(plate_length / mesh_size) + 1))
    # nodes_to_sup_side_bot.append(i)
print(f'supported nodes: {nodes_to_sup_side_bot}')

sm.support_nodes(nodes_to_sup_side_bot, direction='hv')
# sm.support_nodes(nodes_to_sup_sides, direction='h')
lv = LoadVector()

node_bot_left = int(plate_length / mesh_size)
# нагрузка на правый нижний узел
lv.add_concentrated_force(F, node_bot_left * 2 + 1)  # нагрузка
lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])  # собственный вес

# Variable load
lv_var = None
# lv_var = LoadVector(vectors_amount=1)
# lv_var.add_concentrated_force(force=-F, degree_of_freedom=305*2+1)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

autorun = True

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_var,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=20, autorun=autorun, force_incrementation=False)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(len(graph.lemke.zn_anim)-2, len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    print(mytable)

print(f'xn sum: {sum(graph.lemke.xn)}, '
      f'xt sum: {sum(graph.lemke.xt)}')



# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)


