import math
import numpy as np
import time

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

# elements variables
E_plate = 3.5e10  # Pa
mu_plate = 0.2  #
t_plate = 0.1  # m

plate_height = 4
plate_length = 8
mesh_size = 0.5
gap_left = 0.0

F = 2.5e7


start = time.time()

nodes = NodeContainer()
# add nodes for ME 1 lower
nodes.add_node(0, 0)
nodes.add_node(plate_length, 0)
nodes.add_node(plate_length, plate_height)
nodes.add_node(0, plate_height)
# add nodes for ME 2 upper left
nodes.add_node(0, 0 + plate_height + gap_left)
nodes.add_node(plate_length/2, 0 + plate_height)
nodes.add_node(plate_length/2, plate_height + plate_height)
nodes.add_node(0, plate_height + plate_height)
# add nodes for ME 3 upper right
nodes.add_node(plate_length/2, 0 + plate_height)
nodes.add_node(plate_length, 0 + plate_height)
nodes.add_node(plate_length, plate_height + plate_height)
nodes.add_node(plate_length/2, plate_height + plate_height)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = None
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=False, frag_size_v=mesh_size) # frag_size_v=mesh_size-gap_left/plate_height
element_macro.add_element(EN=[2*4, 2*4 + 1, 2*4 + 2, 2*4 + 3], frag_size=mesh_size, E=E_plate, mu=mu_plate, t=t_plate,
                          own_weight=0, stitch=True, stitch_list=[1])

element_macro.fragment_all(element_4node, element_frame, element_null)

# add null elements
nodes_bot_contact = []
nodes_top_contact = []
for i in range(int(plate_length / mesh_size) + 1):
    number = i + (plate_length / mesh_size + 1) * (plate_height / mesh_size + 1) - (plate_length / mesh_size + 1)
    nodes_bot_contact.append(int(number))
    if i <= int((plate_length / mesh_size + 1) / 2):
        nodes_top_contact.append(int(number) + int((plate_length / mesh_size + 1)))
first_me3_node = int(((plate_length / mesh_size + 1) * (plate_height / mesh_size + 1)) +
          (plate_length / mesh_size + 2)/2 * (plate_height / mesh_size + 1))
for i in range(int((plate_length / mesh_size) / 2)):
    nodes_top_contact.append(first_me3_node + i)
print('top c_nodes:', nodes_bot_contact)
print('bot c_nodes:', nodes_top_contact)

# ТУТ ДОБАВЛЯЕМ НУЛЬ ЭЛЕМЕНТЫ ПО СЕТКЕ ЛИБО С ЗАДАННЫМ ЗАЗОРОМ
for node1, node2 in zip(nodes_bot_contact, nodes_top_contact):
    element_null.add_element(EN=[node1, node2], cke=E_plate*t_plate, alpha=math.pi/2)

# i = 0
# for node1, node2 in zip(nodes_bot_contact, nodes_top_contact):
#     gap_initial = 0.02
#     gap = gap_initial - gap_initial/8 * i
#     if gap < 0:
#         gap = 0
#     element_null.add_element(EN=[node1, node2], cke=E_plate * t_plate, alpha=math.pi / 2, gap_length=gap)
#     print(gap)
#     i += 1


# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = []
nodes_to_sup_sides = []  # ! nothing
for i in range(int(plate_length / mesh_size) + 1):
    nodes_to_sup_bot.append(i)
print(f'supported nodes: {nodes_to_sup_bot}')

sm.support_nodes(nodes_to_sup_bot, direction='hv')
# sm.support_nodes(nodes_to_sup_sides, direction='h')
first_node_top = int(nodes_top_contact[0]+(plate_length/(2*mesh_size)+1)*(plate_height/mesh_size))  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lv = LoadVector()
lv.add_concentrated_force(-F, first_node_top * 2 + 1)
lv.add_concentrated_force(-F, (first_node_top + 1) * 2 + 1)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, element_frame=element_frame, element_container_obj=element_4node,
                   element_null=element_null, partition=10, scale_def=1, autorun=True, force_incrementation=False)


print(f'xn sum: {sum(graph.lemke.xn)}, force sum = {2*F}'
      f'xt sum: {sum(graph.lemke.xt)}')



# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)

print('gap:', gap_left)
print('Верхний левый узел dx contact', (graph.u_contact_anim[-1])[first_node_top*2])
print('dy', (graph.u_contact_anim[-1])[153*2+1])





# gap: 0.02
# Верхний левый узел dx contact -0.02030243031744137
# dy -0.02006804729807805
# xt max: 8303.187780478103

# gap: 0.015
# Верхний левый узел dx contact -0.02032566925208392
# dy -0.020067642946017164
# xt max: 6225.86004079859

# gap: 0.01
# Верхний левый узел dx contact -0.02034880175982377
# dy -0.020067232507550024
# xt max: 4148.800545073441

# gap: 0.005
# Верхний левый узел dx contact -0.020371827497700116
# dy -0.0200668159904903
# xt max: 2072.0418311963085

# gap: 0.0
# Верхний левый узел dx contact -0.0203947461260796
# dy -0.02006639340283399
# xt max: 4.383578616183513