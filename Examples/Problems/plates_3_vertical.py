import time
import numpy as np
from prettytable import PrettyTable
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from LCP.initial_table import InitialTable  # to form initial table for LCP
from LCP.lemke import Lemke  # to solve LCP
import math

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application
from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

# elements variables
E_plate = 3.5e10  # Pa
mu_plate = 0.27  #
t_plate = 1  # m

plate_height = 8
plate_length = 2
mesh_size = 0.25  #0.25
gap_left = 0.3  #0.01
gap_right = 0.3  #0.01
F = 1e9  # 1e8


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
lv_const = LoadVector()
lv_const.add_concentrated_force(F/100, nodes_left_side[-1] * 2)
# lv_const.add_concentrated_force(F, nodes_left_side[-2] * 2)
# lv_const.add_concentrated_force(F, nodes_left_side[-3] * 2)
lv_variable = LoadVector()
lv_variable.add_concentrated_force(F, nodes_left_side[-1] * 2)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

autorun = True
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=1, autorun=autorun)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    # print(mytable)
# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
