import time
import numpy as np
import math
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from FEM.element_4node import Element4NodeLinear

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application
from prettytable import PrettyTable



# elements variables
E_soil = 4.5e7  # Pa
mu_soil = 0.27  #
t_soil = 1  # m
own_weight_soil = 22661.1  # N/m^3
# own_weight_soil = 0  # N/m^3
E_pile = 2e11  # Pa
A_pile = 1.831e-2  # m^2
I_pile = 2.760e-5  # m^4
# scheme data
area_width = 20
area_height = 5
mesh_size = 0.25

# force
F = 1e3  # Н
F_const = 2.135184824469461 * F
# PLANE_STRAIN = True
# FRICTION_COEFFICIENT = 0.19
# ACCURACY_OF_LCP = 10**(-14)


start = time.time()

nodes = NodeContainer()

# add nodes for soil ME 1 left top
nodes.add_node(0, 0)
nodes.add_node(area_width/2, 0)
nodes.add_node(area_width / 2, area_height)
nodes.add_node(0, area_height)
# nodes for soil ME 2 right top
for node_num in range(4):
    nodes.add_node(nodes[node_num].x + area_width/2, nodes[node_num].y)
# add nodes for frame (pile)
first_node_num_pile = len(nodes)
for i in range(int(area_height / mesh_size) + 1):  # number of nodes for frame
    nodes.add_node(area_width/2, i * mesh_size)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)
# add frame elements
for i in range(first_node_num_pile, first_node_num_pile + int(area_height / mesh_size)):
    element_frame.add_element(EN=[i, i + 1], E=E_pile, A=A_pile, I=I_pile)
# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.fragment_all(element_4node, element_frame, element_null)

# nodes alongside pile
nodes_along_pile = nodes.find_nodes_numbers_along_segment(point1=(area_width/2, 0), point2=(area_width/2, area_height),
                                                          sorted_by_x=False)
list_nodes_left_soil = []
list_nodes_right_soil = []
list_nodes_pile = list(np.arange(1, int(area_height/mesh_size)+1))
amount_of_nodes_left_area = int((area_width / (mesh_size*2) + 1) * (area_height / mesh_size + 1)) + len(list_nodes_pile)
for i, node_num in enumerate(nodes_along_pile[3:]):
    a = nodes[node_num].elements[0]
    if isinstance(nodes[node_num].elements[0], Element4NodeLinear):
        if node_num <= amount_of_nodes_left_area:
            list_nodes_left_soil.append(node_num)
        else:
            list_nodes_right_soil.append(node_num)

# adding null elements
for pile_node, soil_left_node in zip(list_nodes_pile, list_nodes_left_soil):
    element_null.add_element(EN=[pile_node, soil_left_node], cke=E_pile * A_pile, alpha=math.pi, gap_length=0)
for pile_node, soil_right_node in zip(list_nodes_pile, list_nodes_right_soil):
    element_null.add_element(EN=[pile_node, soil_right_node], cke=1, alpha=0, gap_length=0)
#element_null.add_element(EN=[1, 2], cke=1, alpha=0)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(area_width, 0))
nodes_to_sup_sides = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(0, area_height))
nodes_to_sup_sides += (nodes.find_nodes_numbers_along_segment(point1=(area_width, 0),
                                                              point2=(area_width, area_height)))
sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(nodes_to_sup_sides, direction='h')

tip_pile_node = int(area_height / mesh_size)  # int(h_pile / mesh_size)+1 - its top node of the pile
lv_const = LoadVector()
# lv_const.add_concentrated_force(F_const, tip_pile_node * 2)
# lv_const.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])
# lv_const.add_concentrated_force(-F, tip_pile_node * 2 + 1)
lv_variable = LoadVector()
lv_variable.add_concentrated_force(F, tip_pile_node * 2)
# lv_variable.add_concentrated_force(-F, 130*2+1, vector_num=1)


# set to show only first 5 numbers when printing numpy values
np.set_printoptions(formatter={'float': lambda x: "{0:0.9f}".format(x)})
# calculate time
end = time.time()
last = end - start
print("Time scheme form: ", last)

# plot --------------------------------------------------------------------------
# Calculation and plotting object
autorun = True
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=2e3, autorun=autorun)

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


# коэф трения 0.6 плоская деформация есть зона слипания, проскальзывания
# ACCURACY_OF_LCP = 10**(-4)
# variable
# 88.026041068 366.082349865 843.283446894 950.202622545 xn top right
# const
# 88.026041068 366.082349865 843.283446893 950.202622545 xn top right