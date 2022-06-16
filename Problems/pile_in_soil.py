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
area_width = 20
h_pile = 10
h_bot = 5
mesh_size = 1
rock_support = False

# force
F = 1e6  # Н
# PLANE_STRAIN = True
# FRICTION_COEFFICIENT = 0.19
# ACCURACY_OF_LCP = 10**(-14)


start = time.time()

nodes = NodeContainer()

# add nodes for soil ME 1 left top
nodes.add_node(0, h_bot)
nodes.add_node(area_width/2, h_bot)
nodes.add_node(area_width/2, h_bot+h_pile)
nodes.add_node(0, h_bot+h_pile)
# nodes for soil ME 2 right top
for node_num in range(4):
    nodes.add_node(nodes[node_num].x + area_width/2, nodes[node_num].y)
# nodes for soil ME 3 bottom
nodes.add_node(0, 0)
nodes.add_node(area_width, 0)
nodes.add_node(area_width, h_bot)
nodes.add_node(0, h_bot)
# add nodes for frame (pile)
first_node_num_pile = len(nodes)
for i in range(int(h_pile / mesh_size) + 1):  # number of nodes for frame
    nodes.add_node(area_width/2, h_bot + i * mesh_size)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)
# add frame elements
for i in range(first_node_num_pile, first_node_num_pile + int(h_pile / mesh_size)):
    element_frame.add_element(EN=[i, i + 1], E=E_pile, A=A_pile, I=I_pile)
# add MEs
element_macro.add_element(EN=[0*4, 0*4 + 1, 0*4 + 2, 0*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.add_element(EN=[1*4, 1*4 + 1, 1*4 + 2, 1*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=False)
element_macro.add_element(EN=[2*4, 2*4 + 1, 2*4 + 2, 2*4 + 3], frag_size=mesh_size, E=E_soil, mu=mu_soil, t=t_soil,
                          own_weight=own_weight_soil, stitch=True, stitch_list=[0, 1])
element_macro.fragment_all(element_4node, element_frame, element_null)

# adding null elements\
if rock_support:
    i = 0
    i += 1
element_null.add_element(EN=[1, 2], cke=E_pile * A_pile, alpha=0)





# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(area_width, 0))
nodes_to_sup_sides = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(0, h_bot + h_pile))
nodes_to_sup_sides += (nodes.find_nodes_numbers_along_segment(point1=(area_width, 0),
                                                              point2=(area_width, h_bot + h_pile)))

if rock_support:
    nodes_to_sup_rock = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot), point2=(area_width, h_bot))
    sm.support_nodes(nodes_to_sup_rock, direction='hvr')

sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(nodes_to_sup_sides, direction='h')

lv_const = LoadVector()
# lv_const.add_concentrated_force(F/1000, 0)
lv_const.add_concentrated_force(-F, (int(h_pile / mesh_size)+1)*2 + 1)  # int(h_pile / mesh_size)+1 - its top node of the pile
lv_const.add_concentrated_force(-F, (int(h_pile / mesh_size)+1)*2)
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
autorun = False
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=40, autorun=autorun)

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
