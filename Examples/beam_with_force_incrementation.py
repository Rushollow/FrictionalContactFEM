import math
import numpy as np
import time

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

start = time.time()

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# add nodes # for 4 node element
nodes = NodeContainer()
length = 1  # meter
nodes_amount = 7
for i in range(nodes_amount):  # add nodes for frame
    nodes.add_node(i, 0)
# add nodes for one-way supports
nodes.add_node(2*length, 0)
nodes.add_node(3*length, 0)
nodes.add_node(4*length, 0)
nodes.add_node(6*length, 0)

# set inputs
Ar = 1
Er = 1
Ix = 1
E = 1
F = 0.1  # Newtons 0.1,

# add elements
element_4node = None
# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(nodes_amount-1):
    element_frame.add_element(EN=[i, i+1], E=Er, A=Ar, I=Ix)
    print(i, i+1)
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[7, 2], cke=123, alpha=math.pi/2, gap_length=0.1)
element_null.add_element(EN=[8, 3], cke=123, alpha=math.pi/2, gap_length=-0.1)
element_null.add_element(EN=[9, 4], cke=123, alpha=math.pi/2, gap_length=-0.1)
element_null.add_element(EN=[10, 6], cke=123, alpha=math.pi/2, gap_length=0)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_support = [0, 7, 8, 9, 10]
sm.support_nodes(nodes_to_support, direction='hv')
sm.support_nodes(list_of_nodes=[5], direction='v')
# add constant load
lv_const = LoadVector()
lv_const.add_concentrated_force(force=0, degree_of_freedom=3)
# add variable load, for 'force increment' algorithm, if lv_variable = None it means normal Lemke's algorithm
lv_variable = LoadVector(vectors_amount=2)
lv_variable.add_concentrated_force(force=-F, degree_of_freedom=3, vector_num=0)
lv_variable.add_concentrated_force(force=F, degree_of_freedom=3, vector_num=1)



# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=lv_variable,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=7, autorun=False)

# debug
if False:
    from LCP.lemke import Lemke
    from LCP.initial_table import InitialTable
    graph.intl_table = InitialTable(element_null=element_null, stiffness_matrix=sm, lv_const=lv_const, u_linear_const=None,
                 u_linear_variable=None, lv_variable=lv_variable)
    graph.intl_table.table = np.array([
        [1, 0, 0, -4548.5,  3979.9,  1137.1, -1, 0.375,  -1.125,   3.98],
        [0, 1, 0,  3979.9, -4548.5, -1705.7, -1,  -1.813, 5.438,  -4.549],
        [0, 0, 1,  1137.1, -1705.7, -758.08, -1, 1.0001, -3.0004, -1.70588],
    ], dtype=float)
    graph.intl_table.force_inc = True
    graph.intl_table.n_amount = 3
    graph.intl_table.t_amount = 0
    graph.intl_table.rf_const = np.array([1, 1])
    graph.lemke = Lemke(graph.intl_table)
    graph.lemke.lcp_solve()

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
