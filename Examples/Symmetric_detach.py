import math
import numpy as np
import matplotlib as mpl
import time

from input_data import SCALE_DEF
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector, solve_slae
from FEM.element_null import ElementNullContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_frame import ElementFrameContainer
from LCP.initial_table import InitialTable
from LCP.lemke import Lemke
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

mpl.use('TkAgg')

start = time.time()

# set inputs
Ar = 0.2  # Cross-section area
Er = 2e11  # Yong's modulus
Ix = 2.3 * 10 ** (-4)  # Moment of inertia
F = 1e6  # Force external
# add nodes ---------------------------------------------------------------------
nodes = NodeContainer()
for i in range(6):  # adding nodes every 1 m (along x axis)
    nodes.add_node(i, 0)
    nodes.add_node(i, 0)
nodes.add_node(2.5, 0)  # add node in the center
# add frame elements -------------------------------------------------------------
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[0, 2], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[2, 4], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[4, 12], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[12, 6], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[6, 8], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[8, 10], E=Er, A=Ar, I=Ix)
# add 4 node elements
element_4node = None  # there are no such elements
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
for i in range(0, 12, 2):
    element_null.add_element(EN=[i + 1, i], cke=123, alpha=math.pi / 2)


# form R and solve SLAE --------------------------------------------------------
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=[], el_null=element_null)
nodes_to_support = [i for i in range(1, 12, 2)]
sm.support_nodes(nodes_to_support, direction='hv')
print('rank and shape0:', np.linalg.matrix_rank(sm.r), sm.r.shape[0])
print(f'cond: {np.linalg.cond(sm.r)}')
#assert np.linalg.matrix_rank(sm.r) == sm.r.shape[0]  # check if system in construction

lv = LoadVector()
lv.add_concentrated_force(force=-F, degree_of_freedom=25)
#lv.add_concentrated_force(force=-F/100, degree_of_freedom=21)

u_linear = solve_slae(sm, lv)

# set to show only first 5 numbers when printing numpy values
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})


# do stuff about contact SM and LV ---------------------------------------------
intl_table = InitialTable(element_null, sm, lv, u_linear)
intl_table.form_initial_table()
# do lemke / solve LCP
lemke = Lemke(intl_table)
lemke.lcp_solve()

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, element_frame=element_frame, element_container_obj=None,
                   element_null=element_null, partition=10, scale_def=1000)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
