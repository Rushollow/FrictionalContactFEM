import math
import numpy as np
import time
import pyqtgraph as pg
from prettytable import PrettyTable

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application

start = time.time()
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# set inputs

q = 3975  # N/m Uniformly Distributed Load
general_length = 260  # meter
n = 5  # amount of nodes of frame MINIMUM 2
Ar = math.pi / 2 * (1.5 ** 2 - (1.5 - 0.02) ** 2)
Er = 1.95e9  # N/m
Ix = math.pi * 1.5 ** 2 * 0.02 / 8  #
F = q * general_length / (n - 1)  # concentrated force in each node
Fv = 243.75e3  # N (force pulling to the right)



nodes = NodeContainer()
for i in range(n):  # add nodes for frame
    nodes.add_node(i * general_length / (n - 1), 0)
nodes_to_sup = []
for i in range(n - 1):  # add nodes for supports
    nodes.add_node((i + 1) * general_length / (n - 1), 0)
    nodes_to_sup.append(n + i)

# add elements
element_4node = None
# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
for i in range(n - 1):
    element_frame.add_element(EN=[i, i + 1], E=Er, A=Ar, I=Ix)
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
for i, j in enumerate(range(n, len(nodes))):
    element_null.add_element(EN=[j, i + 1], cke=1, alpha=math.pi / 2, add_t_el=True)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.support_nodes(list_of_nodes=nodes_to_sup, direction='hv')  # sup for unilateral
sm.support_nodes(list_of_nodes=[0], direction='hvr')  # sup for zero left node

lv_const = LoadVector()
for node_num in range(1, n - 1):
    lv_const.add_concentrated_force(force=-F, degree_of_freedom=node_num * 2 + 1)
lv_const.add_concentrated_force(force=-F / 2, degree_of_freedom=(n - 1) * 2 + 1)  # last right node (half-length q)
lv_const.add_concentrated_force(force=Fv, degree_of_freedom=(n - 1) * 2)  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lv_variable = LoadVector()
lv_variable.add_concentrated_force(force=Fv, degree_of_freedom=(n - 1) * 2)

# create data for analytical solution
from input_data import FRICTION_COEFFICIENT
Fn = q
def u_func(z_list, lc):
    res = []
    for z in z_list:
        if z > lc:
            res.append(FRICTION_COEFFICIENT * Fn * (z-lc)**2 / (2 * Er*Ar))
        else:
            res.append(0.0)
    return res
def N_func(z_list, lc):
    res = []
    for z in z_list:
        if z > lc:
            res.append(FRICTION_COEFFICIENT * Fn * (z-lc))
        else:
            res.append(0.0)
    return res
lf = Fv / (FRICTION_COEFFICIENT * Fn)  # slippage distance length
lc = general_length - lf

z_list = list(np.linspace(0, general_length, n, endpoint=True))
y_analytical = u_func(z_list, lc)[1:]
y_numerical = []



autorun = True
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv_const, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=200, autorun=autorun)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    # print(mytable)

    # for i in range(len(graph.lemke.zn_anim)):
    #     print(list(graph.lemke.zt_anim[i]))
    y_numerical = list(graph.lemke.zt_anim[-1])
    print(f'Analytical: {y_analytical}, lc:{lc}, lf{lf}, L:{general_length}')
    print("Numerical: ", y_numerical)
    err = 0
    ut, uc = 0, 0
    for i, j in zip(y_analytical, y_numerical):
        ut += i**2
        uc += j**2
    print("Err:", (np.sqrt(ut)-np.sqrt(uc))/np.sqrt(ut)*100, "%")
# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
