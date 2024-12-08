
import math
import numpy as np
import time
from prettytable import PrettyTable
import pyqtgraph as pg

from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application
import xlsxwriter

from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

assert FRICTION_COEFFICIENT == 0.0, 'Friction coef need to be 0.5 НЕПРАВИЛЬНО!'
assert 1e-16 <= ACCURACY_OF_LCP <= 1e-15
print('Starting to calculate...')

start = time.time()

tonne_meters = True
# set inputs
# TONNES AND METERS
Er = 1  # tonne / meter^2
Ix = 44.5  # meter^4
L = 2  # meters
# load
F1 = - 0.70707  # tonne
F2 = 4.3597  # tonne
F3 = - 2.1155
if not tonne_meters:
    Er = 1
    Ix = 436.4  # kN*m^2
    F1 = -6.934  # kN
    F2 = 42.754 # kN
    F3 = -20.746  # kN

autorun = True
on_springs = False

nodes = NodeContainer()
# nodes for frame
nodes.add_node(0, 0)  # 0
nodes.add_node(L, 0)  # 1
nodes.add_node(L*2, 0)  # 2
nodes.add_node(L*3, 0)  # 3
# nodes for sups
nodes.add_node(0, 1)  # 4

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

element_frame.add_element(EN=[0, 1], E=Er, A=Ix*2, I=Ix)
element_frame.add_element(EN=[1, 2], E=Er, A=Ix*2, I=Ix)
element_frame.add_element(EN=[2, 3], E=Er, A=Ix*2, I=Ix)

# NULL ELEMENTS
element_null.add_element(EN=[4, 1], cke=Er, alpha=math.pi*3/2, gap_length=0)
element_null.add_element(EN=[4, 2], cke=Er, alpha=math.pi/2, gap_length=0)
element_null.add_element(EN=[4, 3], cke=Er, alpha=math.pi/2, gap_length=0)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.add_spring(degree_of_freedom=4*2+1, stiffness=1e1)
sm.add_spring(degree_of_freedom=4*2, stiffness=1e1)

sm.support_nodes([0], direction='hvr')
if not on_springs:
    print("Not springs!")
    sm.support_nodes([4], direction='hv')

# load cheme
lv = LoadVector()
lv_v = LoadVector()
lv.add_concentrated_force(force=F1, degree_of_freedom=1*2+1)
lv.add_concentrated_force(force=F2, degree_of_freedom=2*2+1)
lv.add_concentrated_force(force=F3, degree_of_freedom=3*2+1)

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=20, scale_def=10, autorun=autorun)

def maximum(vec: list):
    return max(abs(min(vec)), max(vec))

if autorun:
    xn = graph.lemke.xn_anim[-1]
    zn = graph.lemke.zn_anim[-1]
    for x, z in zip(xn, zn):
        print(f'xn: {x}, zn {z}')

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(0, len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    print(mytable)

    # for i, (p, xn, xt) in enumerate(zip(graph.lemke.p_anim, graph.lemke.xn_anim, graph.lemke.xt_anim)):
    #     print(f'{i}: sum Xn: {sum(xn)-p*len(xn)}, sum Xt: {sum(xt)}')

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
