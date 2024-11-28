
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

from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

assert FRICTION_COEFFICIENT == 0.2, 'Friction coef need to be 0.4 НЕПРАВИЛЬНО!'
assert PLANE_STRAIN is False, 'PLANE STRAIN need to be false! НЕПРАВИЛЬНО!!!!'
assert 1e-15 <= ACCURACY_OF_LCP <= 1e-10
print('Starting to calculate...')

start = time.time()

# set inputs
t = 0.1  # Thickness
E_bot = 2e10  # steel
E_top = 2e10
mu = 0.2
# load
F = 10_000_000 # Н
# sizes
L_bot = 4
L_top = L_bot / 2
h_bot = 0.5
h_top = 0.5
force_inc = False
autorun = True

mesh_size = 0.05

# add nodes  for frame element
nodes = NodeContainer()
# nodes for bot macro element
nodes.add_node(0, 0) # 0
nodes.add_node(L_bot, 0) # 1
nodes.add_node(L_bot, h_bot) # 2
nodes.add_node(0, h_bot) # 3
# nodes for top macro element
nodes.add_node(0, h_bot) # 4
nodes.add_node(L_top, h_bot) # 5
nodes.add_node(L_top, h_bot + h_top) # 6
nodes.add_node(0, h_bot + h_top) # 7

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# make grid
element_macro.add_element(EN=[0, 1, 2, 3], frag_size=mesh_size, E=E_bot, mu=mu, t=t,
                          own_weight=0, stitch=False)
element_macro.add_element(EN=[4, 5, 6, 7], frag_size=mesh_size, E=E_top, mu=mu, t=t,
                          own_weight=0, stitch=False)
element_macro.fragment_all(element_4node=element_4node, element_frame=element_frame, element_null=element_null)

# add null elements
contact_nodes = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot), point2=(L_top, h_bot))
for i in range(2, len(contact_nodes), 2):
    contact_pair = sorted(contact_nodes[i:i+2])
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=111, alpha=math.pi/2)
    # print(f'{contact_pair[0], contact_pair[1]}')

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)

sup_nodes_left = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot), point2=(0, h_bot + h_top))
sup_nodes_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L_bot, 0))
print(sup_nodes_bot)
print(sup_nodes_left)
# sm.support_nodes(sup_nodes_left, direction='hv')
sm.support_nodes(sup_nodes_bot, direction='hv')


# load cheme
lv = LoadVector()
lv_v = LoadVector()
force_node = int(L_bot/mesh_size + 1) * int(h_bot/mesh_size + 1) + int(L_top/mesh_size + 1) * int(h_top/mesh_size + 1) - 1
print(f'node number with force: {force_node}')

if not force_inc:
    lv.add_concentrated_force(force=-F/4, degree_of_freedom=force_node*2 + 1)
    lv.add_concentrated_force(force=-F/2, degree_of_freedom=(force_node-1) * 2 + 1)
    lv.add_concentrated_force(force=-F/4, degree_of_freedom=(force_node-2) * 2 + 1)
    pass

else:
    pass

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=2, scale_def=1, autorun=autorun)

def maximum(vec: list):
    return max(abs(min(vec)), max(vec))

if autorun:
    xn = graph.lemke.xn_anim[-1]
    xt = graph.lemke.xt_anim[-1]
    xnl = graph.lemke.xn_anim[0]
    xtl = graph.lemke.xt_anim[0]
    stress_t = []
    stress_n = []
    stress_tl = []
    stress_nl = []
    for i in range(len(xn)):
        stress_t.append(xt[i]/mesh_size/1e6) # MPa
        stress_n.append(xn[i]/mesh_size/1e6) # MPa
        stress_tl.append(xtl[i]/mesh_size/1e6) #MPa
        stress_nl.append(xnl[i]/mesh_size/1e6) #MPa
    print(f'MaxAbs: {maximum(stress_n)} Stress n: \n{stress_n}')
    print(f'MaxAbs: {maximum(stress_t)} Stress t: \n{stress_t}')
    print(f"MaxAbs: {maximum(stress_nl)} Stress nl: \n{stress_nl}")
    print(f'MaxAbs: {maximum(stress_tl)} Stress tl: \n{stress_tl}')
    print(f'Average Stress n: {sum(stress_n)/len(stress_n)}')
    print(f'Average Stress t: {sum(stress_t)/len(stress_t)}')


# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if autorun:
    mytable = PrettyTable()
    mytable.field_names = ['step', 'p', 'zn', 'xn', 'zt', 'xt']
    for i in range(len(graph.lemke.zn_anim)-1, len(graph.lemke.zn_anim)):
        mytable.add_row([i, graph.lemke.p_anim[i], graph.lemke.zn_anim[i], graph.lemke.xn_anim[i],
                         graph.lemke.zt_anim[i], graph.lemke.xt_anim[i]])
    print(mytable)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
