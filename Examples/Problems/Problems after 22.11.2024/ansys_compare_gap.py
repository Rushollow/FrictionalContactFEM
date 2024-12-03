
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

assert FRICTION_COEFFICIENT == 0.25, 'Friction coef need to be 0.25 НЕПРАВИЛЬНО!'
assert PLANE_STRAIN is False, 'PLANE STRAIN need to be false! НЕПРАВИЛЬНО!!!!'
assert 1e-15 <= ACCURACY_OF_LCP <= 1e-10
print('Starting to calculate...')

start = time.time()

# set inputs
t = 1  # Thickness
E_bot = 2.5e10  # concrete
E_top = 2.5e9
mu = 0.2
# sizes in meters
L1 = 1
L2 = 1
L3_top = 1
L3_bot = 1
h_top = 1
h_bot = 1
gap1 = 0.2
gap2 = gap1*2

# load
Lq = L3_top
q = 100_000  # Н
F = 1_000  # N

mesh_size = 0.5

force_inc = False
autorun = True
force_in_one_node_flag = True

# add nodes  for frame element
nodes = NodeContainer()
# nodes for bot macro element 0
nodes.add_node(0, h_bot)  # 0
nodes.add_node(L1, h_bot)  # 1
nodes.add_node(L1, h_bot+h_top)  # 2
nodes.add_node(0, h_bot+h_top)  # 3
# nodes for top macro element 1
nodes.add_node(L1+L2, h_bot+gap1)  # 4
nodes.add_node(L1+L2, h_bot+h_top)  # 5
# nodes for top macro element 2
nodes.add_node(L1+L2+L3_top, h_bot+gap2)  # 6
nodes.add_node(L1+L2+L3_top, h_bot+h_top)  # 7
# nodes for top macro element 3
nodes.add_node(0, 0)  # 8
nodes.add_node(L1, 0)  # 9
# nodes for top macro element 4
nodes.add_node(L1+L2, 0)  # 10
nodes.add_node(L1+L2, h_bot)  # 11
# nodes for top macro element 5
nodes.add_node(L1+L2+L3_bot, 0)  # 12
nodes.add_node(L1+L2+L3_bot, h_bot)  # 13

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# make grid
element_macro.add_element(EN=[0, 1, 2, 3], frag_amount_h=int(h_top/ mesh_size), frag_amount_v=int(L1 / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False)  # 0
element_macro.add_element(EN=[1, 4, 5, 2], frag_amount_h=int(h_top / mesh_size), frag_amount_v=int(L2 / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[0])  # 1
element_macro.add_element(EN=[4, 6, 7, 5], frag_amount_h=int(h_top / mesh_size), frag_amount_v=int(L3_top / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[1])  # 2

element_macro.add_element(EN=[8, 9, 1, 0], frag_amount_h=int(h_bot / mesh_size), frag_amount_v=int(L1 / mesh_size),
                          E=E_bot, mu=mu, t=t, own_weight=0, stitch=False)  # 3
element_macro.add_element(EN=[9, 10, 11, 1], frag_amount_h=int(h_bot / mesh_size), frag_amount_v=int(L2 / mesh_size),
                          E=E_bot, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[3])  # 4
element_macro.add_element(EN=[10, 12, 13, 11], frag_amount_h=int(h_bot / mesh_size), frag_amount_v=int(L3_bot / mesh_size),
                          E=E_bot, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[4])  # 5
print()

element_macro.fragment_all(element_4node=element_4node, element_frame=element_frame, element_null=element_null)

# add null elements
contact_nodes0 = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot), point2=(L1, h_bot),
                                                        sorted_by_y=False, relative_tolerance=1e-10)
# 2nd zone
contact_nodes1_1 = nodes.find_nodes_numbers_along_segment(point1=(L1, h_bot), point2=(L1 + L2, h_bot + gap1),
                                                          sorted_by_y=False, relative_tolerance=1e-10)  # top 1-4
contact_nodes1_2 = nodes.find_nodes_numbers_along_segment(point1=(L1, h_bot), point2=(L1 + L2, h_bot),
                                                          sorted_by_y=False, relative_tolerance=1e-10)  # bot 1-11
# 3d zone
contact_nodes2_1 = nodes.find_nodes_numbers_along_segment(point1=(L1 + L2, h_bot + gap1), point2=(L1 + L2 + L3_top, h_bot + gap2),
                                                          sorted_by_y=False, relative_tolerance=1e-10)  # top 4-6
contact_nodes2_2 = nodes.find_nodes_numbers_along_segment(point1=(L1 + L2, h_bot), point2=(L1 + L2 + L3_top, h_bot),
                                                          sorted_by_y=False, relative_tolerance=1e-10)  # bot 11-13(6)
print(f'{contact_nodes0=}')
print(f'{contact_nodes1_1=}')
print(f'{contact_nodes1_2=}')
print(f'{contact_nodes2_1=}')
print(f'{contact_nodes2_2=}')
# NULL ELEMENTS
for i in range(0, len(contact_nodes0), 2):
    contact_pair = sorted(contact_nodes0[i:i + 2], reverse=True)
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2, gap_length=0)
    print(f'{contact_pair[0], contact_pair[1]}')
for i in range(2, len(contact_nodes1_1)):
    contact_pair = [contact_nodes1_2[i], contact_nodes1_1[i]]
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2)  # math.atan(gap1/L2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print(f'{contact_pair[0], contact_pair[1]}')
for i in range(1, len(contact_nodes2_1)):
    contact_pair = [contact_nodes2_2[i], contact_nodes2_1[i]]
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2) # math.atan(gap2/L2) !!!!!!!!!!!!!!!!!!!!!!!!!!
    print(f'{contact_pair[0], contact_pair[1]}')

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)

sup_nodes_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L1+L2+L3_bot, 0))
sm.support_nodes(sup_nodes_bot, direction='hv')

# load cheme
lv = LoadVector()
lv_v = LoadVector()
force_nodes_q = nodes.find_nodes_numbers_along_segment(point1=(L1+L2+L3_top-Lq, h_bot+h_top),
                                                       point2=(L1+L2+L3_top, h_bot+h_top))
force_node = force_nodes_q[-1]

print("LOAD ======================================")
print(f'{q=}, {Lq=}, {mesh_size=}, {force_node=}')
nodes_under_load = len(force_nodes_q)  # how many nodes under the load
force_in_one_node = - q * mesh_size
force_sum = 0
if not force_inc:
    if not force_in_one_node_flag:
        for i in range(nodes_under_load):
            degree_of_freedom = (force_node + i)*2 + 1
            if i == 0 or i == nodes_under_load - 1:
                lv.add_concentrated_force(force=force_in_one_node/2, degree_of_freedom=degree_of_freedom)
                force_sum += force_in_one_node/2
                print(f'force: {force_in_one_node/2} in {force_node + i}')
            else:
                lv.add_concentrated_force(force=force_in_one_node, degree_of_freedom=degree_of_freedom)
                force_sum += force_in_one_node
                print(f'force: {force_in_one_node} in {force_node + i}')

        print(f'{nodes_under_load=}, {force_in_one_node=}')
        print(f'Sum of all forces: {force_sum}, should be {-q*Lq}')
        assert force_sum == -q*Lq, 'forces are WRONG!'

    if force_in_one_node_flag:
        print(f'One force')
        degree_of_freedom = force_node*2 + 1
        lv.add_concentrated_force(force=-q, degree_of_freedom=degree_of_freedom)

else:
    pass
print("LOAD ======================================")

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=2, scale_def=1, autorun=autorun)

def maximum(vec: list):
    return max(abs(min(vec)), max(vec))

print("STRESS________________________________________________________________________________________________________")
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
        k = 1
        if i == 0 or i == len(xn)-1: # if elements on the edge
            k = 2
        stress_n.append(xn[i]/mesh_size/t*k) # Pa
        stress_t.append(xt[i]/mesh_size/t*k) # Pa
        stress_nl.append(xnl[i]/mesh_size/t*k) # Pa
        stress_tl.append(xtl[i] / mesh_size / t * k)  # Pa
    print(f'MaxAbs: {maximum(stress_n)}\nStress n: \n{stress_n}')
    print(f'MaxAbs: {maximum(stress_t)}\n Stress t: \n{stress_t}')
    print(f"MaxAbs: {maximum(stress_nl)}\n Stress nl: \n{stress_nl}")
    print(f'MaxAbs: {maximum(stress_tl)}\n Stress tl: \n{stress_tl}')
    print(f'maxabsXn: {maximum((xn))}\nmaxabsXt: {maximum(xt)}')
    print(f'Z top rigt nonlinear {graph.u_contact_anim[-1][force_node*2 + 1]}')
    print(f'Z 1st contact LINEAR: {graph.u_contact_anim[0][contact_nodes0[0] * 2 + 1]} at node: {contact_nodes0[0]}')
    print(f'SUM xn: {sum(xn)}, SUM xt: {sum(xt)}')

print("STRESS________________________________________________________________________________________________________")


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
