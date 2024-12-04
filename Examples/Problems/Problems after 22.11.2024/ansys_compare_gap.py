
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

assert FRICTION_COEFFICIENT == 0.5, 'Friction coef need to be 0.5 НЕПРАВИЛЬНО!'
assert PLANE_STRAIN is False, 'PLANE STRAIN need to be false! НЕПРАВИЛЬНО!!!!'
assert 1e-15 <= ACCURACY_OF_LCP <= 1e-10
print('Starting to calculate...')

start = time.time()

# set inputs
t = 0.1  # Thickness
E_bot = 2.5e10  # concrete
E_top = 2.5e9
mu = 0.2
# sizes in meters
L1 = 1
L2 = 2
L3_top = 1
L3_bot = 2
h_top = 2
h_bot = 1
gap1 = 0.001
gap2 = 0.002

# load
q = 100_000  # Н
F = 100_000  # N

print(f'Equivalent q: {q*L2}, Ultimate friction: {q*L2*FRICTION_COEFFICIENT}, {F=}')
print(f'Equivalent q is MORE or EQUAL than F: {q*L2 >= F}')

mesh_size = 0.1

force_inc = False
autorun = True
one_force_only = False  # TEST FORCE ON the RIGHT
one_F = 100_000
Excel = True

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
element_macro.add_element(EN=[0, 1, 2, 3], frag_amount_h=int(L1/ mesh_size), frag_amount_v=int(h_top / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False)  # 0
element_macro.add_element(EN=[1, 4, 5, 2], frag_amount_h=int(L2 / mesh_size), frag_amount_v=int(h_top / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[0])  # 1
element_macro.add_element(EN=[4, 6, 7, 5], frag_amount_h=int(L3_top / mesh_size), frag_amount_v=int(h_top / mesh_size),
                          E=E_top, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[1])  # 2

element_macro.add_element(EN=[8, 9, 1, 0], frag_amount_h=int(L1 / mesh_size), frag_amount_v=int(h_bot / mesh_size),
                          E=E_bot, mu=mu, t=t, own_weight=0, stitch=False)  # 3
element_macro.add_element(EN=[9, 10, 11, 1], frag_amount_h=int(L2 / mesh_size), frag_amount_v=int(h_bot / mesh_size),
                          E=E_bot, mu=mu, t=t, own_weight=0, stitch=False, stitch_list=[3])  # 4
element_macro.add_element(EN=[10, 12, 13, 11], frag_amount_h=int(L3_bot / mesh_size), frag_amount_v=int(h_bot / mesh_size),
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
# print(f'{contact_nodes0=}')
# print(f'{contact_nodes1_1=}')
# print(f'{contact_nodes1_2=}')
# print(f'{contact_nodes2_1=}')
# print(f'{contact_nodes2_2=}')
# NULL ELEMENTS
for i in range(0, len(contact_nodes0), 2):
    contact_pair = sorted(contact_nodes0[i:i + 2], reverse=True)
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2, gap_length=0)
    # print(f'{contact_pair[0], contact_pair[1]}')
for i in range(2, len(contact_nodes1_1)):
    contact_pair = [contact_nodes1_2[i], contact_nodes1_1[i]]
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2)  # +math.atan(gap1/L2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # print(f'{contact_pair[0], contact_pair[1]}')
for i in range(1, len(contact_nodes2_1)):
    contact_pair = [contact_nodes2_2[i], contact_nodes2_1[i]]
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=E_top, alpha=math.pi/2) # +math.atan(gap2/L3_top) !!!!!!!!!!!!!!!!!!!!!!!!!!
    # print(f'{contact_pair[0], contact_pair[1]}')

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)

sup_nodes_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L1+L2+L3_bot, 0))
sm.support_nodes(sup_nodes_bot, direction='hv')

# load cheme
lv = LoadVector()
lv_v = LoadVector()
force_nodes_q = nodes.find_nodes_numbers_along_segment(point1=(L1, h_bot+h_top),
                                                       point2=(L1+L2, h_bot+h_top))
force_node_F = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot+h_top),
                                                       point2=(0, h_bot+h_top))
force_node_F = force_node_F[0]

print("LOAD START======================================")
print(f'{q=}, {mesh_size=}')
nodes_under_load_amount = len(force_nodes_q)  # how many nodes under the load
force_in_one_node = - q * mesh_size
force_sum = 0
if not force_inc:
    if not one_force_only:
        for i in range(nodes_under_load_amount):
            degree_of_freedom = force_nodes_q[i]*2 + 1
            if i == 0 or i == nodes_under_load_amount - 1:
                lv.add_concentrated_force(force=force_in_one_node/2, degree_of_freedom=degree_of_freedom)
                force_sum += force_in_one_node/2
                print(f'force q: {force_in_one_node/2} in {force_nodes_q[i]}')
            else:
                lv.add_concentrated_force(force=force_in_one_node, degree_of_freedom=degree_of_freedom)
                force_sum += force_in_one_node
                print(f'force q: {force_in_one_node} in {force_nodes_q[i]}')

        print(f'{nodes_under_load_amount=}, {force_in_one_node=}')
        print(f'Sum of all forces: {force_sum}, should be {-q*L2}')
        assert force_sum == -q*L2, 'forces are WRONG!'

        print(f'Force {F=} in node {force_node_F=}')
        degree_of_freedom = force_node_F*2
        lv.add_concentrated_force(force=F, degree_of_freedom=degree_of_freedom)

    if one_force_only:
        node_one_force = nodes.find_nodes_numbers_along_segment(point1=(L1+L2+L3_top, h_bot + h_top),
                                                              point2=(L1+L2+L3_top, h_bot + h_top))
        node_one_force = node_one_force[0]
        print(f'ONLY ONE force {one_F/3=} in node {node_one_force=}')
        degree_of_freedom = node_one_force * 2 + 1
        lv.add_concentrated_force(force=-one_F/3, degree_of_freedom=degree_of_freedom)
        print(f'ONLY ONE force {one_F/3=} in node {node_one_force-1=}')
        degree_of_freedom = (node_one_force-1) * 2 + 1
        lv.add_concentrated_force(force=-one_F/3, degree_of_freedom=degree_of_freedom)
        print(f'ONLY ONE force {one_F/3=} in node {node_one_force-2=}')
        degree_of_freedom = (node_one_force - 2) * 2 + 1
        lv.add_concentrated_force(force=-one_F / 3, degree_of_freedom=degree_of_freedom)

else:
    pass
print("LOAD END======================================")

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=2, scale_def=1, autorun=autorun)

def maximum(vec: list):
    return max(abs(min(vec)), max(vec))

print("STRESS________________________________________________________________________________________________________")
# most right top node of bot plate
right_top_node_bot_plate = node_one_force = nodes.find_nodes_numbers_along_segment(point1=(L1+L2+L3_bot, h_bot),
                                                              point2=(L1+L2+L3_bot, h_bot))[0]
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
    print(f'MaxAbs Sn: {maximum(stress_n)}\n  Stress n: {stress_n}')
    print(f'MaxAbs St: {maximum(stress_t)}\n  Stress t: {stress_t}')
    print("linear |")
    print(f"MaxAbs Snl: {maximum(stress_nl)}\n  Stress nl: {stress_nl}")
    print(f'MaxAbs Stl: {maximum(stress_tl)}\n  Stress tl: {stress_tl}')
    print(f'maxabsXn: {maximum(xn)}\nmaxabsXt: {maximum(xt)}')
    print(f'Z top right nonlinear {graph.u_contact_anim[-1][force_node_F*2 + 1]}')
    print(f'Z vertical linear most right contact node top plate {graph.u_linear_const[contact_nodes2_1[-1]*2 + 1]}\n'
          f'the node number is: {contact_nodes2_1[-1]}')
    print(f'Z vertical NONlinear most right top node of bot plate {graph.u_contact_anim[-1][right_top_node_bot_plate * 2 + 1]}\n'
          f'the node number is: {right_top_node_bot_plate}')
    print(f'SUM xn: {sum(xn)}, SUM xt: {sum(xt)}')

    if Excel:
        workbook = xlsxwriter.Workbook('LCP results.xlsx')
        worksheet = workbook.add_worksheet()
        column = 0
        names = ['xn', 'xt', 'xnl', 'xtl', 'stress_n', 'stress_t', 'stress_nl', 'stress_tl']
        for vec in [xn, xt, xnl, xtl, stress_n, stress_t, stress_nl, stress_tl]:
            worksheet.write(0, column, names[column])
            row = 1
            for item in vec:
                # write operation perform
                worksheet.write(row, column, item)
                row += 1
            column += 1

        workbook.close()

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

    # for i, (p, xn, xt) in enumerate(zip(graph.lemke.p_anim, graph.lemke.xn_anim, graph.lemke.xt_anim)):
    #     print(f'{i}: sum Xn: {sum(xn)-p*len(xn)}, sum Xt: {sum(xt)}')

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
