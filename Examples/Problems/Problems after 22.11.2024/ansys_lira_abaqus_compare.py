
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

assert FRICTION_COEFFICIENT == 0.25, 'Friction coef need to be 0.5 НЕПРАВИЛЬНО!'
assert PLANE_STRAIN is False, 'PLANE STRAIN need to be false! НЕПРАВИЛЬНО!!!!'
assert 1e-15 <= ACCURACY_OF_LCP <= 1e-10
print('Starting to calculate...')

start = time.time()

# set inputs
t = 0.1  # Thickness
E_bot = 2.5e10  # concrete
E_top = 2.5e10
mu = 0.2
# sizes
L_bot = 2
L_top = L_bot
h = 2
h_bot_left = 1.5
h_bot_right = 1
h_top_left = h - h_bot_left
h_top_right = h - h_bot_right
# load
q = 100_000 # Н
# q = q * 40

Lq = L_top/10 # m    L_top / 10

force_inc = False
autorun = True

mesh_size = 0.05
left_top_sup = False
force_in_one_node_flag = False

# add nodes  for frame element
nodes = NodeContainer()
# nodes for bot macro element
nodes.add_node(0, 0) # 0
nodes.add_node(L_bot, 0) # 1
nodes.add_node(L_bot, h_bot_right) # 2
nodes.add_node(0, h_bot_left) # 3
# nodes for top macro element
nodes.add_node(0, h_bot_left) # 4
nodes.add_node(L_top, h_bot_right) # 5
nodes.add_node(L_top, h) # 6
nodes.add_node(0, h) # 7

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
contact_nodes = nodes.find_nodes_numbers_along_segment(point1=(0, h_bot_left), point2=(L_top, h_bot_right),
                                                       sorted_by_y=False, relative_tolerance=1e-2)

for i in range(0, len(contact_nodes), 2):
    contact_pair = sorted(contact_nodes[i:i+2])
    element_null.add_element(EN=[contact_pair[0], contact_pair[1]], cke=111, alpha=(math.pi/2-math.atan(0.5/L_top)))
    # print(f'{contact_pair[0], contact_pair[1]}')

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)

sup_nodes_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(L_bot, 0))
sm.support_nodes(sup_nodes_bot, direction='hv')

if left_top_sup:
    sup_nodes_top_right = nodes.find_nodes_numbers_along_segment(point1=(L_top, h_bot_right + mesh_size/2), point2=(L_top, h))
    print(f'RIGHT SUPPORT IS ENABLED')
    sm.support_nodes(sup_nodes_top_right, direction='hv')


# load cheme
lv = LoadVector()
lv_v = LoadVector()
force_node = int(L_bot/mesh_size + 1) * int((h_bot_left+h_bot_right)/2/mesh_size + 1) + int(L_top/mesh_size + 1) * int((h_top_left+h_top_right)/2/mesh_size + 1)

print("LOAD ======================================")
print(f'{q=}, {Lq=}, {mesh_size=}, {force_node=}')
nodes_under_load = int(Lq / mesh_size) + 1  # how many nodes under the load
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
      assert force_sum==-q*Lq, 'forces are WRONG!'

    if force_in_one_node_flag:
       print(f'One force')
       degree_of_freedom = force_node*2 + 1
       lv.add_concentrated_force(force=-q , degree_of_freedom=degree_of_freedom)

else:
    pass
print("LOAD ======================================")

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=2, scale_def=200, autorun=autorun)

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
    print(f'Z 1st contact LINEAR: {graph.u_contact_anim[0][contact_nodes[0]*2 + 1]} at node: {contact_nodes[0]}')
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
