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

assert FRICTION_COEFFICIENT == 0.4, 'Friction coef need to be 0.4 НЕПРАВИЛЬНО!'
assert PLANE_STRAIN is True, 'PLANE STRAIN need to be true! НЕПРАВИЛЬНО!!!!'
assert ACCURACY_OF_LCP >= 1e-15
print('Starting to calculate...')

start = time.time()

# set inputs
tc = 1  # Thickness
Ec = 26500e6  # бетон stiffness retaining wall
gamma_c = 24e3  # Н/m^3  own weight for concrete wall
mu_c = 0.15
# Двутавр стальной горячекатный
Er = 2e11  # Pa for beam
Ar = 2.68e-3  # m^2 for beam
Ixr = 1.84e-5  # m^4 for beam
# load
q = 2.5e3 # N/m
# gamma_r = 1  # for beam
force_inc = False
autorun = True

Lp = 3  # half of prol`ot length (span)
Lw = 0.3 # thikness of wall
Lg = 0.1  # how much beam is now reaching towards outside wall
Lb = Lp - Lg  # how long is half of the beam
hw = 3 * Lw  # height of the wall in calculations
if Lg > 0.05:
    mesh_size = Lg/4
else:
    mesh_size = 0.05

# add nodes  for frame element
nodes = NodeContainer()
for i in range(int(Lb/mesh_size)):
    nodes.add_node(Lg+i*mesh_size, hw)
# add nodes for ME
nodes.add_node(0, 0)
nodes.add_node(Lw, 0)
nodes.add_node(Lw, hw)
nodes.add_node(0, hw)
# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)

# add frame elements
for i in range(int(Lb/mesh_size)-1):
    element_frame.add_element(EN=[i, i+1], E=Er, A=Ar, I=Ixr)
# add MEs
nn = int(Lb/mesh_size)
element_macro.add_element(EN=[nn, nn+1, nn+2, nn+3], frag_size=mesh_size, E=Ec, mu=mu_c, t=tc,
                          own_weight=gamma_c, stitch=False)

element_macro.fragment_all(element_4node, element_frame, element_null)
# find nodes for contact pairs
n_contact = nodes.find_nodes_numbers_along_segment((Lg, hw), (Lw, hw))

# n null elements and adding t null elements silently
for i in range(0, len(n_contact), 2):
    if n_contact[i+1] > n_contact[i]:
        bot, top = n_contact[i+1], n_contact[i]
    else:
        top, bot = n_contact[i + 1], n_contact[i]
    element_null.add_element(EN=[bot, top], cke=123, alpha=math.pi / 2, gap_length=0)

    print(f'bot node: {n_contact[i+1]}')
    print(f'top node: {n_contact[i]}')

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_sup_bot = nodes.find_nodes_numbers_along_segment(point1=(0, 0), point2=(Lw, 0))
node_sup_frame = [nn-1]

sm.support_nodes(nodes_to_sup_bot, direction='hv')
sm.support_nodes(node_sup_frame, direction='hr')

lv = LoadVector()
lv_v = LoadVector()

lv.add_own_weight_to_rf(nodes, [element_4node])
# add load at the top of the beam
for i in range(nn):
    lv.add_concentrated_force(-q*mesh_size, degree_of_freedom =i * 2 + 1)

if not force_inc:
    pass
    # lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])
else:
    pass
    # lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])

# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=lv_v,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=2, scale_def=10, autorun=autorun)


if autorun:
    xn = graph.lemke.xn_anim[-1]
    xt = graph.lemke.xt_anim[-1]
    stress_x = []
    stress_y = []
    for i in range(len(xn)):
        stress_x.append(-xt[i]/mesh_size/1000000)  # MPa
        stress_y.append(-xn[i]/mesh_size/1000000)   # MPa
    print(stress_x)
    print(stress_y)
    x = np.arange(1, len(stress_x) + 1, 1, dtype=int)
    stress_x = np.array(stress_x)
    stress_y = np.array(stress_y)
    # pg.plot(x, stress_x, pen=None, symbol='o')  # setting pen=None disables line drawing
    # pg.plot(x, stress_y, pen=None, symbol='o')  # setting pen=None disables line drawing

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
    # print(mytable)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
