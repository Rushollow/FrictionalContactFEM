import time
import numpy as np
from matplotlib import pyplot as plt
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector, solve_slae
from FEM.element_null import ElementNullContainer  # to add null-element
from FEM.element_4node import Element4NodeLinearContainer  # to add 4node element
from FEM.element_frame import ElementFrameContainer  # to add frame element
from SchemeForm.macro_element import ElementMacroContainer
from LCP.initial_table import InitialTable  # to form initial table for LCP
from LCP.lemke import Lemke  # to solve LCP
from input_data import FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP

from Visualize.plot_data_scheme import PlotScheme  # for visualizing
from GUI.tkinter_gui import ContactFEM
from constants_data import Arch
from input_data import ARCH_SPAN, ARCH_HEIGHT

# set variables scheme
num_me_in_quarter = 10  # 15 on how many parts quarter with ME will be divided. Better be more than frag_amount_v
frag_amount_v = 10  # 14 how much each quarter section will be fragmented in vertical direction
frag_amount_h = 1  # how much each quarter section will be fragmented in horizontal direction (USE num_me_in_quarter instead!!)
frag_amount_v_under_arch = frag_amount_v  # better be equal to frag_amount_v
right_border = ARCH_SPAN + 20
top_border = ARCH_HEIGHT + 20
bottom_border = -10
left_border = -right_border
# elements variables
E_soil = 4.5e7  # Pa
mu_soil = 0.27  #
t_soil = 1  # m
own_weight_soil = 22661.1  # N/m^3
E_pile = 2e10  # Pa
A_pile = 0.25  # m^2
I_pile = 0.005208333333333333  # m^4
E_arch = 2e11
A_arch = 0.02366322
I_arch = 2.110901e-4

start = time.time()
arch = Arch()
def outline(x):
    return arch.outline_circle(x)
def linear_fnc(x, alpha):
    return np.tan(alpha) * x
def top_border_fnc(x):
    return top_border
def bisection(outline, line, alpha, tol, a, b):
    i = 0
    xc = None
    xl = a
    xr = b
    while np.abs(xl - xr) > tol:
        xc = (xl + xr) / 2.0
        prod = (outline(xl) - line(xl, alpha)) * (outline(xc) - line(xc, alpha))
        if prod > 0:
            xl = xc
        else:
            xr = xc
        i += 1
        if i > 100000:
            xc = None
    if xc == None:
        raise ValueError('HEY WHAT JUST HAPPENED?!')
    return xc

nodes = NodeContainer()
nodes_frame = NodeContainer()
# form nodes for ME SOIL
# around arch for soil on right
nodes.add_node(0, ARCH_HEIGHT)
nodes.add_node(0, top_border)
for i in range(num_me_in_quarter - 1, 0, -1):
    alpha_step = np.pi / 2 / num_me_in_quarter
    alpha = alpha_step * i
    x1 = bisection(outline, linear_fnc, alpha, 1e-3, a=0, b=ARCH_SPAN)
    y1 = outline(x1)
    nodes.add_node(x1, y1)  # add soil node
    y2 = linear_fnc(right_border, alpha)
    x2 = right_border
    if y2 > top_border:
        x2 = bisection(top_border_fnc, linear_fnc, alpha, tol=1e-3, a=0, b=right_border)
        y2 = top_border
    nodes.add_node(x2, y2)
nodes.add_node(ARCH_SPAN/2, 0)
nodes.add_node(right_border, 0)
# add nodes underneath right side of soil
nodes.add_node(ARCH_SPAN/2, 0)
nodes.add_node(ARCH_SPAN/2, bottom_border)
nodes.add_node(right_border, bottom_border)
nodes.add_node(right_border, 0)
# copy all nodes from right side to left
for i in range(len(nodes)):
    node = nodes[i]
    nodes.add_node(-node.x, node.y)
# make nodes for central bottom ME
nodes.add_node(-ARCH_SPAN / 2, bottom_border)
nodes.add_node(ARCH_SPAN / 2, bottom_border)
nodes.add_node(ARCH_SPAN / 2, 0)
nodes.add_node(-ARCH_SPAN / 2, 0)

# FRAME NODES
# add right pile nodes
for i in range(frag_amount_v_under_arch-1, 0, -1):
    vertical_step = bottom_border/frag_amount_v_under_arch
    nodes_frame.add_node(ARCH_SPAN / 2, i * vertical_step)
alpha_arch_null_el = []  # angle for the arch outline and horizontal axis
# add arch node on the right bottom of the arch
nodes_frame.add_node(ARCH_SPAN / 2, 0)
for i in range(1, num_me_in_quarter, 1):
    alpha_step = np.pi / 2 / num_me_in_quarter
    alpha = alpha_step * i
    x1 = bisection(outline, linear_fnc, alpha, tol=1e-3, a=0, b=ARCH_SPAN)
    x_a = bisection(outline, linear_fnc, alpha - 0.01, tol=1e-3, a=0, b=ARCH_SPAN)
    y_a = outline(x_a)
    x_b = bisection(outline, linear_fnc, alpha + 0.01, tol=1e-3, a=0, b=ARCH_SPAN)
    y_b = outline(x_b)
    angle = np.arctan((y_a - y_b) / (x_b - x_a))  # find alpha angle from derivative
    alpha_arch_null_el.append(angle)
    y1 = outline(x1)
    nodes_frame.add_node(x1, y1)  # add arch node
# add mode alpha(s) to list
alpha_arch_null_el_left = sorted([np.pi - i for i in alpha_arch_null_el])
alpha_arch_null_el.sort()
alpha_arch_null_el.append(np.pi/2)
alpha_arch_null_el += alpha_arch_null_el_left
# add central node for arch
nodes_frame.add_node(0, ARCH_HEIGHT)
# copy nodes of the arch and pile to the left
for i in range(len(nodes_frame) - 2, -1, -1):
    node = nodes_frame[i]
    nodes_frame.add_node(-node.x, node.y)
amount_of_nodes_before_fragment = len(nodes)
# add frame nodes to scheme nodes
nodes.nodes_list += nodes_frame.nodes_list
nodes_amount_frame = len(nodes_frame)  # how many nodes used for frame elements

# Find node numbers for null elements
nodes_right_pile = []
nodes_right_pile_right_soil = []
nodes_right_pile_left_soil = []
nodes_left_pile = []
nodes_left_pile_left_soil = []
nodes_left_pile_right_soil = []
nodes_soil_above_arch = []
nodes_arch = []
for i in range(len(nodes)):
    if i < frag_amount_v_under_arch:
        nodes_right_pile.append(i)
    elif i >= frag_amount_v_under_arch + num_me_in_quarter * 2 - 1 and\
            i < frag_amount_v_under_arch * 2 + num_me_in_quarter * 2 - 1:
        nodes_left_pile.append(i)
    elif i < nodes_amount_frame:
        nodes_arch.append(i)
nodes_right_pile_right_soil.append(len(nodes) - amount_of_nodes_before_fragment +
                                   (frag_amount_v + 1) * (num_me_in_quarter + 1) - frag_amount_v - 1)
for i in range(len(nodes) - amount_of_nodes_before_fragment + (frag_amount_v + 1) * (num_me_in_quarter + 1),
               len(nodes) - amount_of_nodes_before_fragment + (frag_amount_v + 1) * (num_me_in_quarter + 1) + frag_amount_v_under_arch - 1):
    nodes_right_pile_right_soil.append(i)
nodes_right_pile_right_soil.sort(reverse=True)
for i in range(nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) * 2 +
               (frag_amount_v + 1) * num_me_in_quarter + (frag_amount_v_under_arch + 1)*2,
               nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) * 2 +
               (frag_amount_v + 1) * num_me_in_quarter + (frag_amount_v_under_arch + 1) * 2 +
               (frag_amount_v_under_arch + 1) * (frag_amount_v_under_arch - 1) + 1,
               frag_amount_v_under_arch + 1):
    nodes_right_pile_left_soil.append(i)
nodes_left_pile_left_soil.append(nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) +
                                 frag_amount_v_under_arch * (frag_amount_v + 1) + (frag_amount_v + 1) * num_me_in_quarter - frag_amount_v)
for i in range(nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) +
               (frag_amount_v + 1) * num_me_in_quarter + 1,
               nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) +
               (frag_amount_v + 1) * num_me_in_quarter+ frag_amount_v_under_arch):
    nodes_left_pile_left_soil.append(i)
for i in range(nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) * 2 +
               (frag_amount_v + 1) * num_me_in_quarter + frag_amount_v_under_arch + 2,
               nodes_amount_frame - 1 + (num_me_in_quarter + 1) * (frag_amount_v + 1) + frag_amount_v_under_arch * (frag_amount_v + 1) * 2 +
               (frag_amount_v + 1) * num_me_in_quarter + frag_amount_v_under_arch + 2 + frag_amount_v_under_arch * (frag_amount_v_under_arch + 1),
               frag_amount_v_under_arch + 1):
    nodes_left_pile_right_soil.append(i)
nodes_left_pile_right_soil.sort(reverse=True)
# nodes for soil above arch
nodes_soil_above_arch.append(nodes_amount_frame)
nodes_soil_above_arch.append(nodes_amount_frame + 1)
for i in range(nodes_amount_frame + 2 + frag_amount_v * 2,
               nodes_amount_frame + (frag_amount_v + 1) * (num_me_in_quarter - 1) + 1,
               frag_amount_v + 1):
    nodes_soil_above_arch.append(i)
nodes_soil_above_arch.sort(reverse=True)
for i in range(nodes_amount_frame + (frag_amount_v + 1)*(num_me_in_quarter + 1) + (frag_amount_v + 1)*(frag_amount_v_under_arch),
               nodes_amount_frame + (frag_amount_v + 1)*(num_me_in_quarter + 1) + (frag_amount_v + 1)*(frag_amount_v_under_arch) +
               (num_me_in_quarter - 2) * (frag_amount_v + 1) + 1,
               frag_amount_v + 1):
    nodes_soil_above_arch.append(i)

# Set elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
# add arch and pile frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
# add frame elements
for i in range(amount_of_nodes_before_fragment, len(nodes) - 1):
    k = i - amount_of_nodes_before_fragment  # number of node if frame nodes would be first
    if k in nodes_right_pile[:-1] or k in nodes_left_pile[:-1]:
        element_frame.add_element(EN=[i, i+1], E=E_pile, A=A_pile, I=I_pile)
    else:
        element_frame.add_element(EN=[i, i+1], E=E_arch, A=A_arch, I=I_arch)

element_null = ElementNullContainer(nodes_scheme=nodes)
element_macro = ElementMacroContainer(nodes_scheme=nodes)
# add ME for right side arch soil
for i in range(num_me_in_quarter):
    n1 = 2 * i
    n2 = n1 + 2
    n3 = n1 + 3
    n4 = n1 + 1
    element_macro.add_element(EN=[n1, n2, n3, n4], E=E_soil, mu=mu_soil, t=t_soil, own_weight=own_weight_soil,
                              frag_amount_h=frag_amount_h, frag_amount_v=frag_amount_v, stitch=True)
# add ME right side bottom soil
n1 = num_me_in_quarter * 2 + 2
element_macro.add_element(EN=[n1, n1+1, n1+2, n1+3], E=E_soil, mu=mu_soil, t=t_soil, own_weight=own_weight_soil,
                          frag_amount_h=frag_amount_v_under_arch, frag_amount_v=frag_amount_v, stitch=True)
# # add ME for right side arch soil
for i in range(num_me_in_quarter):
    nodes_amount_right = num_me_in_quarter * 2 + 2 + 4  # *2 +2 - amount of nodes for quarter, + 4 - nodes for bottom ME
    n1 = 2 * i + nodes_amount_right
    n2 = n1 + 2
    n3 = n1 + 3
    n4 = n1 + 1
    element_macro.add_element(EN=[n1, n2, n3, n4], E=E_soil, mu=mu_soil, t=t_soil, own_weight=own_weight_soil,
                              frag_amount_h=frag_amount_h, frag_amount_v=frag_amount_v, stitch=True)
# add ME left side bottom soil
n1 = (num_me_in_quarter * 2 + 2) * 2 + 4
element_macro.add_element(EN=[n1, n1+1, n1+2, n1+3], E=E_soil, mu=mu_soil, t=t_soil, own_weight=own_weight_soil,
                          frag_amount_h=frag_amount_v_under_arch, frag_amount_v=frag_amount_v, stitch=True)
# add ME under the arch in the center
n1 = (num_me_in_quarter * 2 + 2) * 2 + 4 * 2
element_macro.add_element(EN=[n1, n1+1, n1+2, n1+3], E=E_soil, mu=mu_soil, t=t_soil, own_weight=own_weight_soil,
                          frag_amount_h=frag_amount_v_under_arch, frag_amount_v=frag_amount_v_under_arch, stitch=False)

element_macro.fragment_all(element_4node=element_4node, element_frame=element_frame, element_null=element_null)

# add NULL-ELEMENTS after fragmentation
for i in range(len(nodes_right_pile)):
    # add null_el to the right pile
    element_null.add_element(EN=[nodes_right_pile[i], nodes_right_pile_right_soil[i]], cke=E_arch*A_arch, alpha=0)
# add nulls to the arch and soil
for i in range(len(nodes_arch)):
    element_null.add_element(EN=[nodes_arch[i], nodes_soil_above_arch[i]], cke=E_arch*A_arch, alpha=alpha_arch_null_el[i])
for i in range(len(nodes_right_pile)):
    # add null_el to the left pile
    element_null.add_element(EN=[nodes_left_pile[i], nodes_left_pile_left_soil[i]], cke=E_arch*A_arch, alpha=np.pi)
# add nulls between central soil under arch and piles
for i in range(len(nodes_right_pile)):
    # add null_el to the right pile
    element_null.add_element(EN=[nodes_right_pile_left_soil[i], nodes_right_pile[i]], cke=E_arch*A_arch, alpha=0)
for i in range(len(nodes_right_pile)):
    # left pile
    element_null.add_element(EN=[nodes_left_pile_right_soil[i], nodes_left_pile[i]], cke=E_arch*A_arch, alpha=np.pi)

nodes_to_support_bottom = []
nodes_to_support_side = []
external_force_nodes = []
for i, node in enumerate(nodes):
    if np.abs(node.y - bottom_border) < 0.01:
        nodes_to_support_bottom.append(i)
    if np.abs(node.x - left_border) < 0.01 or np.abs(node.x - right_border) < 0.01:
        nodes_to_support_side.append(i)
    if np.abs(node.y - top_border) < 0.01:
        external_force_nodes.append(i)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
sm.support_nodes(nodes_to_support_bottom, direction='v')
sm.support_nodes(nodes_to_support_side, direction='h')
# sm.support_nodes([0, len(nodes_frame)-1], direction='hv')  # support both piles
lv = LoadVector()
F = 0
i = 1
# for node_num in external_force_nodes:
#     if i == 10 or i == 11 or i == 12:
#         lv.add_concentrated_force(-F, node_num * 2 + 1)
#         print(f'load:{-F} to {node_num}')
    # lv.add_concentrated_force(-F, node_num * 2 + 1)
    # print(f'load:{-F} to {node_num}')
    # i += 1
lv.add_own_weight_to_rf(nodes_scheme=nodes, element_container_list=[element_4node])

u_linear = solve_slae(sm, lv)

# set to show only first 5 numbers when printing numpy values
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# do stuff about contact SM and LV
intl_table = InitialTable(element_null, sm, lv, u_linear)
intl_table.form_initial_table()
# do lemke ot solve LCP
lemke = Lemke(intl_table)
lemke.lcp_solve()

# # plot info
# print(f'nodes_right_pile:           {nodes_right_pile}\n'
#       f'nodes_right_pile_right_soil:{nodes_right_pile_right_soil}\n'
#       f'nodes_right_pile_left_soil: {nodes_right_pile_left_soil}\n'
#       f'nodes_left_pile:               {nodes_left_pile}\n'
#       f'nodes_left_pile_left_soil:     {nodes_left_pile_left_soil}\n'
#       f'nodes_left_pile_right_soil:    {nodes_left_pile_right_soil}\n'
#       f'nodes_arch:             {nodes_arch}\n'
#       f'nodes_soil_above_arch:  {nodes_soil_above_arch}\n'
#       f'alpha_arch_null_el      {alpha_arch_null_el}')
# for i, el in enumerate(element_frame):
#     print(f'frame {i} EN:{el.EN}, MI:{el.MI}')
# for i, el in enumerate(element_4node):
#     print(f'4node {i} EN:{el.EN}, MI:{el.MI}')
# for i, el in enumerate(element_null):
#     print(f'null {i} EN:{el.EN}, MI:{el.MI}, alpha:{el.alpha}')


# plot --------------------------------------------------------------------------
# calculate data to plot
graph = PlotScheme(nodes, element_null, sm, lv, u_linear, lemke,
                   element_container_obj=element_4node, element_frame=element_frame,
                   partition=10, scale_def=10, text=False)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)
print(f'Right Border:{right_border}m\n'
      f'Maximum lemke zn:{max(lemke.xn_anim[-1])},\n'
      f'Minimum u_linear{min(u_linear)},\n'
      f'Minimum u_contact:{min(graph.u_contact_anim[-1])}')
# min_u_frame = 0
# for el in element_frame:
#     for node_num in el.EN:
#         if graph.u_contact_anim[-1][nodes[node_num].indices[1]] < min_u_frame:
#             min_u_frame = graph.u_contact_anim[-1][nodes[node_num].indices[1]]
# print(f'Minimum u for frame nodes:{min_u_frame}')
# for i, node in enumerate(nodes):
#     print(f'{i} '
#           f'linear h:{str("%.4f" % u_linear[node.indices[0]])} v:{str("%.4f" % u_linear[node.indices[1]])}'
#           f'contact h:{str("%.4f" % graph.u_contact_anim[-1][node.indices[0]])}, v:{str("%.4f" % graph.u_contact_anim[-1][node.indices[1]])},'
#           f' load:{lv.rf[node.indices[1]]}')
# smt_sub_sm = sm.r.T - sm.r
# print(f'SM rank:{np.linalg.matrix_rank(sm.r)}, rows:{sm.r.shape[1]}, cond:{np.linalg.cond(sm.r)} smT-sm max:{np.max(smt_sub_sm)}, min:{np.min(smt_sub_sm)}')

# Contact info
print(f'zn: {lemke.zn}\nzt:{lemke.zt}')
print(f'xn: {lemke.xn}\nxt:{lemke.xt}')
print(f'min zn:{min(lemke.zn)}, max zn: {max(lemke.zn)}, min zt:{min(lemke.zt)}, max zt:{max(lemke.zt)}')
print(f'min xn:{min(lemke.xn)}, max xn: {max(lemke.xn)}, min xt:{min(lemke.xt)}, max xt:{max(lemke.xt)}')
print(f'sum xn:{sum(lemke.xn)}, sum xt{sum(lemke.xt)}')
print(f'max_u_linear: {max(u_linear)}, min_u_linear: {min(u_linear)}')
print(f'FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP:{FRICTION_COEFFICIENT, PLANE_STRAIN, ACCURACY_OF_LCP}')

app = ContactFEM(graph=graph)
app.mainloop()

ps = PlotScheme(nodes, element_null, sm, lv, u_linear=u_linear, lemke=lemke, element_container_obj=element_4node,
                element_frame=element_frame, partition=10, scale_def=200)

# Scheme with nodes numbers
fig1, ax1 = plt.subplots(1, 1)
plt.gca().set_aspect('equal', adjustable='box')
for i, node in enumerate(nodes.nodes_list):
    ax1.scatter(node.x, node.y, color='lightgreen', zorder=3, marker='s')
ps.plot_list_of_elements(ax1)
for i, node in enumerate(nodes.nodes_list):
    ax1.scatter(node.x, node.y, color='lightgreen', zorder=3)
    dx = 0.2
    dy = 0.2
    if node.x < 0:
        dx *= -1
    if node.y < 0:
        dy *= -1
    if i < nodes_amount_frame:
        ax1.text(node.x - dx, node.y - dy, str("%.0f" % i), color='blue')
    elif i > len(nodes) - (frag_amount_v_under_arch + 1) * (frag_amount_v_under_arch + 1) - 1:
        ax1.text(node.x - dx, node.y + dy, str("%.0f" % i), color='red')
    else:
        dx = 0.5
        dy = 0.5
        if node.x < 0:
            dx *= -1
        if node.y < 0:
            dy *= -1
        ax1.text(node.x + dx/5, node.y + dy/5, str("%.0f" % i), color='black')

plt.show()

# 1) symmetry bottom border hv
# Normal solution of LCP in 57 steps
# Time:  54.269999504089355
# Right Border:30m
# Maximum lemke zn:2496589.0097933696,
# Minimum u_linear-3.7550786252912123,
# Minimum u_contact:-3.850372123179188
# SM rank:3283, rows:3283, cond:53230498148.392235 smT-sm max:1.9073486328125e-06, min:-1.9073486328125e-06

# 2) symmetry bottom border v only
# Normal solution of LCP in 49 steps
# Time:  22.337992429733276
# Right Border:30m
# Maximum lemke zn:4090641.441111044,
# Minimum u_linear-3.77568821268562,
# Minimum u_contact:-3.8654057975860923


