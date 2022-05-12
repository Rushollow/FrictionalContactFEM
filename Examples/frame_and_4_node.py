import comtypes.client as com  # for mathcad tests

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
MATHCAD_OPEN = False

######### Mathcad TEST
if MATHCAD_OPEN:
    mathcad15 = com.CreateObject("Mathcad.Application")
    mathcad15.Visible = True
    WS = mathcad15.Worksheets
    WK = WS.Open(r"C:\1MyFiles\Python\Arch\test.xmcd")
#####################
###### EXAMPLE:
###### b = WK.GetValue('v_output').Real
###### WK.SetValue("val1", val1)

# add nodes # for 4 node element
nodes = NodeContainer()
nodes.add_node(0, 0)
nodes.add_node(1.75, 0)
nodes.add_node(1.75, 1.75)
nodes.add_node(0, 1.75)
# for frame
nodes.add_node(0, 1.75)
nodes.add_node(1.75, 1.75)
nodes.add_node(2.5, 2.75)
nodes.add_node(4, 2.75)
# for 2nd 4node element
nodes.add_node(2.5, 2.75)
nodes.add_node(4, 2.75)
nodes.add_node(4, 4)
nodes.add_node(2.5, 4)
XY1, XY2 = [], []
i = 0
for node in nodes:
    XY1.append(node.x)
    XY2.append(node.y)
    i += 1
if MATHCAD_OPEN:
    WK.SetValue("XY1", XY1)
    WK.SetValue('XY2', XY2)

# set inputs
Ar = 0.2
Er = 2e11
Ix = 2.3 * 10**(-4)
E = 5e7
F = 1e6
gamma = 5
t = 1
mu = 0.3

# add elements
element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_4node.add_element(EN=[0, 1, 2, 3], E=E, mu=mu, t=t)
element_4node.add_element(EN=[8, 9, 10, 11], E=E, mu=mu, t=t)

# add frame elements
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[4, 5], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[5, 6], E=Er, A=Ar, I=Ix)
element_frame.add_element(EN=[6, 7], E=Er, A=Ar, I=Ix)
# n null elements and adding t null elements silently
element_null = ElementNullContainer(nodes_scheme=nodes)
element_null.add_element(EN=[3, 4], cke=123, alpha=math.pi / 2)
element_null.add_element(EN=[2, 5], cke=123, alpha=math.pi / 2)
element_null.add_element(EN=[6, 8], cke=123, alpha=math.pi / 2)
element_null.add_element(EN=[7, 9], cke=123, alpha=math.pi / 2)

# form R, RF and solve SLAE
sm = StiffnessMatrix(nodes=nodes, el_frame=element_frame, el_4node=element_4node, el_null=element_null)
nodes_to_support = [0, 1, 10, 11]
sm.support_nodes(nodes_to_support, direction='hv')
lv = LoadVector()
lv.add_concentrated_force(force=F, degree_of_freedom=14)
lv.add_concentrated_force(force=-F, degree_of_freedom=8)
lv.add_concentrated_force(force=F * 1.4, degree_of_freedom=13)
lv.add_concentrated_force(force=-F * 2, degree_of_freedom=11)


# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=lv, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=25, autorun=True)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

for i in range(len(graph.lemke.zn_anim)):
    print(f'{i} zn:{graph.lemke.zn_anim[i]} xn:{graph.lemke.zt_anim[i]}')

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
