from FEM.element_null import ElementNullContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_frame import ElementFrameContainer
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from SchemeForm.macro_element import ElementMacroContainer
from LCP.initial_table import InitialTable
from LCP.lemke import Lemke
from GUI.tkinter_gui import ContactFEM
from Visualize.plot_data_scheme import PlotScheme


from matplotlib import pyplot as plt
import math
import numpy as np
import time

# TIME
start = time.time()

nodes = NodeContainer()
nodes.add_node(0, -2)

nodes.add_node(0, 0)
nodes.add_node(4, 0)
nodes.add_node(4, 4)
nodes.add_node(1, 3)

nodes.add_node(7, 1)
nodes.add_node(7, 3)

nodes.add_node(1, -2)

E = 1e8
mu = 0.3
t = 1
Ef = 1e12
Af = 0.05
If = 1e-4

element_4node = Element4NodeLinearContainer(nodes_scheme=nodes)
element_4node.add_element(EN=[0, 7, 2, 1], E=1, mu=0.1, t=3)
#
element_frame = ElementFrameContainer(nodes_scheme=nodes)
element_frame.add_element(EN=[7, 5], E=Ef, A=Af, I=If)
element_frame.add_element(EN=[0, 1], E=Ef, A=Af, I=If)

element_null = ElementNullContainer(nodes_scheme=nodes)
# element_null.add_element(EN=[1, 3], cke=1, alpha=math.pi * 3 / 4)

macro_element = ElementMacroContainer(nodes_scheme=nodes)
macro_element.add_element(EN=[1, 2, 3, 4], frag_size=3, frag_size_h=1, frag_size_v=1,
                          E=E, mu=mu, t=t)
macro_element.add_element(EN=[2, 5, 6, 3], frag_size=3, frag_size_h=0.5, frag_size_v=0.7,
                          E=E, mu=mu, t=t, stitch=True)
# macro_element.add_element(EN=[1, 2, 3, 4], frag_size=0.2, frag_size_h=0.3, frag_size_v=1,
#                           E=E, mu=mu, t=t)
# macro_element.add_element(EN=[2, 5, 6, 3], frag_size=0.2, frag_size_h=0.7, frag_size_v=0.7,
#                           E=E, mu=mu, t=t)

macro_element.fragment_all(element_4node, element_frame, element_null)

for i, node in enumerate(nodes):
    print(f'№{i} h:{node.indices[0]}, v:{node.indices[1]}')
# TIME
end = time.time()
print(f'Time spend: {end - start} seconds')
# _--------------------------------------------------------------------

for i, node in enumerate(nodes.nodes_list):
    plt.scatter(node.x, node.y)
    plt.text(node.x, node.y, str("%.0f" % i), color='black')

def _modify_coordinates_for_plotting_4node(x, y):
    """
    Moving lines of sides of the element to the center for good looking
    :param x: vector of horizontal x coordinates
    :param y: vector of vertical y coordinates
    :return: x, y new vectors
    """
    avg_x = np.mean(x)
    avg_y = np.mean(y)
    x = [i + (avg_x - i) / 50 for i in x]
    y = [i + (avg_y - i) / 50 for i in y]
    return x, y

for element in element_4node.elements_list:
    x, y = element.nodes_coordinates(nodes)
    x, y = _modify_coordinates_for_plotting_4node(x, y)
    x.append(x[0])  # to encircle the element
    y.append(y[0])  # to encircle the element
    plt.plot(x, y, color='grey', zorder=1)

for element in element_frame:
    x, y = element.nodes_coordinates(nodes)
    plt.plot(x, y, color='blue', zorder=1)



plt.show()
