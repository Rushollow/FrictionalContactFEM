from FEM.element_null import ElementNullContainer
from FEM.element_4node import Element4NodeLinearContainer
from FEM.element_frame import ElementFrameContainer
from FEM.scheme import NodeContainer, StiffnessMatrix, LoadVector
from SchemeForm.macro_element import ElementMacroContainer

from Visualize.plot_data_qt import PlotScheme  # for visualizing
from GUI.PyQt.contactFEM import application
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

sm = StiffnessMatrix(nodes, el_frame=None, el_4node=element_4node, el_null=None)
# plot --------------------------------------------------------------------------
# Calculation and plotting object
graph = PlotScheme(nodes=nodes, sm=sm, lv_const=None, lv_variable=None,
                   element_frame=element_frame, element_container_obj=element_4node, element_null=element_null,
                   partition=10, scale_def=25, autorun=False)

# calculate time
end = time.time()
last = end - start
print("Time: ", last)

if __name__ == "__main__":
    graph.fill_arrays_scheme()  # form info for plot at UI
    application(graph)
