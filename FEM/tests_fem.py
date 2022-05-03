import time
from FEM.scheme import NodeContainer, Node
from FEM.element_4node import Element4NodeLinearContainer

tests_passed = 0
tests_failed = 0


def decorator_info(func):
    """
    Decorator for printing information about finished tests
    :param func: test function
    :return:
    """
    def wrap(tests_passed, tests_failed):
        start = time.time()
        passed = func
        if passed:
            tests_passed += 1  # TODO: tests
        else:
            tests_failed += 1
        end = time.time()
        time_spent = end - start
        print(f'{func.__name__} - OK, time:{time_spent}')
    return wrap(tests_passed, tests_failed)


nodes = NodeContainer()


@decorator_info
def create_node_container():
    assert isinstance(nodes, NodeContainer)
    return True
#create_node_container()


@decorator_info
def node_in_scheme():
    x = 1
    y = 1
    nodes.add_node(x, y)
    assert isinstance(nodes[0], Node)
    node = nodes[0]
    assert node.x == x and node.y == y and node.rotation == None
    assert node.indices == [0, 1]
    node.rotation = True
    return True

print('tests_passed:', tests_passed)
print('tests failed:', tests_failed)




