# With variable __all__ in __init__.py file we can choose what modules we need to import from "some" directory
# Ex: We have directory "FEM" and packejes theme called: "Elements", "Stiffness_matrix", so we can import only "FEM"
# using __all__ = ["FEM"] and importing such way: "from FEM import *". If __all__ variable is default by importing
# "*" we will have import all modules.

__all__ = ["FEM"]

