# Contact_FEM
Constructive nonlinearity structural mechanics FEM program  
OOP based FEM program with Frame, 4Node, null - elements.  
Contact problems solving with LCP - linear complementary problem  
work in progress...  
----------------------------
Python v3.10 
Libraries:  
>   external:  
>>        Numpy (arrays, concatenation, linalg.solve...)  
>>        PyQt 5 v5.15.6 (for UI)  
>>        pyqtgraph (for graphs and animations)
>>        XlsxWriter 3.0.3 (for writing LCP tables to Excel if you need)
>>        prettytable 3.3.0 (for some good looking tables in the console)
>    build in:  
>>        weakref (collecting instances of elements classes)  
>>        gc (garbage collector)  
>>        math (for math)

Usage:
> 1. Clone the repository
> 2. Install the dependencies (libraries)
> 3. Run any .py file in Examples folder for the test
> 4. You can edit input data global parameters in .../`input_data.py` and in .../`constants_data.py` files
> 5. After problem solved and interface is loaded press button `run` in the up left corner
> 6. Now you can zoom graphs and press LEFT and RIGHT to see the solutions steps in  UI

This is my PhD project so...
there will be no guides or documentation about this project

Also, I was experimenting with data transfer from python to MathCad 15. Some code in examples could be there.

---------------------------
Variables and expressions:  
R or r - stiffness matrix (R - global stiffness matrix)  
RF of rf - load vector  
U or u - displacement vector (the displacements by dof in the system that we get by solving SLAE)  
SLAE - system of linear algebraic equations  
MI - Matrix of Indices (element.MI)  
EN - Element Nodes (element.EN)  

dof - degree of freedom  
SSS - stress strain state  
N - longitudinal effort  
M - bending moment  



