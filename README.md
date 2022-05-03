# Contact_FEM
Constructive nonlinearity structural mechanics FEM program  
OOP based FEM program with Frame, 4Node, null - elements.  
Contact problems solving with LCP - linear complementary problem  
work in progress...  
----------------------------
pycharm project  
Python v3.7  
Libraries:  
>   external:  
>>        Numpy (arrays, concatenation, linalg.solve...)  
>>        Matplotlib (fig, ax)  
>>        Tkinter
>    build in:  
>>        weakref (collecting instances of elements classes)  
>>        gc (garbage collector)  
>>        math (for math)
>>        tkinter (GUI)  

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



