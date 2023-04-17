from cvxopt import matrix, solvers
import numpy as np
# M = matrix([ [ .3, -.4,  -.2,  -.4,  1.3 ],
#                  [ .6, 1.2, -1.7,   .3,  -.3 ],
#                  [-.3,  .0,   .6, -1.2, -2.0 ] ])
# q = matrix([1.5, .0, -1.2, -.7, .0])
#
M = matrix([[4,2,0,3],
           [-1,4,-3,-6],
           [1,-1,1,1],
           [0,1,0,5]])
# q = matrix([-1,2,-1,-1])
q = matrix([1.5, -1.2, -.7, .0])

m, n = M.size
I = matrix(0.0, (n,n))
I[::n+1] = 1.0
G = matrix([-I, matrix(0.0, (1,n)), I])
h = matrix(n*[0.0] + [1.0] + n*[0.0])
dims = {'l': n, 'q': [n+1], 's': []}
x = solvers.coneqp(M.T * M, -M.T * q, G, h, dims)['x']
print(x)