### Use the FEM_P2_solver class here ###
from FEM_P2 import FEM_P2_solver
import numpy as np
import sympy as sp

print("\n")
print("Solving Partial Differential Equation")
x,C,D, u, d,dx = sp.symbols('x C D u d dx')
print(-d**2*u/dx**2)
print("=", 2*x-1)
print("with solution")
p3 = -sp.Rational(1,3)*x**3 + sp.Rational(1,2)*x**2 + C*x + D-C-sp.Rational(1,6)
print(p3)
print("applying the Finite Element Method. \n")


C = float(input("Input Neumann Boundary Condition, C:"));
D = float(input("Input Dirichlet Boundary Condition, D:"));
Ne = int(input("Choose number of elements (int):"));

x = np.linspace(0,1, 10000)
pdenum = FEM_P2_solver(C,D,Ne)
pdenum.get_numerical_solution(x)

print("Press 1 to plot");
print("Press 2 to run a convergence test")
arg = input("Input number:")

if int(arg) == 1:
    pdenum.plot()
elif int(arg) == 2:
    N = int(input("Input number of experiments (int):"))

    # setup arrays to store h and l2
    l2_mesh = np.zeros(N); h_mesh = np.zeros(N)

    for i in range(N):
        norm = FEM_P2_solver(C,D,Ne)
        norm.get_numerical_solution(x)
        l2,h = norm.L2_norm()

        # Update Ne
        Ne = 4*Ne

        # Store l2 and h
        l2_mesh[i] = l2
        h_mesh[i] = h

    # Calculate convergence rates
    r = np.zeros(N)
    r = [np.log(l2_mesh[i]/l2_mesh[i-1])/np.log(h_mesh[i]/h_mesh[i-1]) \
    for i in range(1,N)]
    # print out table
    print('h            l2-norm     convergence rate')
    print("{:.4e}   {:.4e}  ---------".format(h_mesh[0],l2_mesh[0],))

    for i in range(1,N):
        print("{:.4e}   {:.4e}  {:.4e}".format(h_mesh[i],l2_mesh[i],r[i-1]))
