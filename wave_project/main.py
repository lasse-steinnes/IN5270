#### Main programme for running wave project ####
import numpy as np
#from visualize import Visualize
from wave2D import Wave_2D_solver
from test_funcs import TestClass

def f(x,y,t):
    return 0

def g(x,y):
    return x

def I(x,y):
    return 1

def V(x,y):
    return x

Nx = 4; Ny = 4; T = 1; dt = 0.1
Lx = 1; Ly = 1; b = 0;


"""
 make some input options here:
 """
do = input('What do you want to do with the 2D Wavesolver \n (solve/test_constant/test_converge)?')

if do == 'solve':
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt)       # Initialize and create instance
    u, t = solver.solve(I,b,f,V)                  # Solve the equation

if do == 'test_constant':
    test = TestClass()
    test.test_constant_solution()

# Want to calculate convergence rates
if do == 'test_converge':
    A = 1; mx = 1; my = 1; omega = 0.2;
    k_x = mx*np.pi/Lx; k_y = my*np.pi/Ly;
    b = 0

    def u_e(x,y,t):
        u_emat = np.zeros((x.shape[0],y.shape[0]))

        for i in range(len(x)):
            for j in range(len(y)):
                u_emat[i,j] = A*np.cos(k_x*x[i])*np.cos(k_y*y[j])*np.cos(omega*t)
        return u_emat

    I = lambda x,y: A*np.cos(k_x*x)*np.cos(k_y*y)
    V = lambda x,y: 0

    n_exp = int(input('Number of experiments:'))
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt, test_convergence = True)
    r = solver.test_convergence_rate(n_exp,u_e,I,b,f,V)

print('Note: You can also test the implementation of constant solution by calling python -q test_funcs.py')
