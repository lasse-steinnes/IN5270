#### Main programme for running wave project ####
import numpy as np
import sympy as sym

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

Nx = 5; Ny = 5; T = 1; dt = 0.1
Lx = 10; Ly = 10; b = 0;


"""
 make some input options here:
 """
do = input('What do you want to do with the 2D Wavesolver \n (solve/test_constant/converge_standing/converge_damped)?')

if do == 'solve':
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt)       # Initialize and create instance
    u, t = solver.solve(I,b,f,V)                  # Solve the equation

if do == 'test_constant':
    test = TestClass()
    test.test_constant_solution()
    print('Note: You can also test the implementation of constant solution by calling python -q test_funcs.py')

# Want to calculate convergence rates
if do == 'converge_standing':
    A = 1; mx = 1; my = 1;
    k_x = mx*np.pi/Lx; k_y = my*np.pi/Ly;
    b = 0; omega = np.sqrt(k_x**2 + k_y**2)  # for g = 1
    dt = 0.5

    def u_e(x,y,t):
        u_emat = np.zeros((x.shape[0],y.shape[0]))

        for i in range(len(x)):
            for j in range(len(y)):
                u_emat[i,j] = A*np.cos(k_x*x[i])*np.cos(k_y*y[j])*np.cos(omega*t)
        return u_emat

    I = lambda x,y: A*np.cos(k_x*x)*np.cos(k_y*y)
    V = lambda x,y: 0
    g = lambda x,y: 1
    f = lambda x,y,t: 0

    n_exp = int(input('Number of experiments (int):'))
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt, test_convergence = True)
    r = solver.test_convergence_rate(n_exp,u_e,I,b,f,V)

if do == 'converge_damped':
    b, A, B, d, kx, ky, omega, x, y, t = sym.symbols('b, A B d kx ky omega x y t')
    g = x**2
    u_e = (A*sym.cos(omega*t) + \
        B*sym.sin(omega*t))*sym.exp(-d*t)*sym.cos(kx*x)*sym.cos(ky*y)
    u_e_yy = sym.simplify(sym.diff(sym.diff(u_e,y),y))
    u_e_xx = sym.simplify(sym.diff(sym.diff(u_e,x),x))
    gu_x_x = sym.simplify(sym.diff(g*sym.diff(u_e,x),x))
    gu_y_y = sym.simplify(sym.diff(g*sym.diff(u_e,y),y))
    u_e_t = sym.simplify(sym.diff(u_e,t))
    u_e_tt = sym.simplify(sym.diff(sym.diff(u_e,t),t))

    Vt = sym.lambdify([x,y,t],u_e_t, modules = ['numpy'])
    It = sym.lambdify([x,y,t],u_e, modules = ['numpy'])


    def V(x,y):
        return Vt(x,y,0)

    def I(x,y):
        return It(x,y,0)

    div = input('derivate (yes/no):')
    if div == 'yes':
        print('\n')
        #print('du/dt:', u_e_t,'\n')
        #print('d^2u/dy^2:', u_e_yy,'\n')
        #print('d^2u/dx^2:', u_e_xx,'\n')
        #print('d^2u/dt^2:', u_e_tt,'\n')
        #print('d/dx (g du/dx):', gu_x_x,'\n')
        #print('d/dy (g du/dy):', gu_y_y,'\n')
        sym.solveset(u_e_tt + b*u_e_t - gu_x_x - gu_y_y,omega)
        #print('waveeq:', sym.simplify(u_e_tt + b*u_e_t),'\n')
        #print('=', sym.simplify(gu_x_x + gu_y_y))

#    A = 1; mx = 1; my = 1; omega = 0.2;
#    k_x = mx*np.pi/Lx; k_y = my*np.pi/Ly;
