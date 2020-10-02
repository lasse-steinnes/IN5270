#### Main programme for running wave project ####
import numpy as np
from sympy import exp, sin, cos, symbols, simplify, diff, lambdify
from sympy.utilities.lambdify import lambdastr

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

Nx = 100; Ny = 100; T = 3; dt = 0.1
Lx = 10; Ly = 10; b = 0;


"""
 make some input options here:
 """
do = input('What do you want to do with the 2D Wavesolver \n (solve/test_constant/converge_standing/converge_damped/plotwaves)?')

if do == 'solve':
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt,'none')       # Initialize and create instance
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
    Nx = 5; Ny = 5; T = 5; Lx = 10;
    Ly = 10; b = 0; dt = 0.5

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
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt, 'standing')
    r = solver.test_convergence_rate(n_exp,u_e,I,b,f,V)

if do == 'converge_damped':
        # Define constants:  A,B, omega and d here.
    dt = 0.5
    mx = 1; my = 1;
    kx = mx*np.pi/Lx; ky = my*np.pi/Ly;
    b = 0.4; A = 1 ; B = 1 ;
    omega = 0.2
    d = 2
    x, y, t = symbols('x y t')
    g = x
    u_e = (A*cos(omega*t) + \
        B*sin(omega*t))*exp(-d*t)*cos(kx*x)*cos(ky*y)
    u_e_yy = simplify(diff(diff(u_e,y),y))
    u_e_xx = simplify(diff(diff(u_e,x),x))
    gu_x_x = simplify(diff(g*diff(u_e,x),x))
    gu_y_y = simplify(diff(g*diff(u_e,y),y))
    u_e_t = simplify(diff(u_e,t))
    u_e_tt = simplify(diff(diff(u_e,t),t))

    Vt = lambdify((x,y,t),u_e_t,'numpy') # modules = ['numpy']) # Maybe not need to lambdify this
    It = lambdify((x,y,t),u_e,'numpy') #, modules = ['numpy'])
    fsym = simplify(u_e_tt + b*u_e_t - gu_x_x - gu_y_y)

    # define functions
    V = lambda x,y: Vt(x,y,0)
    I = lambda x,y: It(x,y,0)
    g = lambdify((x,y),g,'numpy')
    f = lambdify((x,y,t),fsym, 'numpy')
    u_e = lambdify((x,y,t),u_e, 'numpy')

    n_exp = int(input('Number of experiments (int):'))
    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt,'damped')
    r = solver.test_convergence_rate(n_exp,u_e,I,b,f,V)

if do == 'plotwaves':
    # Define parameters
    b = 0;

    B0 = 1; Ba = B0*2; Bs = Lx*0.3; Bmx = 0.5*Lx; Bmy = 0.5*Ly; b = 1.1;
    H0 = 3*B0; I0 = 0; Ia = 6*B0; g_a = 9.81;

    speeds = input('Choose subsea surface (smooth/steep/rectangle)')
    if speeds == 'smooth':
        def g(x,y):
            B = I0 + Ia*np.exp(-((x-Bmx)/Bs)**2 -(y-Bmy/(b*Bs))**2)
            speed = g_a*(H0 - B)
            return speed

    elif speeds == 'steep':
        def g(x,y):
            B = B0 + Ba*np.cos((np.pi*x-Bmx)/(2*Bs))*np.cos((np.pi*y-Bmy)/(2*Bs))
            speed = g_a*(H0-B)
            return speed

    elif speeds == 'rectangle':
        def g(x,y):
            if Bmx - Bs <= x <= Bmx + Bs  and Bmy - b*Bs <= y <= Bmy + b*Bs:
                B = B0 + Ba
            else:
                B = B0
            speed = g_a*(H0-B)
            return speed

    # Define functions
    V = lambda x,y: 0
    I = lambda x,y: I0 +  Ia*np.exp(-((x-Bmx)/Bs)**2)

    solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt,'plot')       # Initialize and create instance
    u, t = solver.solve(I,b,f,V)                             # Solve the equation
