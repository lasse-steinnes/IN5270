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

Nx = 40; Ny = 40; T = 10; dt = 0.1
Lx = 10; Ly = 10; b = 0;


solver = Wave_2D_solver(Lx,Ly,Nx,Ny,T,dt)       # Initialize and create instance
u, t = solver.solve(I,b,g,f,V)                  # Solve the equation

test = TestClass()
test.test_constant_solution()
