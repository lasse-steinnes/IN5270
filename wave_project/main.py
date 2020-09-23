#### Main programme for running wave project ####
import numpy as np
#from visualize import Visualize
from wave2D import Wave_2D_solver

def f(x,y,t):
    return 0

def g(x,y):
    return x

def I(x,y):
    return y

def V(x,y):
    return x

Nx = 40; Ny = 40; T = 10; dt = 0.1
Lx = 10; Ly = 10; b = 1;


solver = Wave_2D_solver(Lx,Ly,Nx,Ny,T,dt)       # Initialize and create instance
solver.set_inital(I)                           # Set t = 0
u, t = solver.solve(b,g,f,V)                    # Solve the equation
