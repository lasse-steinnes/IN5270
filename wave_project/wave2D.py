import numpy as np
import matplotlib.pyplot as plt

class Wave_2D_solver():
    """
    Provides a solver to various 2D wave equations with Neumann
    boundary conditions (BC's) for a square domain (0,Lx)x(0,Ly).
    The method used is a finite (centered difference scheme.
    """
    def __init__(self, b, g, f, I, V, Lx, Ly, Nx, Ny, T, dt):
        """
        Initalizes parameters
        - g(x,y)
        - f(x,y)
        """
        self.x = np.linspace(0,Lx,Nx+1);                # +1 because then we can use index Nx for last element
        self.y = np.linspace(0,Ly,Ny+1);                # same here
        self.dx = x[1] -x[0]; self.dy = y[1] - y[0];    # Spacing homogenous
        self.Nt = int(round(T/float(dt)))               # Integer number of time iterations
        self.t = np.linspace(0,Nt*dt,N+1)               # time mesh

        # Initialize vectors (2 dimensional mesh [i,j])
        u_nn = np.zeros((Nx+3,Ny+3))           # time n-1, t -dt
        u_n = np.zeros((Nx+3,Ny+3))            # time n, t
        u = np.zeros((Nx+3,Ny+3))           # time n+1, t + dt

        # Index sets, only iterates over INNER points
        Ix = range(1,u.shape[0]-1); Iy = range(1,u.shape[1]-1);
        It = range(0,t.shape[0])

        # Note: x[i-Ix[0]] Is the right index. Same for y

    def set_inital(self): # Set inital condition for t = 0
        for i in Ix:
            for j in Iy:
                u_n[i,j] = I(x[i-Ix[0],(y[j-Iy[0]]))

    def advance(self,step1 = False):
        """
        Advances the solver
        """
        if step1:                   # write first step here (t = 1)
        for i in Ix:
            for j in Iy:
                u[i,j] =

        else:
            for i in Ix:
                for j in Iy:
                    # write other steps here
        return u


        # Boundary condition which enforces Neuman condition du/dn = 0
        i = Ix[0]; u[i-1,:] = u[i+1,:]         # x = 0
        i = Ix[-1]; u[i+1,:] = u[i-1,:]        # x = Lx
        j = Iy[0]; u[:,j-1] = u[:,j+1]         # y = 0
        j = Iy[-1]; u[:,j+1] = u[:,j-1]        # y = Ly
