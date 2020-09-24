import numpy as np
import matplotlib.pyplot as plt

class Wave_2D_solver():
    """
    Provides a solver to various 2D wave equations with Neumann
    boundary conditions (BC's) for a square domain (0,Lx)x(0,Ly).
    The method used is a finite (centered difference scheme.
    """
    def __init__(self, Lx, Ly, Nx, Ny, T, dt):
        """
        Initalizes parameters
        - g(x,y)
        - f(x,y)
        """
        self.Nx, self.Ny = Nx, Ny                                     # x,y physical points only
        self.x = np.linspace(0,Lx,Nx+1);                # +1 because then we can use index Nx for last element
        self.y = np.linspace(0,Ly,Ny+1);                # same here
        self.dx = self.x[1] - self.x[0];
        self.dy = self.y[1] - self.y[0];                          # Spacing homogenous
        Nt = int(round(T/float(dt)))                    # Integer number of time iterations
        self.t = np.linspace(0,Nt*dt,Nt+1)              # time mesh
        self.dt = dt

        # Initialize vectors (2 dimensional mesh [i,j]) using ghost cells i = -1,0 ... Nx+1 and so on
        self.u_nn = np.zeros((Nx+3,Ny+3))           # time n-1, so that t -dt
        self.u_n = np.zeros((Nx+3,Ny+3))            # time n, so that
        self.u = np.zeros((Nx+3,Ny+3))           # time n+1, so that t + dt

        # Index sets, only iterates over INNER points
        u_nn, u_n, u = self.u,self.u_n, self.u
        self.Ix = range(1,u.shape[0]-1); self.Iy = range(1,u.shape[1]-1);
        self.It = range(0,self.t.shape[0])

        # Note: x[i-Ix[0]] Is the right index. Same for y

    def set_inital(self,I): # Set inital condition for t = 0
        Ix = self.Ix; x = self.y; Iy = self.Iy; y = self.y
        for i in Ix:
            for j in Iy:
                ix = i-Ix[0]; jy = j-Iy[0]
                self.u_n[i,j] = I(x[ix],y[jy])

        # Boundary condition which enforces Neuman condition du/dn = 0
        i = Ix[0]; self.u_n[i-1,:] = self.u_n[i+1,:]         # x = 0
        i = Ix[-1]; self.u_n[i+1,:] = self.u_n[i-1,:]        # x = Lx
        j = Iy[0]; self.u_n[:,j-1] = self.u_n[:,j+1]         # y = 0
        j = Iy[-1]; self.u_n[:,j+1] = self.u_n[:,j-1]        # y = Ly


    def var_coeff(self,g,x,y,i,j):
            dx, dy = self.dx, self.dy; u_n = self.u_n; Ix, Iy = self.Ix, self.Iy
            ix = i-Ix[0]; jy = j-Iy[0]

            if ix < self.Nx:
                gx_rhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix+1],y[jy]));
            else:
                gx_rhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix-1],y[jy]))

            if jy < self.Ny:
                gy_rhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix],y[jy+1]))
            else:
                gy_rhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix],y[jy-1]))

            gx_lhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix-1],y[jy]))
            gy_lhalf = 1/2*(g(x[ix],y[jy]) + g(x[ix],y[jy-1]))

            self.g_ux_x = 1/dx**2*(gx_rhalf*(u_n[i+1,j]-u_n[i,j]) \
                            - gx_lhalf*(u_n[i,j]-u_n[i-1,j]))
            self.gy_u_y = 1/dy**2*(gy_rhalf*(u_n[i,j+1]-u_n[i,j]) \
                            - gy_lhalf*(u_n[i,j]-u_n[i,j-1]))

    def set_boundary(self):
        Ix = self.Ix;  Iy = self.Iy;
        # Boundary condition which enforces Neuman condition du/dn = 0
        i = Ix[0]; self.u[i-1,:] = self.u[i+1,:]         # x = 0
        i = Ix[-1]; self.u[i+1,:] = self.u[i-1,:]        # x = Lx
        j = Iy[0]; self.u[:,j-1] = self.u[:,j+1]         # y = 0
        j = Iy[-1]; self.u[:,j+1] = self.u[:,j-1]        # y = Ly

    def advance(self,n,b,g,f,V):
        """
        Advances the solver
        """
        Ix = self.Ix; x = self.y; Iy = self.Iy; y = self.y
        t = self.t
        # first step(t = 1), <-- u_t initial condition V
        if n==1:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(g,x,y,i,j)
                    g_ux_x = self.g_ux_x; gy_u_y = self.gy_u_y

                    self.u[i,j] = 1/(2*b*self.dt)*(4*self.u_n[i,j] - 2*self.dt*V(x[ix],y[jy])*(b*self.dt-2)) \
                            + (self.dt/b)*(g_ux_x + gy_u_y + f(x[ix],y[jy],t[n]))
        # All other times steps
        else:
            for i in Ix:
                for j in Iy:

                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(g,x,y,i,j)
                    g_ux_x = self.g_ux_x; g_uy_y = self.gy_u_y

                    self.u[i,j] = 1/(2 + b*self.dt)*(self.u_nn[i,j]*(b*self.dt -2) \
                    + 4*self.u_n[i,j] + (2*self.dt**2)*(g_ux_x + g_uy_y + f(x[ix],y[jy],t[n])))

        self.set_boundary()

    def advance_b0(self,n,g,f,V):
        """
        Advances the solver for b = 0
        """
        Ix = self.Ix; x = self.y; Iy = self.Iy; y = self.y
        t = self.t
        # first step(t = 1), <-- u_t initial condition V
        if n==1:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(g,x,y,i,j)
                    g_ux_x = self.g_ux_x; gy_u_y = self.gy_u_y
                    self.u[i,j] = self.u_n[i,j] - V(x[ix],y[jy])*self.dt \
                                +(self.dt**2)/2*(g_ux_x + gy_u_y + f(x[ix],y[jy],t[n]))
        # All other times steps
        else:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(g,x,y,i,j)
                    g_ux_x = self.g_ux_x; g_uy_y = self.gy_u_y
                    self.u[i,j] = 2*self.u_n[i,j] - self.u_nn[i,j] + \
                                self.dt**2*(g_ux_x + g_uy_y + f(x[ix],y[jy],t[n]))
        self.set_boundary()


    def solve(self,I,b,g,f,V):
        # Make solver for all time steps
        self.set_inital(I)

        if b != 0:
            for n in self.It[1:-1]:
                self.advance(n,b,g,f,V)
                self.u_nn, self.u_n, self.u = self.u_n, self.u, self.u_nn  # Update time vectors

        else:
            for n in self.It[1:-1]:
                self.advance_b0(n,g,f,V)
                self.u_nn, self.u_n, self.u = self.u_n, self.u, self.u_nn  # Update time vectors

        return self.u, self.t
