import numpy as np
import matplotlib.pyplot as plt
from visualize import plot_u

class Wave_2D_solver():
    """
    Provides a solver to various 2D wave equations with Neumann
    boundary conditions (BC's) for a square domain (0,Lx)x(0,Ly).
    The method used is a finite centered difference scheme.
    """
    def __init__(self, g, Lx, Ly, Nx, Ny, T, dt, do):
        """
        Initalizes parameters
        - If test_convergence, compute error and test convergence rates
            for a constant g
        """
        self.do = do
        self.g = g

        self.T = T; self.Lx = Lx; self.Ly = Ly
        self.Nx, self.Ny = Nx, Ny                       # x,y physical points only
        self.x = np.linspace(0,Lx,Nx+1);                # +1 because then we can use index Nx for last element
        self.y = np.linspace(0,Ly,Ny+1);                # same here
        self.dx = self.x[1] - self.x[0];
        self.dy = self.y[1] - self.y[0];                # Spacing homogenous
        Nt = int(round(T/float(dt)))                    # Integer number of time iterations
        self.t = np.linspace(0,Nt*dt,Nt+1)              # time mesh
        self.dt = dt

        # Initialize vectors (2 dimensional mesh [i,j]) using ghost cells i = -1,0 ... Nx+1 and so on
        self.u_nn = np.zeros((self.Nx+3,self.Ny+3))           # time n-1, so that t -dt
        self.u_n = np.zeros((self.Nx+3,self.Ny+3))            # time n, so that
        self.u = np.zeros((self.Nx+3,self.Ny+3))           # time n+1, so that t + dt

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

        ## parameters for converge test initalized here for convenience
        if not self.do == 'damped' and not self.do =='standing':
            self.error = 0; self.residual = 0; self.u_e = lambda x,y,t: 0

    def var_coeff(self,x,y,i,j):
            g = self.g

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

    def advance(self,n,b,f,V):
        """
        Advances the solver
        """
        Ix = self.Ix; x = self.y; Iy = self.Iy; y = self.y
        t = self.t; g = self.g;
        # first step(t = 1), <-- u_t initial condition V
        if n==1:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(x,y,i,j)
                    g_ux_x = self.g_ux_x; gy_u_y = self.gy_u_y

                    self.u[i,j] = 1/4*(4*self.u_n[i,j] - 2*self.dt*V(x[ix],y[jy])*(b*self.dt-2) \
                           + (2*self.dt**2)*(g_ux_x + gy_u_y + f(x[ix],y[jy],t[n])))

        # All other times steps
        else:
            for i in Ix:
                for j in Iy:

                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(x,y,i,j)
                    g_ux_x = self.g_ux_x; g_uy_y = self.gy_u_y

                    self.u[i,j] = 1/(2 + b*self.dt)*(self.u_nn[i,j]*(b*self.dt -2) \
                    + 4*self.u_n[i,j] + (2*self.dt**2)*(g_ux_x + g_uy_y + f(x[ix],y[jy],t[n])))

        self.set_boundary()

    def advance_b0(self,n,f,V):
        """
        Advances the solver for b = 0
        """
        Ix = self.Ix; x = self.y; Iy = self.Iy; y = self.y
        t = self.t; g = self.g;
        # first step(t = 1), <-- u_t initial condition V
        if n==1:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(x,y,i,j)
                    g_ux_x = self.g_ux_x; gy_u_y = self.gy_u_y
                    self.u[i,j] = self.u_n[i,j] - V(x[ix],y[jy])*self.dt \
                                +(self.dt**2)/2*(g_ux_x + gy_u_y + f(x[ix],y[jy],t[n]))

        # All other times steps
        else:
            for i in Ix:
                for j in Iy:
                    ix = i-Ix[0]; jy = j-Iy[0]
                    self.var_coeff(x,y,i,j)
                    g_ux_x = self.g_ux_x; g_uy_y = self.gy_u_y
                    self.u[i,j] = 2*self.u_n[i,j] - self.u_nn[i,j] + \
                                self.dt**2*(g_ux_x + g_uy_y + f(x[ix],y[jy],t[n]))
        self.set_boundary()


    def solve(self,I,b,f,V):
        # Make solver for all time steps
        g = self.g
        self.set_inital(I)

        if self.do == 'plot':
            fig = True
        else:
            fig = False

        if b != 0:
            self.error = 0
            for self.n in self.It[1:-1]:
                n = self.n
                self.advance(n,b,f,V)
                self.residual = np.subtract(self.u_e(self.x,self.y,self.t[self.n]),self.u[1:-1,1:-1])
                self.error = max(self.error,np.abs(self.residual.max()))

                plot_u(self.u[1:-1], self.x,self.y, self.t,n,fig, save_plot = False)
                self.u_nn, self.u_n, self.u = self.u_n, self.u, self.u_nn  # Update time vectors

        else:
            self.error = 0
            for self.n in self.It[1:-1]:
                n = self.n
                self.advance_b0(n,f,V)
                self.residual = np.subtract(self.u_e(self.x,self.y,self.t[self.n]),self.u[1:-1,1:-1])
                self.error = max(self.error,np.abs(self.residual.max()))

                plot_u(self.u[1:-1], self.x,self.y, self.t,n,fig, save_plot = False)
                self.u_nn, self.u_n, self.u = self.u_n, self.u, self.u_nn  # Update time vectors
        return self.u[1:-1,1:-1], self.t


    def test_convergence_rate(self,n_experiments,u_e,I,b,f,V):
        # write code to do a convergence test here
        """
        compare numerical and an analytically exact solution
        (truncation and round off error in numerical solution)
        to see if converge == 2 for finite difference method
        input u_e
        """
        self.Es = []                                      # initialize error
        self.hs = []
        self.u_e = u_e
        self.initialize_convergence()
        for i in range(n_experiments):
            self.error = 0
            self.solve(I,b,f,V)
            # calulate errors
            self.Es.append(self.error)
            self.hs.append(self.dt)

            ### initialize correct in each case ###
            self.dt /= 2
            self.initialize_convergence()
        print('E',self.Es)                                                 # half time step for next simulation
        # now get convergence rate
        r = [np.log(self.Es[i]/self.Es[i-1])/np.log(self.hs[i]/self.hs[i-1]) \
            for i in range(1,n_experiments)]
        print('r',r)
        print('\nExpected converge rate: 2. tol: 0.05')
        for i in range(n_experiments-1):
            if abs(r[i] - 2) < 0.05:
                print('Expected convergence rate after experiment: {:d} , with h: {:.4e}'.format(i+1,self.hs[i]))
                break
            else:
                print('Not expected convergence rate after experiment: {:d}, with h: {:.4e}'.format(i+1,self.hs[i]))

        return r, self.Es, self.hs

    def initialize_convergence(self):
        g, Lx, Ly = self.g, self.Lx, self.Ly

        if self.do == 'standing':
            C = 0.6
            c_max = np.sqrt(g(1,1))

        elif self.do == 'damped':
            beta =  1
            C = float(0.3/beta)
            x_mock, y_mock = np.linspace(0,Lx,100), np.linspace(0,Ly,100)
            c_max = np.sqrt(np.max(g(x_mock,y_mock)))

        self.Nt = int(round(self.T/self.dt))
        self.t = np.linspace(0, self.Nt*self.dt, self.Nt+1)
        self.dx = self.dt*c_max/C
        self.dy = self.dt*c_max/C
        self.Nx = int(round(Lx/self.dx))
        self.Ny = int(round(Ly/self.dy))
        self.x = np.linspace(0, self.Nx*self.dx, self.Nx+1)
        self.y = np.linspace(0, self.Ny*self.dx , self.Ny+1)

        # print stability criteria
        print("stability criteria")

        print('dt:',self.dt,'dx:',self.dx,'dy',self.dy,'dt less than',1/c_max*1/(np.sqrt(1/self.dx**2 + 1/self.dy**2)))

        # Initialize vectors (2 dimensional mesh [i,j]) using ghost cells i = -1,0 ... Nx+1 and so on
        self.u_nn = np.zeros((self.Nx+3,self.Ny+3))           # time n-1, so that t -dt
        self.u_n = np.zeros((self.Nx+3,self.Ny+3))            # time n, so that
        self.u = np.zeros((self.Nx+3,self.Ny+3))           # time n+1, so that t + dt

        # Index sets, only iterates over INNER points
        u_nn, u_n, u = self.u,self.u_n, self.u
        self.Ix = range(1,u.shape[0]-1); self.Iy = range(1,u.shape[1]-1);
        self.It = range(0,self.t.shape[0])
