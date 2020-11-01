import sympy as sp
import numpy as np

class FEM_P2_solver():
    """
    Provides an ODE solver using the finite element method with
    P2 elements.
    """
    def __init__(self,C,D,Ne):
        self.C = C          # Neumann condition
        self.D = D          # Dirichlet condition
        self.Ne = Ne        # Number of elements, each with three nodes.
        self.oth = 2        # P2 elements

        # setup interval and number of nodes
        n_nodes = 3*Ne - (Ne-1)
        b = 1; a = 0; interval = b-a

        # step size
        self.step =  interval/(n_nodes-1)

        # Get midpoint and nodes mesh
        self.nodes = [i*self.step for i in range(n_nodes)] # setup physical mesh point
        self.two_n_m = [self.nodes[i*self.oth+2]  + \
                        self.nodes[i*self.oth] for i in range(self.Ne)]
        self.nodes_m = 0.5*np.array(self.two_n_m)
        self.nodes = np.array(self.nodes)

        # Get step size for element endpoints for a uniform mesh
        self.h = self.step*self.oth

        # Set up element
        elements = [[i*self.oth,i*self.oth+1,i*self.oth+2] for i in range(self.Ne)]
        #print(self.nodes)
        #print(self.h)


    def setup_matrix(self):
        """
        sets up the element matrix
        """
        oth = self.oth
        self.A = np.zeros((oth + 1,oth + 1))
        self.A[0,0] = self.A[self.oth,oth] = 7;
        self.A[oth,1] = self.A[1,oth] = -8
        self.A[1,0] = self.A[0,1] = -8
        self.A[oth,0] = self.A[0,oth] = 1

        self.A = 1/(3*h)*self.A     # The final A matrix check
                                   # if constant can be factored out in the end
    def setup_matrixN(self):
        """
        Set up special element matrix to enforce Dirichlet condition
        """
        oth = self.oth
        self.A_n = np.zeros((oth + 1,oth + 1))
        self.A_n[0,0] = 7;
        self.A_n[1,oth] = -8
        self.A_n[1,0] = self.A_n[0,1] = -8
        self.A_n[0,oth] = 1

        # Special alteration
        self.A_n[oth,0] = 0; self.A_n[oth,1] = 0;
        self.A_n[self.oth,oth] = 3*h

        self.A_n = 1/(3*h)*self.A_n


    def setup_basis(self):
        """
        sets up the reference basis
        """

    def assemble(self):
        """
        Map to global system/assemble
        element matrices and element vectors
        """

    def find_b_ref(self):
         x, x_m, h, X = sp.symbols('x x_m h X')         # make symbols
         x = x_m + h/2*X                                # Apply transformation
         wrt = (X, -1, 1)                               # Interval and variable
         b_0 = sp.integrate(h*sp.Rational(1,4)*(x-sp.Rational(1,2))*(X-1)*X, wrt)   # b0
         b_1 = sp.integrate(h*sp.Rational(1,4)*(x-sp.Rational(1,2))*(1-X**2), wrt)  # b1
         b_2 = sp.integrate(h*sp.Rational(1,4)*(x-sp.Rational(1,2))*(X+1)*X, wrt)   # b2

         print("Integrating with numerical python element vector in reference")
         print("b0",b_0)
         print("b1",b_1)
         print("b2",b_2)

         self.b0_e = sp.lambdify((x_m),b_0,'numpy')
         self.b1_e = sp.lambdify((x_m),b_1,'numpy')
         self.b2_e = sp.lambdify((x_m),b_2,'numpy')

    def b_first(self):
        self.b0first_e = self.b0_e - self.C

    def b_last(self):
        self.b2last_e = self.D

    def solve(self):
        """
        solves the assembled system
        """
    def plot(self):
        """
        Plot the numerical and analytical solution
        """
