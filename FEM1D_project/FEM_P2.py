import sympy as sp
import numpy as np
import pentapy as pp
import matplotlib.pyplot as plt

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
        self.n_nodes = 3*Ne - (Ne-1)
        b = 1; a = 0; interval = b-a

        # step size
        self.step =  interval/(self.n_nodes-1)

        # Get midpoint and nodes mesh (nodes are physical mesh points)
        self.nodes = [i*self.step for i in range(self.n_nodes)] # setup physical mesh point
        self.two_n_m = [self.nodes[i*self.oth+2]  + \
                        self.nodes[i*self.oth] for i in range(self.Ne)]
        self.nodes_m = 0.5*np.array(self.two_n_m)
        self.nodes = np.array(self.nodes)

        # Get step size for element endpoints for a uniform mesh
        self.h = self.step*self.oth

        # Set up element
        self.elements = [[i*self.oth,i*self.oth+1,i*self.oth+2] for i in range(self.Ne)]
        self.local_nodes = len(self.elements[0])
        # Set up final matrix A, numerical solution u and right hand side b
        self.A = np.zeros((self.n_nodes,self.n_nodes))
        self.b = np.zeros(self.n_nodes);
        #print(self.A)#print(self.nodes_m)#print(self.h) #print(elements)
        #print(self.nodes_m)

    def setup_element_matrix(self):
        """
        sets up the element matrix
        """
        oth,h = self.oth,self.h;
        self.A_e = np.zeros((oth + 1,oth + 1))
        self.A_e[0,0] = self.A_e[self.oth,oth] = 7;
        self.A_e[1,1] = 16;
        self.A_e[oth,1] = self.A_e[1,oth] = -8
        self.A_e[1,0] = self.A_e[0,1] = -8
        self.A_e[oth,0] = self.A_e[0,oth] = 1

        self.A_e = 1/(3*h)*self.A_e     # The final A matrix check
                                   # if constant can be factored out in the end
    def setup_element_matrixN(self):
        """
        Set up special element matrix to enforce Dirichlet condition
        """
        oth,h = self.oth, self.h
        self.A_n = np.zeros((oth + 1,oth + 1))
        self.A_n[0,0] = 7;
        self.A_n[1,1] = 16;
        self.A_n[1,oth] = -8
        self.A_n[1,0] = self.A_n[0,1] = -8
        self.A_n[0,oth] = 1

        # Special alteration
        self.A_n[oth,0] = 0; self.A_n[oth,1] = 0;
        self.A_n[self.oth,oth] = 3*h

        self.A_n = 1/(3*h)*self.A_n



    def evaluate_num(self,x):
        """
        Evaluates the numerical solution at mesh points
        """
        ## Decide what local node a globale node lies within with use of modulo
        ## but here is the easier way, but maybe less efficient way chosen
        ## Set up basis functions
        self.psi = lambda x,a,b,c: (x-b)/(a-b)*(x-c)/(a-c);
        ## numerical solution
        self.u = np.zeros(len(x))
        j = 0; i = 0
        while j <= len(self.nodes)-1:
            while i < len(x) and x[i] <= self.nodes[j+2]:
                #print('j:',j)
                #print('xi:',x[i])
                #print('node:',self.nodes[j+2])
                self.u[i] = self.c[j]*self.psi(x[i],self.nodes[j],self.nodes[j+1],self.nodes[j+2]) + \
                self.c[j+1]*self.psi(x[i],self.nodes[j+1],self.nodes[j],self.nodes[j+2]) + \
                self.c[j+2]*self.psi(x[i],self.nodes[j+2],self.nodes[j],self.nodes[j+1])
                i += 1
            j += 2

    def assemble(self):
        """
        Map to global system/assemble
        element matrices and element vectors
        """
        local_nodes = self.local_nodes;

        # Integrate on referance
        self.find_b_ref();
        b0_e, b1_e, b2_e = self.b0_e,self.b1_e,self.b2_e
        b_e = [b0_e,b1_e,b2_e]

        # setup element matrix
        self.setup_element_matrix();

         # Assemble inner elements
        for e in range(1,len(self.elements)-1):
            for r in range(local_nodes):
                for s in range(local_nodes):
                     self.A[self.elements[e][r], self.elements[e][s]] += self.A_e[r,s]
                self.b[self.elements[e][r]] += b_e[r](self.nodes_m[e],self.h)

        ### Handle boundary conditions in first and last element ###
        # change rhs index 0 for Neumann boundary
        e = 0;
        for r in range(local_nodes):
            for s in range(local_nodes):
                self.A[self.elements[e][r], self.elements[e][s]] += self.A_e[r,s]
            if r != 0:
                self.b[self.elements[e][r]] += b_e[r](self.nodes_m[e],self.h)
            else:
                self.b[self.elements[e][r]] += b_e[r](self.nodes_m[e],self.h) \
                    - self.C


        # last element
        self.setup_element_matrixN();
        e = self.Ne-1;
        for r in range(local_nodes):
            for s in range(local_nodes):
                self.A[self.elements[e][r], self.elements[e][s]] += self.A_n[r,s]
            self.b[self.elements[e][r]] += b_e[r](self.nodes_m[e],self.h)

        # Alter last index right hand side
        r = local_nodes-1
        self.b[self.elements[e][r]] = self.D
        #print(self.b[self.elements[e][r]])


    def find_b_ref(self):
         x, x_m, h, X = sp.symbols('x x_m h X')         # make symbols
         x = x_m + h/2*X                                # Apply transformation
         wrt = (X, -1, 1)                               # Interval and variable
         b_0 = sp.integrate(h*sp.Rational(1,4)*(2*x-1)*(X-1)*X, wrt)   # b0
         b_1 = sp.integrate(h*sp.Rational(1,2)*(2*x-1)*(1-X**2), wrt)  # b1
         b_2 = sp.integrate(h*sp.Rational(1,4)*(2*x-1)*(X+1)*X, wrt)   # b2

         #print("Integrating with numerical python element vector in reference")
         #print("b0",b_0)
         #print("b1",b_1)
         #print("b2",b_2)

         self.b0_e = sp.lambdify((x_m,h),b_0,'numpy')
         self.b1_e = sp.lambdify((x_m,h),b_1,'numpy')
         self.b2_e = sp.lambdify((x_m,h),b_2,'numpy')

    def solve_linear(self, option = 0):
        """
        solves the assembled system using either pentapy or numpy linalg solve
        """
        if option == 0:
            self.c = pp.solve(self.A,self.b, is_flat = False)
        else:
            self.c = np.linalg.solve(self.A,self.b)

    def get_numerical_solution(self,x):
        self.x = x;
        self.assemble()
        self.solve_linear()
        self.evaluate_num(x)
        return self.u

    def exact(self,x):
        constant = (self.D-self.C-1/6)
        p3 = -1/3*x**3 + 1/2*x**2 + self.C*x
        return p3 + constant

    def L2_norm(self):
        """
        Calculate the L2-norm by approximating the integral as a sum
        """                          # Interval and variable
        L2 = np.sqrt(self.h*sum((self.exact(self.x)-self.u)**2))
        return L2, self.h

    def plot(self):
        """
        Plot the numerical and analytical solution
        """
        plt.plot(self.x,self.exact(self.x),'-',label = "analytical")
        plt.plot(self.x,self.u,'--', label = "numerical")
        plt.title("Finite Element Method")
        plt.xlabel("x",fontsize = 13)
        plt.ylabel("u",fontsize = 13)
        plt.legend(loc = "upper left",fontsize = 15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.show()
