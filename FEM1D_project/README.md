# FYS4150-Project2: Eigenvalue Problems
Git for Project 2 in Numerical methods for partial differential equations (IN5270).

### Main overview
* The programmes in this repository aim at solving the questions about a 1D poisson equation posed in the project 2 description: [Project 2 - Finite Element Methods on a 1D PDE ](https://github.com/lasse-steinnes/IN5270/blob/master/FEM1D_project/Report/in5270_project2_description.pdf). The final report can be found at: [Solving a 1D Poisson Equation with the Galerkin Finite Element Method](https://github.com/lasse-steinnes/IN5270/blob/master/FEM1D_project/Report/SteinnesLSolving-a-1D-Poisson-Equation-with-the-Galerkin-Finite-Element-Method.pdf).

* The main challenge was to create a Finite Element Method solver using the Galerkin method, also known as the Projection method, and tailor the algorithm to solve a 1D poisson equation. Here, P2-elements were used. The problem boils down to linear algebra and solving a linear system.

* Another central task was to compute the L2-norm of the error, and find the converge rate.

* Example results can be found in the folder figures.

### Code: Link and description of programmes
- [main.py](https://github.com/lasse-steinnes/IN5270/blob/master/FEM1D_project/main.py) : Runs the methods provided by the class FEM_P2_solver.


- [FEM_P2.py](https://github.com/lasse-steinnes/IN5270/blob/master/FEM1D_project/FEM_P2.py) : Class methods provided are given in following order
  1. initialize : Set up the elements, boundary condtions, nodes and other essential parameters
  2. setup_element_matrix: Sets up element matrix for reference domain
  3. setup_element_matrixN: Sets up special element for last reference cell.
  4. Evaluate_num: Evaluates the numerical solution at mesh points
  5. Assemble: Assembles the global matrix
  6. find_b_ref: Find the right hand side element vector on the reference domain
  7. solve_linear: Solves the global linear system
  8. get_numerical_solution: Uses the other methods and returns u_num
  9. exact: Gives the exact solution as a function of x
  10. l2-norm: Calculates the l2 norm of e = u_e - u_num
  11. plot: Plots the numerical and analytical solution in the same figure

How to run the programmes to reproduce the results discussed in the article: The menu gives input options for C, D, number of elements to use (Ne) and if you want to plot or calculate the convergence rates. 

To calculate convergence rates, a number of experiments must be provided. For each experiment Ne = 4*Ne and h = 1/N_e. To reproduce the convergence rates described in the paper, use [C:1, D:8, 5 elements 7 experiments] and  [C:0, D:0.5, 5 elements 7 experiments].

### Links and packages
  - Documentation for numpy [here](https://numpy.org/doc/)
  - Documentation for matplotlib [here](https://matplotlib.org/)
  - Documentation for sympy [here](https://docs.sympy.org/latest/index.html)
  - Documentation for pentapy [here](https://geostat-framework.readthedocs.io/projects/pentapy/en/stable/)
  

###  FEEDBACK and necessary fixes:
You have written a very good report, showing that you have grasped the mathematical and numerical details related to the finite element method when solving a 1D Poisson equation. The handling of the Neumann and Dirichlet conditions is also correct. (My only comment in this regard concerns page 5, right above formula (35), "the extended matrix" is an unusual notation.) The use of the standardized reference element is excellent. Your implementation seems to be (mostly) correct because you have achieved convergence with you numerical solutions.

However, your achieved convergence rate (3.5) is better than expected, thus a "mystery". The possible reasons are listed below:
1. You have used a fixed array of sampling points, specifically,
x = np.linspace(0,1, 10000)
which is independent of the number of elements actually used. This MAY lead to insufficient accuracy of the numerical integration (essentially used in the L2_norm function) when the number of cells is really large. By the way, there is no need to increase the number of cells by a factor of 4 every time. (An increase factor of 2 is standard.)
2. The L2_norm function is incorrect in the sense that you multiply with "self.h", which is mathematically correct (if you had sampled a fixed number of points per element), but it does NOT match with your actual fixed sampling rate (always 10000 points over the entire domain no matter the actual number of elements used).
3. The function "evaluate_num" may require a closer look. I must confess that I find the current implementation somewhat unusual, but cannot quite pinpoint any bug (if at all).

Your project is approved, although the convergence rate achieved is too good (better than what the theory expects).
- Xing Cai
  
