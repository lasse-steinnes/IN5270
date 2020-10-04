# IN5270 Project 1
Git for project 1 in IN5270 - Numerical Methods for Partial Differential Equations: Solving Partial Differential Equations using Finite Difference Methods.

### Main overview
* The programs in this repository aim at solving the questions about a general 2D wave equation posed in the project 1 description: [Project 1 - Finite difference simulation of 2D waves ](https://github.com/lasse-steinnes/IN5270/blob/master/wave_project/Report/in5270_oblig_project1.pdf). The final report can be found at: [Solving a 2D Wave Equation using Finite Difference Methods](https://github.com/lasse-steinnes/IN5270/blob/master/wave_project/Report/SteinnesL-Solving-a-2D-wave-eq-using-finite-difference-methods.pdf).

* The main challenge was to create a numerical scheme using finite difference method to solve a general 2D-wave equation, then apply it to differential equations known from Physics.

* Another central task was to check crucial functions of the algorithm, using unit tests for verification. In addition, convergence test for a damped and undamped 2D wave is performed.

* Finally, to explore numerical artefacts, a 3D animation time series is made for a wave propagating over subsea hills of 3 different shapes.

### Code - link and description of programmes
- [main.py](https://github.com/lasse-steinnes/FYS4150-Project2/blob/master/code-and-results/main.cpp) : Runs the other programmes and provide user options through terminal.
- [wave2D.py](https://github.com/lasse-steinnes/IN5270/blob/master/wave_project/wave2D.py) : A class for the 2D numerical scheme to solve a 2D wave equation. The methods given in the class is
  1. initialize : Sets up class parameters
  2. set_inital : Set initial condition
  3. var_coeff : Uses a numerical approximation to the composite derivative (including wave velocity g)
  4. set_boundary : Set the boundary conditions (Neumann)
  5.  advance: Advances the numerical scheme in space
  6. advance_b0: Advances the numerical scheme in space for b = 0 (no damping).
  7. solve : Advances numerical solution for all time steps
  8. test_convergence_rate: Find convergence rate for given exact solution
  9. intialize convergence: sets up initial conditions for each experiment used in test_convergence_rate
- [vizualize.py](https://github.com/lasse-steinnes/IN5270/blob/master/wave_project/visualize.py) : Visualizes the numerical solution in 3D for a given time step.
- [seabed.py](https://github.com/lasse-steinnes/IN5270/blob/master/wave_project/seabed.py) : Explores the different subsea surfaces and the wave velocity they produce.
