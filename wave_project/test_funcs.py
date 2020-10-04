from wave2D import Wave_2D_solver
import numpy as np
import pytest

class TestClass():
    def test_constant_solution(self):
        global b, U                  # make constant variables globally accessable
        U = 1.0;

        global f, g, V, I, t, x, y   # make functions globally accessable
        f = lambda x,y,t: 0
        I = lambda x,y: U
        V = lambda x,y: 0
        g = lambda x,y: 1

        dt = 0.1; T = 2; Nx = 100; Ny = 100;
        Lx = 4; Ly = 4;  b = 0

        u_e = np.full((Nx+1,Ny+1),U)
        # solve and compare with analytical
        solver = Wave_2D_solver(g,Lx,Ly,Nx,Ny,T,dt,'none')
        u, t = solver.solve(I,b,f,V)

        tol = 1E-12
        residual = np.subtract(u_e,u)
        print('\nMaximum error in computing a constant solution {:.3e}\n'.format(np.amax(residual)))
        assert pytest.approx(u_e, abs = tol) == u
