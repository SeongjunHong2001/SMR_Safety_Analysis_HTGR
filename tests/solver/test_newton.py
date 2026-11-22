import numpy as np
from solvers.newton import newton_solve

def test_nonlinear_root():
    sol = newton_solve(lambda x: np.array([x[0]**2-4]), np.array([3.0]))
    assert abs(sol.x[0]-2) < 1e-8
