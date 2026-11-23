import numpy as np
import pytest
from solvers.newton import newton_solve

def test_nonlinear_root():
    sol = newton_solve(lambda x: np.array([x[0]**2-4]), np.array([3.0]))
    assert abs(sol.x[0]-2) < 1e-8


def test_newton_rejects_nonvector_or_changing_residual():
    with pytest.raises(ValueError, match="1-D vector"):
        newton_solve(lambda x: np.array([[x[0] - 1.0]]), [0.0])

    def changing(x):
        return np.array([x[0] - 1.0]) if x[0] == 0.0 else np.array([x[0] - 1.0, 0.0])

    with pytest.raises(ValueError, match="shape changed"):
        newton_solve(changing, [0.0])
