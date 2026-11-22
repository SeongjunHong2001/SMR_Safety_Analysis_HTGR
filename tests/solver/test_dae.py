import numpy as np
import pytest

from solvers import DASSL, DASSLPrototype, IntegrationFailure
from solvers.initialization import consistent_initialize

def test_linear_index1_dae():
    # x'=-2x, z=3x; exact x=e^-2t, z=3e^-2t
    def residual(t,y,ydot):
        x,z=y
        return np.array([ydot[0]+2*x, z-3*x])
    y0, yp0 = consistent_initialize(residual, 0.0,
                                     np.array([1.0, 0.0]),
                                     np.array([True,False]))
    assert np.allclose(y0, [1.0,3.0], atol=1e-8)
    assert abs(yp0[0]+2) < 1e-8
    res=DASSLPrototype(dt=0.025,max_order=2).solve(
        residual,t_span=(0,1),y0=y0,ydot0=yp0)
    expected=np.exp(-2)
    assert abs(res.y[-1,0]-expected) < 0.004
    assert abs(res.y[-1,1]-3*res.y[-1,0]) < 1e-8
    assert np.all(np.diff(res.t)>0)


def test_adaptive_stiff_dae_with_exact_solution_and_analytic_jacobian():
    stiffness = 1e5
    jacobian_calls = []

    def residual(t, y, yp):
        x, z = y
        return np.array([yp[0] + stiffness * (x - z) + np.sin(t), z - np.cos(t)])

    def jacobian(t, y, yp, cj):
        jacobian_calls.append(cj)
        return np.array([[stiffness + cj, -stiffness], [0.0, 1.0]])

    def solve(rtol):
        return DASSL(dt=0.3, rtol=rtol, atol=[1e-8, 1e-9], max_step=0.4).solve(
            residual, t_span=(0.0, 5.037), y0=[1.0, 0.0],
            differential_mask=[True, False], jacobian=jacobian,
        )

    loose = solve(1e-3)
    tight = solve(1e-6)
    steps = np.diff(tight.t)
    assert tight.t[-1] == 5.037
    assert np.all(steps > 0)
    assert steps.max() <= 0.4 + 1e-14
    assert not np.allclose(steps, steps[0])
    assert set(tight.orders) == {1, 2, 3, 4, 5}
    assert tight.error_test_failures > 0
    assert np.all(tight.error_norms <= 1.0)
    assert jacobian_calls
    assert np.allclose(tight.y[:, 1], np.cos(tight.t), atol=1e-9)

    def max_error(result):
        return np.max(np.abs(result.y[:, 0] - np.cos(result.t)))

    assert max_error(tight) < min(5e-6, max_error(loose) / 5)
    for t, y, yp in zip(tight.t, tight.y, tight.ydot):
        assert np.max(np.abs(residual(t, y, yp))) < 2e-6


def test_retries_when_newton_predictor_leaves_residual_domain():
    def residual(t, y, yp):
        if abs(y[0] - np.cos(t)) > 0.15:
            raise FloatingPointError("trial state outside domain")
        return yp + 1e4 * (y - np.cos(t)) + np.sin(t)

    solver = DASSL(dt=0.8)
    result = solver.solve(residual, t_span=(0.0, 2.0), y0=[1.0])
    assert result.newton_failures > 0
    assert result.t[-1] == 2.0
    assert np.allclose(result.y[:, 0], np.cos(result.t), atol=5e-5)
    with pytest.raises(IntegrationFailure, match="trial state outside domain"):
        DASSL(dt=0.8, min_step=0.8).solve(residual, t_span=(0.0, 2.0), y0=[1.0])


def test_small_step_with_nonlinear_derivative_uses_scaled_difference():
    result = DASSL(dt=1e-12, atol=1e-15).solve(
        lambda t, y, yp: np.exp(yp) - 2.0,
        t_span=(0.0, 1e-12), y0=[0.0], ydot0=[np.log(2)],
    )
    assert result.ydot[-1, 0] == pytest.approx(np.log(2), abs=1e-9)
    assert result.y[-1, 0] == pytest.approx(1e-12 * np.log(2), abs=1e-20)


def test_backward_integration_with_tiny_residual_scale():
    result = DASSL(rtol=1e-7).solve(
        lambda t, y, yp: 1e-20 * (yp - y), t_span=(0.0, -3.0), y0=[1.0],
    )
    assert np.all(np.diff(result.t) < 0)
    assert result.t[-1] == -3.0
    assert result.ydot[0, 0] == pytest.approx(1.0)
    assert np.allclose(result.y[:, 0], np.exp(result.t), atol=5e-7, rtol=0.0)
