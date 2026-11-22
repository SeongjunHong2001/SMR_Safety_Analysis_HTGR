import numpy as np
import pytest
from solvers.bdf import bdf_derivative, bdf_predictor

@pytest.mark.parametrize('order', range(1, 6))
def test_bdf_derivative_linear(order):
    # y(t)=1+2t이면 모든 BDF order가 정확히 2를 반환
    dt=0.2
    t_new = 1.0
    prev = [np.array([1.0+2*(t_new-dt*j)]) for j in range(1,order+1)]
    got = bdf_derivative(np.array([1.0+2*t_new]), prev, dt, order)
    assert np.allclose(got, [2.0], atol=1e-12)

def test_bdf_history_validation():
    with pytest.raises(ValueError):
        bdf_derivative(np.array([1.0]), [np.array([0.0])], .1, 2)


@pytest.mark.parametrize("order", range(1, 6))
def test_variable_step_predictor_is_exact_for_polynomial_of_same_order(order):
    t_new = 1.0
    times = np.array([0.9, 0.75, 0.56, 0.5, 0.25, 0.1])[:order + 1]
    states = (2.0 + times ** order)[:, None]
    y, yp, cj, error_factor = bdf_predictor(t_new, times, states, order)
    assert y[0] == pytest.approx(3.0)
    assert yp[0] == pytest.approx(order)
    assert cj == pytest.approx(sum(1 / j for j in range(1, order + 1)) / (t_new - times[0]))
    assert error_factor > 0


@pytest.mark.parametrize("order", range(1, 6))
def test_fixed_leading_corrector_matches_classical_bdf_on_uniform_grid(order):
    times = np.arange(0, -order - 1, -1, dtype=float)
    states = (times ** (order + 1))[:, None]
    y_new = np.array([1.0])
    prediction, slope, cj, error_factor = bdf_predictor(1.0, times, states, order)
    assert np.allclose(slope + cj * (y_new - prediction), bdf_derivative(y_new, states, 1.0, order))
    assert error_factor == pytest.approx(1.0 / (order + 1))
