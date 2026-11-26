"""지수 1 미분대수방정식의 일치 초기조건과 정상상태 계산."""
import numpy as np
from .newton import newton_solve


def consistent_initialize(residual, t0, y_guess, differential_mask, *,
                          ydot_guess=None, options=None, residual_scale=None,
                          correction_scale=None, jacobian_sparsity=None):
    """미분 상태는 고정하고 대수 상태의 변화율은 0으로 두어 초기 잔차식 F=0을 푼다."""
    y0 = np.asarray(y_guess, dtype=float).copy()
    mask = np.asarray(differential_mask, dtype=bool)
    if y0.ndim != 1 or not y0.size or mask.shape != y0.shape:
        raise ValueError("differential_mask shape mismatch")
    yp0 = np.zeros_like(y0) if ydot_guess is None else np.asarray(ydot_guess, dtype=float).copy()
    if yp0.shape != y0.shape:
        raise ValueError("ydot_guess shape must match y_guess")
    if not np.all(np.isfinite(y0)) or not np.all(np.isfinite(yp0)):
        raise ValueError("initial values must be finite")
    yp0[~mask] = 0.0
    n_algebraic = np.count_nonzero(~mask)
    # 뉴턴 미지수의 앞부분에는 대수 상태를, 뒷부분에는 미분 상태의 변화율을 넣는다.
    z0 = np.concatenate((y0[~mask], yp0[mask]))

    def fun(z):
        yy, yyp = y0.copy(), yp0.copy()
        yy[~mask] = z[:n_algebraic]
        yyp[mask] = z[n_algebraic:]
        return residual(t0, yy, yyp)

    sol = newton_solve(fun, z0, options=options, residual_scale=residual_scale,
                       correction_scale=correction_scale,
                       jacobian_sparsity=jacobian_sparsity)
    y0[~mask] = sol.x[:n_algebraic]
    yp0[mask] = sol.x[n_algebraic:]
    return y0, yp0


def solve_steady_state(residual, t0, y_guess, *, options=None, residual_scale=None,
                       jacobian_sparsity=None):
    """주어진 시각에서 모든 변화율을 0으로 두고 평형 상태를 구해 (y, ydot)를 반환한다."""
    y_guess = np.asarray(y_guess, dtype=float)
    if y_guess.ndim != 1 or not y_guess.size or not np.all(np.isfinite(y_guess)):
        raise ValueError("y_guess must be a nonempty finite state vector")
    ydot = np.zeros_like(y_guess)
    solution = newton_solve(lambda y: residual(t0, y, ydot), y_guess,
                            options=options, residual_scale=residual_scale,
                            jacobian_sparsity=jacobian_sparsity)
    return solution.x, ydot
