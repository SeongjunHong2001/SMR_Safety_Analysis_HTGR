"""잔차식의 크기를 조정하고 후퇴 선 탐색으로 보정량을 줄이는 뉴턴법."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from .jacobian import finite_difference_jacobian
from .error_control import rms_norm
from .linear_solver import LinearSolveFailure, solve_linear

class NewtonFailure(RuntimeError):
    """뉴턴 반복에서 유효한 보정량이나 수렴한 해를 얻지 못할 때 발생한다."""
    pass

@dataclass
class NewtonOptions:
    """최대 반복 횟수, 잔차 허용오차, 선 탐색에서 허용할 최소 감쇠 계수."""
    max_iter: int = 20
    rtol: float = 1e-9
    atol: float = 1e-10
    min_damping: float = 1/1024

@dataclass
class NewtonResult:
    """근의 추정값, 수행한 반복 횟수, 척도로 나눈 최종 잔차의 무한대 노름."""
    x: np.ndarray
    iterations: int
    residual_norm: float

def newton_solve(fun, x0, *, options=None, residual_scale=None,
                 jacobian=None, correction_scale=None, fd_scales=1.0,
                 jacobian_sparsity=None) -> NewtonResult:
    """잔차 벡터의 크기와 미지수 수가 같은 비선형 연립방정식 fun(x) = 0을 푼다."""
    options = options or NewtonOptions()
    x = np.asarray(x0, dtype=float).copy()
    if x.ndim != 1 or not x.size or not np.all(np.isfinite(x)):
        raise ValueError("Newton initial state must be a nonempty finite 1-D vector")
    r = np.asarray(fun(x), dtype=float)
    if r.shape != x.shape:
        raise ValueError("Newton residual must be a 1-D vector matching the state")
    scale = np.ones_like(r) if residual_scale is None else np.asarray(residual_scale, dtype=float)
    if scale.shape != r.shape or np.any(scale <= 0) or not np.all(np.isfinite(scale)):
        raise ValueError("residual_scale must be positive, finite and match residual dimension")
    if options.max_iter < 0 or not 0 < options.min_damping <= 1:
        raise ValueError("Newton requires max_iter >= 0 and 0 < min_damping <= 1")
    if not np.isfinite(options.atol) or not np.isfinite(options.rtol):
        raise ValueError("Newton tolerances must be finite")
    if options.atol < 0 or options.rtol < 0:
        raise ValueError("Newton tolerances must be nonnegative")
    norm = lambda rr: float(np.linalg.norm(rr/scale, ord=np.inf)) if rr.size else 0.0
    err = norm(r)
    if not np.isfinite(err):
        raise NewtonFailure("Newton initial residual is non-finite")
    target = options.atol + options.rtol * max(1.0, err)
    for k in range(options.max_iter + 1):
        if err == 0 or (correction_scale is None and err <= target):
            return NewtonResult(x=x, iterations=k, residual_norm=err)
        if k == options.max_iter and (correction_scale is None or k == 0):
            break
        jac = (finite_difference_jacobian(fun, x, step_scales=fd_scales,
                                          sparsity=jacobian_sparsity)
               if jacobian is None else jacobian(x))
        jac = jac.astype(float) if sparse.issparse(jac) else np.asarray(jac, dtype=float)
        if jac.shape != (r.size, x.size):
            raise ValueError("Newton Jacobian shape must match the system")
        entries = jac.data if sparse.issparse(jac) else jac
        if not np.all(np.isfinite(entries)):
            raise NewtonFailure("Newton Jacobian is non-finite")
        try:
            # 단위가 다른 잔차식을 비교할 수 있도록 야코비안의 각 행과 우변을 같은 척도로 나눈다.
            scaled_jac = (jac.multiply(1.0 / scale[:, None]) if sparse.issparse(jac)
                          else jac / scale[:, None])
            dx = solve_linear(scaled_jac, -r / scale)
        except LinearSolveFailure as exc:
            raise NewtonFailure(str(exc)) from exc
        # 보정 척도가 주어지면 잔차와 함께 보정량의 가중 RMS도 수렴 판정에 사용한다.
        small_correction = correction_scale is not None and rms_norm(dx / correction_scale) <= 0.03
        if k > 0 and small_correction and err <= target:
            return NewtonResult(x=x, iterations=k, residual_norm=err)
        if k == options.max_iter:
            break
        # 전체 보정량부터 시험하고 잔차가 감소할 때까지 보정량을 절반씩 줄인다.
        damping = 1.0
        while damping >= options.min_damping:
            candidate = x + damping * dx
            try:
                r_try = np.asarray(fun(candidate), dtype=float)
                if r_try.shape != r.shape:
                    raise ValueError("Newton residual shape changed during line search")
                trial_error = norm(r_try)
            except FloatingPointError:
                trial_error = np.inf
            if np.isfinite(trial_error) and trial_error < err:
                x, r, err = candidate, r_try, trial_error
                if small_correction and damping == 1.0 and err <= target:
                    return NewtonResult(x=x, iterations=k + 1, residual_norm=err)
                break
            damping *= 0.5
        else:
            if small_correction and err <= target:
                return NewtonResult(x=x, iterations=k, residual_norm=err)
            raise NewtonFailure(f"Newton line search failed after iteration {k+1}; norm={err:.3e}")
    raise NewtonFailure(f"Newton did not converge in {options.max_iter} iterations; norm={err:.3e}")
