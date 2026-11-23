"""잔차 벡터를 유한차분해 야코비안 행렬을 구한다."""
import numpy as np

def finite_difference_jacobian(fun, x, *, rel_step=None, step_scales=1.0):
    """입력 성분을 하나씩 바꿔 잔차의 변화율을 야코비안의 각 열에 저장한다."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not x.size or not np.all(np.isfinite(x)):
        raise ValueError("Jacobian state must be a nonempty finite 1-D vector")
    f0 = np.asarray(fun(x), dtype=float)
    if f0.ndim != 1:
        raise ValueError("Jacobian residual must be a 1-D vector")
    m, n = f0.size, x.size
    J = np.empty((m, n), dtype=float)
    step_factor = np.sqrt(np.finfo(float).eps) if rel_step is None else float(rel_step)
    steps = step_factor * np.maximum(step_scales, np.abs(x))
    if np.shape(steps) != x.shape or np.any(steps <= 0) or not np.all(np.isfinite(steps)):
        raise ValueError("Finite-difference steps must be positive, finite and match the state")
    for j, h in enumerate(steps):
        trial = x.copy()
        trial[j] += h
        try:
            fj = np.asarray(fun(trial), dtype=float)
        except FloatingPointError:
            fj = np.full_like(f0, np.nan)
        # 양의 방향에서 잔차를 계산할 수 없으면 원래 값의 반대쪽에서 다시 평가한다.
        if not np.all(np.isfinite(fj)):
            trial[j] = x[j] - h
            fj = np.asarray(fun(trial), dtype=float)
        if fj.shape != f0.shape:
            raise ValueError("Jacobian residual shape changed during finite differences")
        if trial[j] == x[j]:
            raise FloatingPointError("Finite-difference perturbation cannot change the state")
        # 반올림과 차분 방향을 반영하도록 실제로 적용된 변화량을 분모에 쓴다.
        J[:, j] = (fj - f0) / (trial[j] - x[j])
    return J
