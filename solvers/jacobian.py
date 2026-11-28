"""잔차 벡터를 유한차분해 야코비안 행렬을 구한다."""
import numpy as np
from scipy import sparse

def finite_difference_jacobian(fun, x, *, rel_step=None, step_scales=1.0,
                               sparsity=None):
    """Compute a dense Jacobian, or a colored sparse Jacobian when a pattern is given.

    ``sparsity[i, j]`` must be true for every possible nonzero derivative.
    Columns with disjoint residual rows are perturbed in one function call.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not x.size or not np.all(np.isfinite(x)):
        raise ValueError("Jacobian state must be a nonempty finite 1-D vector")
    f0 = np.asarray(fun(x), dtype=float)
    if f0.ndim != 1:
        raise ValueError("Jacobian residual must be a 1-D vector")
    m, n = f0.size, x.size
    step_factor = np.sqrt(np.finfo(float).eps) if rel_step is None else float(rel_step)
    steps = step_factor * np.maximum(step_scales, np.abs(x))
    if np.shape(steps) != x.shape or np.any(steps <= 0) or not np.all(np.isfinite(steps)):
        raise ValueError("Finite-difference steps must be positive, finite and match the state")
    if sparsity is not None:
        return _sparse_jacobian(fun, x, f0, steps, sparsity)
    J = np.empty((m, n), dtype=float)
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


def _sparse_jacobian(fun, x, f0, steps, sparsity):
    # copy=True: 일부 SciPy(1.15.x)는 dtype 변환에서도 입력 행렬의 인덱스 배열을 공유한다.
    # 아래 sum_duplicates/eliminate_zeros가 호출자의 행렬을 바꾸지 않도록 사본을 만든다.
    pattern = sparse.csc_matrix(sparsity, dtype=bool, copy=True)
    if pattern.shape != (f0.size, x.size):
        raise ValueError("Jacobian sparsity must match residual and state dimensions")
    pattern.sum_duplicates()
    pattern.eliminate_zeros()
    row_colors = [set() for _ in range(f0.size)]
    groups = {}
    for column in range(x.size):
        rows = pattern.indices[pattern.indptr[column]:pattern.indptr[column + 1]]
        if not rows.size:
            continue
        unavailable = set().union(*(row_colors[row] for row in rows))
        color = 0
        while color in unavailable:
            color += 1
        groups.setdefault(color, []).append(column)
        for row in rows:
            row_colors[row].add(color)

    def evaluate(trial):
        try:
            values = np.asarray(fun(trial), dtype=float)
        except FloatingPointError:
            return None
        if values.shape != f0.shape:
            raise ValueError("Jacobian residual shape changed during finite differences")
        return values if np.all(np.isfinite(values)) else None

    data, row_indices, column_indices = [], [], []
    for columns in groups.values():
        trial = x.copy()
        trial[columns] += steps[columns]
        if np.any(trial[columns] == x[columns]):
            raise FloatingPointError("Finite-difference perturbation cannot change the state")
        values = evaluate(trial)
        if values is None:
            # A grouped perturbation can cross a model boundary. Retry each
            # column separately and use the negative direction when needed.
            candidates = []
            for column in columns:
                single = x.copy()
                single[column] += steps[column]
                result = evaluate(single)
                if result is None:
                    single[column] = x[column] - steps[column]
                    result = evaluate(single)
                if result is None or single[column] == x[column]:
                    raise FloatingPointError(
                        f"Finite-difference perturbation failed for column {column}"
                    )
                candidates.append((column, single[column] - x[column], result))
        else:
            candidates = [(column, trial[column] - x[column], values)
                          for column in columns]
        for column, delta, result in candidates:
            rows = pattern.indices[pattern.indptr[column]:pattern.indptr[column + 1]]
            row_indices.extend(rows)
            column_indices.extend([column] * len(rows))
            data.extend((result[rows] - f0[rows]) / delta)
    return sparse.csc_matrix((data, (row_indices, column_indices)), shape=pattern.shape)
