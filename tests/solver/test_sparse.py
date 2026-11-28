import numpy as np
import pytest
from scipy import sparse

from solvers import DASSL
from solvers.jacobian import finite_difference_jacobian
from solvers.linear_solver import LinearSolveFailure, solve_linear
from solvers.newton import newton_solve


def test_sparse_coloring_matches_matrix_with_fewer_residual_evaluations():
    size = 12
    matrix = sparse.diags(
        (np.ones(size - 1), 4 * np.ones(size), -2 * np.ones(size - 1)),
        offsets=(-1, 0, 1), format="csc",
    )
    calls = []

    def residual(x):
        calls.append(1)
        return matrix @ x + 1.0

    jacobian = finite_difference_jacobian(
        residual, np.linspace(1, 2, size), sparsity=matrix,
    )
    assert sparse.issparse(jacobian)
    assert np.allclose(jacobian.toarray(), matrix.toarray(), atol=1e-7)
    assert len(calls) <= 4  # baseline plus three nonoverlapping column groups


def test_sparse_coloring_retries_individual_columns_at_a_domain_boundary():
    def residual(x):
        if x[0] > 0:
            raise FloatingPointError("outside domain")
        return np.array([x[0], 2 * x[1]])

    jacobian = finite_difference_jacobian(
        residual, np.zeros(2), sparsity=sparse.eye(2, format="csc"),
    )
    assert np.allclose(jacobian.toarray(), np.diag([1.0, 2.0]))


def test_sparse_linear_and_analytic_newton_paths():
    matrix = sparse.diags([2.0, 3.0, 4.0], format="csc")
    rhs = np.array([2.0, 6.0, 12.0])
    assert np.allclose(solve_linear(matrix, rhs), [1.0, 2.0, 3.0])
    result = newton_solve(lambda x: matrix @ x - rhs, np.zeros(3),
                          jacobian=lambda x: matrix)
    assert np.allclose(result.x, [1.0, 2.0, 3.0])
    with pytest.raises(LinearSolveFailure, match="singular"):
        solve_linear(sparse.csc_matrix([[1.0, 1.0], [2.0, 2.0]]), [1.0, 2.0])


def test_dassl_sparse_pattern_matches_dense_and_calls_preflight():
    size = 8
    residual = lambda t, y, ydot: ydot + y
    options = dict(dt=0.1, max_order=2, rtol=1e-4, atol=1e-6)
    kwargs = dict(t_span=(0.0, 0.5), y0=np.ones(size),
                  differential_mask=np.ones(size, dtype=bool))
    dense = DASSL(**options).solve(residual, **kwargs)
    checks = []

    def preflight(t, y, ydot, *, cj):
        checks.append((t, cj, np.max(np.abs(residual(t, y, ydot)))))

    colored = DASSL(**options).solve(
        residual, jacobian_sparsity=sparse.eye(size, format="csc"),
        preflight=preflight, **kwargs,
    )
    assert np.array_equal(colored.t, dense.t)
    assert np.allclose(colored.y, dense.y, atol=1e-7)
    assert len(checks) == 1
    assert checks[0][0] == 0.0 and checks[0][1] == pytest.approx(10.0)
    assert checks[0][2] < 1e-8


def test_sparse_pattern_input_is_not_modified():
    # 패턴 행렬에 명시적 0이 있어도 호출자의 행렬(값·인덱스)은 그대로여야 한다 (SciPy 1.15.x 인덱스 공유 회귀)
    pattern = sparse.csc_matrix(np.array([[1.0, 0.0], [0.0, 1.0]]))
    pattern.data[0] = 0.0                      # 저장된 명시적 0
    before = (pattern.data.copy(), pattern.indices.copy(), pattern.indptr.copy(), pattern.nnz)
    finite_difference_jacobian(lambda x: 2.0 * x, np.ones(2), sparsity=pattern)
    assert pattern.nnz == before[3]
    assert np.array_equal(pattern.data, before[0]) and np.array_equal(pattern.indices, before[1])
    assert np.array_equal(pattern.indptr, before[2])
