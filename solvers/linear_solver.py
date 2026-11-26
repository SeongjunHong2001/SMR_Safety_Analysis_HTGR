"""밀집 행렬로 주어진 선형 연립방정식 풀이."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


class LinearSolveFailure(RuntimeError):
    """뉴턴 보정에 필요한 선형 연립방정식을 풀 수 없을 때 발생한다."""
    pass


def solve_linear(A, b):
    """A x = b를 직접 풀고, 분해가 실패하면 상위 뉴턴 반복에 실패를 알린다."""
    try:
        rhs = np.asarray(b, dtype=float)
        if sparse.issparse(A):
            solution = splu(A.tocsc()).solve(rhs)
        else:
            solution = np.linalg.solve(np.asarray(A, dtype=float), rhs)
    except (np.linalg.LinAlgError, RuntimeError) as exc:
        raise LinearSolveFailure("Newton Jacobian is singular") from exc
    if not np.all(np.isfinite(solution)):
        raise LinearSolveFailure("Newton correction is non-finite")
    return solution
