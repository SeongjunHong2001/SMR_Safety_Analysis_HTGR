"""밀집 행렬로 주어진 선형 연립방정식 풀이."""
import numpy as np


class LinearSolveFailure(RuntimeError):
    """뉴턴 보정에 필요한 선형 연립방정식을 풀 수 없을 때 발생한다."""
    pass


def solve_linear(A, b):
    """A x = b를 직접 풀고, 분해가 실패하면 상위 뉴턴 반복에 실패를 알린다."""
    try:
        return np.linalg.solve(np.asarray(A, dtype=float), np.asarray(b, dtype=float))
    except np.linalg.LinAlgError as exc:
        raise LinearSolveFailure("Newton Jacobian is singular") from exc
