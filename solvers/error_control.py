"""가중 RMS 오차와 BDF 스텝 크기 제어."""
import numpy as np


def error_weights(state, *, atol, rtol, previous=None):
    """상태별 절대·상대 허용오차로 오차를 나눌 척도를 구한다."""
    magnitude = np.abs(state)
    # 두 시점 중 큰 절댓값을 써서 상태가 0에 가까워질 때도 직전 크기를 반영한다.
    if previous is not None:
        magnitude = np.maximum(magnitude, np.abs(previous))
    weights = np.asarray(atol, dtype=float) + rtol * magnitude
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ValueError("error weights must be positive and finite")
    return weights


def weighted_rms_norm(error, state, *, atol, rtol):
    """각 성분의 오차를 허용오차로 나눈 뒤 제곱 평균의 제곱근(RMS)을 구한다."""
    return rms_norm(np.asarray(error, dtype=float) / error_weights(state, atol=atol, rtol=rtol))


def rms_norm(values):
    # 가장 큰 절댓값으로 먼저 나누어 제곱 과정에서 큰 값이 넘치는 것을 피한다.
    magnitude = np.max(np.abs(values))
    if magnitude == 0 or not np.isfinite(magnitude):
        return float(magnitude)
    return float(magnitude * np.sqrt(np.mean((values / magnitude) ** 2)))


def step_factor(error, order):
    """국소 오차가 간격의 (차수 + 1)제곱에 비례한다고 보고 안전 계수 0.9를 적용한 간격 배율을 구한다."""
    if not np.isfinite(error):
        return 0.0
    return 0.9 * max(error, 1e-12) ** (-1.0 / (order + 1))
