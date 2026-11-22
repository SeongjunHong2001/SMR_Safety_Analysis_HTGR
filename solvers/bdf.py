"""등간격 후진차분공식(BDF)과 가변 간격 선두 계수 고정형 BDF의 예측값 계산."""
import numpy as np

# 등간격 미분 공식의 계수를 새 상태부터 과거 상태 순서로 저장한다.
_COEFF = {
    1: np.array([1.0, -1.0]),
    2: np.array([1.5, -2.0, 0.5]),
    3: np.array([11/6, -3.0, 1.5, -1/3]),
    4: np.array([25/12, -4.0, 3.0, -4/3, 1/4]),
    5: np.array([137/60, -5.0, 5.0, -10/3, 5/4, -1/5]),
}


def coefficients(order: int) -> np.ndarray:
    """지정한 차수의 등간격 BDF 계수를 복사해 반환한다."""
    if order not in _COEFF:
        raise ValueError("BDF order must be between 1 and 5")
    return _COEFF[order].copy()


def bdf_derivative(y_new: np.ndarray, history, dt: float, order: int) -> np.ndarray:
    """최신 상태부터 담긴 이력(history[0] = y_n)으로 새 시각의 등간격 BDF 미분값을 구한다."""
    if dt <= 0:
        raise ValueError("dt must be positive")
    coeff = coefficients(order)
    if len(history) < order:
        raise ValueError(f"BDF{order} needs {order} previous states")
    y_new = np.asarray(y_new, dtype=float)
    derivative = np.zeros_like(y_new)
    # 계수의 합이 0이므로 각 과거 상태와 새 상태의 차이만으로 미분값을 계산할 수 있다.
    for weight, previous in zip(coeff[1:], history):
        derivative += weight * (np.asarray(previous, dtype=float) - y_new)
    return derivative / dt


def divided_differences(nodes, states):
    """최신 시각부터 정렬한 무차원 절점에서 뉴턴 보간의 분할차분 계수를 구한다."""
    coeff = np.array(states, dtype=float, copy=True)
    for degree in range(1, len(nodes)):
        distances = nodes[degree:] - nodes[:-degree]
        coeff[degree:] = (coeff[degree:] - coeff[degree - 1:-1]) / distances[:, None]
    return coeff


def bdf_predictor(t_new, past_times, past_states, order):
    """상태 예측값, 예측 시간미분, 보정 계수 cj, 국소 절단 오차 추정용 계수를 반환한다."""
    step = t_new - past_times[0]
    # 새 시각을 0, 직전 시각을 -1로 두면 가변 간격 이력도 같은 좌표에서 보간할 수 있다.
    nodes = (np.asarray(past_times[:order + 1]) - t_new) / step
    coeff = divided_differences(nodes, past_states[:order + 1])
    value = coeff[-1]
    slope = np.zeros_like(value)
    for node, term in zip(nodes[-2::-1], coeff[-2::-1]):
        slope = value - node * slope
        value = term - node * value

    # 상태 보정량에 cj를 곱해 예측 시간미분에 더한다. 등간격에서는 통상적인 BDF 미분식과 같다.
    harmonic = _COEFF[order][0]
    alpha = -1.0 / nodes
    error_factor = max(abs(np.sum(alpha) - harmonic), alpha[-1])
    return value, slope / step, harmonic / step, error_factor
