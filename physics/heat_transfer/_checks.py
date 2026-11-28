"""열전달 커널 공통 입력 검사 (Role 3 내부용).

정책
    - 기하·상수 입력(면적, 길이, 수력직경, sigma)이 음수·비유한이면 ValueError. 설정 오류다.
    - 상태에서 계산한 계수(k, h_conv, Nu, emissivity)가 NaN이면 그대로 NaN을 전달하고,
      유한한 음수면 ValueError. 상관식/물성 구현 오류를 숨기지 않는다.
    - 온도가 0 이하이거나 비유한이면 그 원소만 NaN. Newton 선탐색이 비유한 잔차를 거부한다.
"""
import numpy as np


def as_float(x):
    return np.asarray(x, dtype=float)


def out(x):
    """0-d 배열은 NumPy 스칼라로."""
    return x[()] if isinstance(x, np.ndarray) and x.ndim == 0 else x


def require_geometry(name, value, *, allow_zero=False):
    value = as_float(value)
    bad = ~np.isfinite(value) | ((value < 0) if allow_zero else (value <= 0))
    if np.any(bad):
        rule = ">= 0" if allow_zero else "> 0"
        raise ValueError(f"{name} must be finite and {rule}; got {value[bad].ravel()[:3]!r}")
    return value


def require_coefficient(name, value, *, upper=None):
    value = as_float(value)
    finite = np.isfinite(value)
    if np.any(finite & (value < 0)):
        raise ValueError(f"{name} must be >= 0; got {value[finite & (value < 0)].ravel()[:3]!r}")
    if upper is not None and np.any(finite & (value > upper)):
        raise ValueError(f"{name} must be <= {upper}; got {value[finite & (value > upper)].ravel()[:3]!r}")
    return value


def bad_temperature(*temps):
    bad = None
    for T in temps:
        b = ~(np.isfinite(T) & (T > 0))
        bad = b if bad is None else (bad | b)
    return bad


def mask(value, bad):
    value = np.array(value, dtype=float, copy=True)
    if bad is not None:
        value[np.broadcast_to(bad, value.shape)] = np.nan
    return out(value)
