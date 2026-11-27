"""1차원 평판 전도 열률 (Role 3).

Q_ab = k A (T_a − T_b) / L   [W], a → b 방향이 양수.
k [W/(m·K)], A [m²], L [m], T [K]. 입력은 NumPy broadcasting.
원통 벽 등 다른 기하는 호출자가 등가 A/L을 넘긴다 (예: 원통 A/L = 2πH / ln(r_o/r_i)).
"""
import numpy as np

from ._checks import as_float, bad_temperature, mask, require_coefficient, require_geometry


def evaluate(T_a, T_b, k, A, L):
    T_a, T_b = as_float(T_a), as_float(T_b)
    k = require_coefficient("conductivity k", k)
    A = require_geometry("area A", A, allow_zero=True)
    L = require_geometry("length L", L)
    with np.errstate(invalid="ignore"):
        q = k * A * (T_a - T_b) / L
    return mask(q, bad_temperature(T_a, T_b))
