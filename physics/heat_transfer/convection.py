"""대류 열률과 Nusselt → 열전달계수 변환 (Role 3).

Q_ab = h_conv A (T_a − T_b)   [W], a → b 방향이 양수.
h_conv = Nu k / D_h           [W/(m²·K)]
"""
import numpy as np

from ._checks import as_float, bad_temperature, mask, require_coefficient, require_geometry


def evaluate(T_a, T_b, h_conv, A):
    T_a, T_b = as_float(T_a), as_float(T_b)
    h_conv = require_coefficient("h_conv", h_conv)
    A = require_geometry("area A", A, allow_zero=True)
    with np.errstate(invalid="ignore"):
        q = h_conv * A * (T_a - T_b)
    return mask(q, bad_temperature(T_a, T_b))


def coefficient_from_nusselt(Nu, k, D_h):
    Nu = require_coefficient("Nusselt number", Nu)
    k = require_coefficient("conductivity k", k)
    D_h = require_geometry("hydraulic diameter D_h", D_h)
    with np.errstate(invalid="ignore"):
        return mask(Nu * k / D_h, None)
