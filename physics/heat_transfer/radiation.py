"""회색체 복사 열률 (Role 3).

Q_ab = epsilon sigma A (T_a⁴ − T_b⁴)   [W], a → b 방향이 양수.
epsilon은 두 면 사이의 유효(교환) 방사율이다. 두 면의 방사율에서 유효값을 만드는 도구를 함께 둔다.
"""
import numpy as np

from ._checks import as_float, bad_temperature, mask, require_coefficient, require_geometry

# CODATA 2018 정확값. 변수표 sigma_SB 5.6704e-8은 이 값의 반올림.
STEFAN_BOLTZMANN = 5.670374419e-8   # W/(m²·K⁴)


def evaluate(T_a, T_b, epsilon, A, sigma):
    T_a, T_b = as_float(T_a), as_float(T_b)
    epsilon = require_coefficient("effective emissivity", epsilon, upper=1.0)
    A = require_geometry("area A", A, allow_zero=True)
    sigma = require_geometry("Stefan-Boltzmann constant sigma", sigma)
    with np.errstate(invalid="ignore", over="ignore"):
        q = epsilon * sigma * A * (T_a ** 4 - T_b ** 4)
    return mask(q, bad_temperature(T_a, T_b))


def effective_emissivity_gray_enclosure(eps_a, eps_b, area_ratio_a_to_b):
    """닫힌 두 회색 확산면 (a가 b에 완전히 둘러싸임: 동심 원통·구, 평행 평판은 면적비 1).

    epsilon_eff = 1 / (1/eps_a + (A_a/A_b) (1/eps_b − 1)),  Q_ab = epsilon_eff sigma A_a (T_a⁴ − T_b⁴)
    """
    eps_a = require_geometry("emissivity eps_a", eps_a)
    eps_b = require_geometry("emissivity eps_b", eps_b)
    ratio = require_geometry("area ratio A_a/A_b", area_ratio_a_to_b)
    if np.any(eps_a > 1) or np.any(eps_b > 1):
        raise ValueError("emissivities must be in (0, 1]")
    if np.any(ratio > 1):
        raise ValueError("A_a/A_b must be <= 1 for an inner surface a enclosed by b")
    return mask(1.0 / (1.0 / eps_a + ratio * (1.0 / eps_b - 1.0)), None)
