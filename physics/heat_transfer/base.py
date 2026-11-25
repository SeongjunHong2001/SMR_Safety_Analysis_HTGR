"""Stable scalar/array heat-transfer signatures; physics implementation is pending.

All heat rates are positive from side a to side b (W). Temperatures are K.
Input arrays follow NumPy broadcasting, and providers must validate their
physical application range before returning a coefficient.
"""


def conduction(T_a, T_b, k, A, L):
    """Conductive W rate; k W/(m K), A m², L m."""
    from .conduction import evaluate
    return evaluate(T_a, T_b, k, A, L)


def convection(T_a, T_b, h_conv, A):
    """Convective W rate; h_conv W/(m² K)."""
    from .convection import evaluate
    return evaluate(T_a, T_b, h_conv, A)


def radiation(T_a, T_b, epsilon, A, sigma):
    """Radiative W rate; absolute temperatures K, sigma W/(m² K⁴)."""
    from .radiation import evaluate
    return evaluate(T_a, T_b, epsilon, A, sigma)


def heat_transfer_coefficient(Nu, k, D_h):
    """Return h_conv in W/(m² K) from Nu, k W/(m K), D_h m."""
    from .convection import coefficient_from_nusselt
    return coefficient_from_nusselt(Nu, k, D_h)
