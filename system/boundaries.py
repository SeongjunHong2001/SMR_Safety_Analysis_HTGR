"""개방 계산 영역의 경계 Component (Role 5 builder 전용).

폐루프가 아닌 영역(S0 core_benchmark)은 하부 플레넘이 p_sys를 소유하지 않으므로 열역학 기준압과
잘린 연결부의 유량·압력·유입 stream을 경계로 준다. 폐루프 builder는 이 클래스를 쓰지 않는다.
"""
from __future__ import annotations

import numpy as np

from interfaces.component import Component


class FixedSystemPressure(Component):
    """개방 영역의 열역학 기준압: ṗ_sys + (p_sys − p_ref)/τ = 0.

    정상해에서 p_sys = p_ref이고, 과도에서도 p_sys가 p_ref에 붙어 있다 (τ = 1 s, 영역 밖 계통이
    압력을 잡아 준다고 보는 경계 [가정]). p_sys를 대수식으로 고정하면 미분 상태가 대수적으로
    묶여 일치 초기화가 실패하므로 이완식으로 둔다.
    """

    def __init__(self, p_ref, *, name="bc_p_sys", tau_s=1.0, p_sys_key="p_sys"):
        self.name = name
        self.p_ref = float(p_ref)
        self.tau = float(tau_s)
        self.p_sys_key = p_sys_key
        if not self.p_ref > 0 or not self.tau > 0:
            raise ValueError("FixedSystemPressure needs positive p_ref and tau_s")

    def register(self, registry):
        registry.register(self.p_sys_key, kind="differential", unit="Pa", initial=self.p_ref, scale=1e5,
                          owner="primary_loop")

    def residual(self, t, y, ydot, context):
        p = float(context.get(y, self.p_sys_key))
        return np.array([float(context.get(ydot, self.p_sys_key)) + (p - self.p_ref) / self.tau])

    def residual_ids(self):
        return (f"{self.name}.relax",)

    def residual_scales(self):
        return np.array([1e5])


class PortBoundary(Component):
    """잘린 연결부의 포트 하나에 유량 또는 압력, 그리고 들어오는 stream 엔탈피를 지정한다 (잔차 2개).

    ``FluidStreamPort``의 포트 잔차 2개(압력/유량 중 하나, h_instream)를 StreamJunction 대신 준다.
    """

    def __init__(self, name, port, *, h_in, m_dot=None, p=None):
        if (m_dot is None) == (p is None):
            raise ValueError(f"{name}: give exactly one of m_dot or p")
        self.name, self.port, self.h_in = name, port, float(h_in)
        self.m_dot = None if m_dot is None else float(m_dot)
        self.p = None if p is None else float(p)

    def register(self, registry):
        pass

    def connected_ports(self):
        return (self.port,)

    def residual(self, t, y, ydot, context):
        g = lambda k: float(context.get(y, k))  # noqa: E731
        first = g(self.port.m_dot) - self.m_dot if self.m_dot is not None else g(self.port.p) - self.p
        return np.array([first, g(self.port.h_instream) - self.h_in])

    def residual_ids(self):
        return (f"{self.name}.{'mass_flow' if self.m_dot is not None else 'pressure'}", f"{self.name}.instream")

    def residual_scales(self):
        return np.array([100.0 if self.m_dot is not None else 1e3, 1e6])
