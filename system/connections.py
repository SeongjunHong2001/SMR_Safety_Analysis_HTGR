"""간단한 2-Port 이상적 접속 예시 (모델 선정 후 수정 필요).

연결식은 단독 시간적분을 하지 않고 전역 Residual에 포함되어야 한다.
다분기 유량 혼합 및 stream 변수는 이 파일에서 별도 구현 필요.
"""
import numpy as np
from interfaces.ports import FluidPort, ThermalPort

def ideal_fluid_two_port_residual(y, registry, a: FluidPort, b: FluidPort):
    """손실 없는 두 포트: 압력 동일, 유량합 0, 엔탈피 동일.

    *정확히 두 포트*가 연결되는 경우의 초기 예시. 실제 펌프/배관의
    압력손실은 Component 방정식으로 계산할 것.
    """
    get = lambda n: float(registry.get(y, n))
    return np.array([get(a.p)-get(b.p),
                     get(a.m_dot)+get(b.m_dot),
                     get(a.h)-get(b.h)])

def ideal_thermal_two_port_residual(y, registry, a: ThermalPort, b: ThermalPort):
    get = lambda n: float(registry.get(y, n))
    return np.array([get(a.T)-get(b.T), get(a.Q_dot)+get(b.Q_dot)])
