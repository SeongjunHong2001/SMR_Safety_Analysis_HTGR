"""유체·열 포트 접속식. 연결식은 전역 Residual에 포함된다."""
import numpy as np
from interfaces.ports import FluidPort, FluidStreamPort, ThermalPort

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


def ideal_fluid_stream_junction_residual(
    y, registry, ports: tuple[FluidStreamPort, ...], *,
    mass_flow_epsilon: float,
):
    """무손실 다분기 접속의 압력·질량·유입 stream 조건을 반환한다.

    포트 유량은 구성요소 유입이 양수다. 따라서 음수 유량인 포트가
    분기점에 유체를 공급한다. 각 포트의 유입 엔탈피는 자기 포트를
    제외한 공급 포트들의 출구 엔탈피를 유량으로 가중한 값이다.
    epsilon은 정체·역류 시 혼합식을 매끄럽게 하는 kg/s 척도다.
    """
    ports = tuple(ports)
    if len(ports) < 2 or any(not isinstance(port, FluidStreamPort) for port in ports):
        raise ValueError("stream junction needs at least two FluidStreamPort instances")
    if not np.isfinite(mass_flow_epsilon) or mass_flow_epsilon <= 0:
        raise ValueError("mass_flow_epsilon must be positive and finite")

    get = lambda name: float(registry.get(y, name))
    pressures = np.array([get(port.p) for port in ports])
    flows = np.array([get(port.m_dot) for port in ports])
    outflow = np.array([get(port.h_outflow) for port in ports])
    instream = np.array([get(port.h_instream) for port in ports])

    # max(-m_dot, 0)의 매끄러운 근사. 양의 유량에서는 상쇄 오차를 피한다.
    magnitude = np.hypot(flows, mass_flow_epsilon)
    incoming = np.empty_like(flows)
    positive = flows >= 0
    incoming[positive] = (
        0.5 * mass_flow_epsilon
        * (mass_flow_epsilon / (magnitude[positive] + flows[positive]))
    )
    incoming[~positive] = 0.5 * (magnitude[~positive] - flows[~positive])
    mixed = np.empty(len(ports))
    for index in range(len(ports)):
        other = np.arange(len(ports)) != index
        weights = incoming[other]
        largest = weights.max()
        if largest == 0:
            mixed[index] = np.mean(outflow[other])
        else:
            normalized = weights / largest
            mixed[index] = np.dot(normalized, outflow[other]) / normalized.sum()

    return np.concatenate((
        pressures[1:] - pressures[0],
        [flows.sum()],
        instream - mixed,
    ))
