"""유체·열 포트 접속식. 연결식은 전역 Residual에 포함된다."""
from dataclasses import dataclass
import numpy as np
from interfaces.ports import FluidPort, FluidStreamPort, ThermalPort, ThermalSurfacePort


def validate_port(port, registry):
    """Resolve symbolic keys only after all components have registered."""
    if isinstance(port, FluidStreamPort):
        if not isinstance(port.medium, str) or not port.medium:
            raise ValueError("Fluid medium must be a nonempty name")
        fields = (("p", "Pa"), ("m_dot", "kg/s"),
                  ("h_outflow", "J/kg"), ("h_instream", "J/kg"))
    elif isinstance(port, ThermalPort):
        fields = (("T", "K"), ("Q_dot", "W"))
    elif isinstance(port, ThermalSurfacePort):
        fields = (("T", "K"),)
    elif isinstance(port, FluidPort):  # old two-port demonstration API
        fields = (("p", None), ("h", None), ("m_dot", None))
    else:
        raise TypeError(f"Unsupported port type: {type(port).__name__}")
    for field, unit in fields:
        key = getattr(port, field)
        try:
            spec = registry.spec(key)
        except KeyError as exc:
            raise ValueError(f"Port {field} references missing Registry key {key}") from exc
        if isinstance(port, ThermalSurfacePort):
            index = port.index
            if index is None and spec.shape:
                raise ValueError(f"Port {field} key {key} needs an array index")
            if index is not None and (len(index) != len(spec.shape)
                                      or any(i >= n for i, n in zip(index, spec.shape))):
                raise ValueError(f"Port {field} key {key} has an invalid array index")
        elif spec.shape:
            raise ValueError(f"Port {field} key {key} must be scalar")
        if unit is not None and spec.unit != unit:
            raise ValueError(f"Port {field} key {key}: expected {unit}, got {spec.unit}")
        if ((isinstance(port, FluidStreamPort) and field == "h_instream") or
                (isinstance(port, ThermalPort) and field == "Q_dot")) and spec.kind != "algebraic":
            raise ValueError(f"Port {field} key {key} must be algebraic")


def validate_stream_junction(ports, registry):
    ports = tuple(ports)
    if len(ports) < 2 or any(not isinstance(port, FluidStreamPort) for port in ports):
        raise ValueError("stream junction needs at least two FluidStreamPort instances")
    media = [port.medium for port in ports]
    if any(not isinstance(medium, str) or not medium for medium in media):
        raise ValueError("Fluid medium must be a nonempty name")
    if len(set(media)) != 1:
        raise ValueError(f"Cannot connect different fluid media: {sorted(media)}")
    for port in ports:
        validate_port(port, registry)
    for field in ("p", "m_dot", "h_instream"):
        keys = [registry.resolve(getattr(port, field)) for port in ports]
        if len(set(keys)) != len(keys):
            raise ValueError(f"Junction has duplicate {field} state keys")
    instream_keys = {registry.resolve(port.h_instream) for port in ports}
    outflow_keys = {registry.resolve(port.h_outflow) for port in ports}
    if instream_keys & outflow_keys:
        raise ValueError("Junction h_instream must be independent of h_outflow")
    return ports


@dataclass(frozen=True)
class StreamJunction:
    """Named 2N-equation junction with explicit physical residual scales."""

    name: str
    ports: tuple[FluidStreamPort, ...]
    mass_flow_epsilon: float
    pressure_scale: float
    flow_scale: float
    enthalpy_scale: float

    def __post_init__(self):
        object.__setattr__(self, "ports", tuple(self.ports))
        if not self.name or len(self.ports) < 2:
            raise ValueError("Junction needs a name and at least two ports")
        for value in (self.mass_flow_epsilon, self.pressure_scale,
                      self.flow_scale, self.enthalpy_scale):
            if not np.isfinite(value) or value <= 0:
                raise ValueError("Junction epsilon and residual scales must be positive")

    def validate(self, registry):
        validate_stream_junction(self.ports, registry)

    def __call__(self, y, registry):
        return ideal_fluid_stream_junction_residual(
            y, registry, self.ports, mass_flow_epsilon=self.mass_flow_epsilon,
        )

    def residual_ids(self):
        n = len(self.ports)
        return tuple([f"{self.name}.pressure[{i}]" for i in range(1, n)]
                     + [f"{self.name}.mass"]
                     + [f"{self.name}.instream[{i}]" for i in range(n)])

    def residual_scales(self):
        n = len(self.ports)
        return np.array([self.pressure_scale] * (n - 1)
                        + [self.flow_scale]
                        + [self.enthalpy_scale] * n)

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
    ports = validate_stream_junction(ports, registry)
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
