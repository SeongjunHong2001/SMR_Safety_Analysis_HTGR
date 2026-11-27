"""S0 core_benchmark 조립 원형 (시험 전용 builder).

config/scenarios/steady_state.yaml ``conditions.core_benchmark``의 계산 영역(노드 130·140·150·200·210·220·230,
HS-01~HS-08)과 경계(입구 J03 유량·온도, 출구 J12 압력, RCCS-0 30 °C)를 실제 Role 2·3·4 Component와
StreamJunction·ThermalLink로 손으로 조립한다. Role 5 builder가 YAML에서 같은 조립을 하기 전까지의 원형이다.

[시험 가정] 우회 유로 수력직경 4 mm 평행평판 (geometry null), 상·하부 헤드(HS-09·10) 단열,
p_sys는 출구 압력 6.39 MPa로 고정(개방 영역이라 하부 플레넘이 p_sys를 소유하지 않음).
"""
from pathlib import Path
import math

import numpy as np
import yaml

from interfaces.component import Component
from models.core import Core
from models.core.coupling import UnitCellCoreCoupling
from models.core.fluid import ParallelChannelCoreFluid
from models.core.power import SpecifiedCorePower
from models.core.solid import UnitCellCoreSolid
from models.piping import Piping
from models.plenum import Plenum
from models.rccs.rccs import RCCS0Sink
from models.vessel import Vessel
from models.vessel.vessel import SurfaceFluidConvection, helium_gap_conductance
from physics.heat_transfer.correlations import make_convection_correlation
from physics.heat_transfer.solid_properties import R1SolidProperties
from physics.properties import make_helium_properties
from system.assembly import SystemAssembly
from system.connections import StreamJunction

ROOT = Path(__file__).resolve().parents[2]
BYPASS_OVERRIDE = {"bypass": {"hydraulic_diameter_m": 0.004, "shape": "parallel_plates"}}


def load_inputs():
    load = lambda p: yaml.safe_load((ROOT / p).open(encoding="utf-8"))  # noqa: E731
    return (load("config/geometry.yaml"), load("config/operation.yaml"),
            load("config/scenarios/steady_state.yaml")["conditions"]["core_benchmark"])


class FixedSystemPressure(Component):
    """개방 영역의 열역학 기준압: ṗ_sys + (p_sys − p_out)/τ = 0 (정상 p_sys = 출구 압력, 과도 일정)."""

    name = "bc_p_sys"

    def __init__(self, p):
        self.p = float(p)

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=self.p, scale=1e5, owner="primary_loop")

    def residual(self, t, y, ydot, context):
        return np.array([float(context.get(ydot, "p_sys")) + (float(context.get(y, "p_sys")) - self.p)])

    def residual_ids(self):
        return ("bc_p_sys.relax",)

    def residual_scales(self):
        return np.array([1e5])


class PortBoundary(Component):
    """포트 하나에 유량 또는 압력, 그리고 들어오는 stream 엔탈피를 지정."""

    def __init__(self, name, port, *, h_in, m=None, p=None):
        self.name, self.port, self.h_in, self.m, self.p = name, port, float(h_in), m, p

    def register(self, registry):
        pass

    def residual(self, t, y, ydot, context):
        g = lambda k: float(context.get(y, k))  # noqa: E731
        first = g(self.port.m_dot) - self.m if self.m is not None else g(self.port.p) - self.p
        return np.array([first, g(self.port.h_instream) - self.h_in])

    def residual_ids(self):
        return (f"{self.name}.flow_or_pressure", f"{self.name}.instream")

    def residual_scales(self):
        return np.array([100.0 if self.m is not None else 1e3, 1e6])


def build(power_W=None, *, geo=None, op=None, s0=None):
    geo, op, s0 = (geo, op, s0) if geo is not None else load_inputs()
    he = make_helium_properties()
    sp = R1SolidProperties()
    corr = make_convection_correlation("mceligot")
    nodes = geo["flow_network"]["nodes"]
    m_in = float(s0["inlet_boundary"]["mass_flow_kg_s"])
    T_in = float(s0["inlet_boundary"]["temperature_K"])
    p_out = float(s0["outlet_boundary"]["pressure_Pa"])
    Q = float(s0["reactor"]["thermal_power_W"]) if power_W is None else power_W
    T_out = T_in + (Q if not callable(Q) else Q(0.0)) / (m_in * 5195.0)
    ini = lambda T: dict(m_dot=m_in, T=T, p=p_out)  # noqa: E731

    rv_in = Plenum(nodes["rv_inlet_mcss"], he, name="rv_inlet_mcss", initial=ini(T_in),
                   port_labels=("from_cold_annulus", "to_riser"))
    upper = Plenum(nodes["upper_plenum"], he, name="upper_plenum", initial=ini(T_in),
                   port_labels=("from_riser", "to_core"))
    lower = Plenum(nodes["lower_plenum"], he, name="lower_plenum", initial=ini(T_out),
                   port_labels=("from_core", "to_hot_duct"))

    def riser_heat(t, y, context, pipe):
        return context.heat_rates["barrel_riser.to_b"] + context.heat_rates["riser_rpv.to_b"]

    riser = Piping(nodes["riser"], he, name="riser", initial=ini(T_in), inlet=("plenum", 0.5, 1.0),
                   outlet=("plenum", 1.0, 0.5), heat_source=riser_heat)

    fluid = ParallelChannelCoreFluid(geo, he, initial=dict(m_dot=m_in, T_inlet=T_in, T_outlet=T_out, p_inlet=p_out),
                                     mode_bypass=s0["bypass_model"], channel_overrides=BYPASS_OVERRIDE)
    z_core = np.asarray(geo["core"]["elevation_profile_m"], dtype=float)
    solid = UnitCellCoreSolid(geo, sp, helium=he, initial=dict(T=np.linspace(T_in, T_out, 17)),
                              outer_link="barrel_riser")
    power = SpecifiedCorePower(solid, op["reactor"]["power_distribution"], Q)
    coupling = UnitCellCoreCoupling(fluid, solid, corr,
                                    bypass_regions=geo["core"]["channels"]["bypass"]["type_I_fraction"])
    core = Core(fluid=fluid, solid=solid, power=power, coupling=coupling)

    rpv = Vessel(geo["rpv"], sp, name="rpv", z_faces=riser.z_faces, initial=dict(T=T_in - 50.0),
                 links=(("riser_rpv", "to_a"), ("rpv_rccs", "to_a")))
    rad = geo["core"]["radial_m"]
    r_bo, r_bi = float(rad["core_barrel_outer"]), float(rad["core_barrel_inner"])
    r_riser_o, r_rpv_i = float(rad["riser_outer"]), float(rad["rpv_inner"])
    dz_core = np.abs(np.diff(z_core))
    k_barrel = sp.conductivity("alloy_800H")
    barrel_half = 2 * math.pi * k_barrel * dz_core / math.log(r_bo / (0.5 * (r_bi + r_bo)))
    barrel_key = solid._state_key(("core_barrel",))
    barrel_riser = SurfaceFluidConvection(
        "barrel_riser", solid_T=lambda y, c: c.get(y, barrel_key), z_solid=z_core,
        area=2 * math.pi * r_bo * dz_core, pipe=riser, correlation=corr,
        series=lambda Ts, Tf, p: barrel_half)
    gap = helium_gap_conductance(he, r_riser_o, r_rpv_i, rpv.dz, sp.emissivity("alloy_800H"),
                                 sp.emissivity("pressure_vessel"))
    riser_rpv = SurfaceFluidConvection(
        "riser_rpv", solid_T=lambda y, c: c.get(y, rpv.T_key), z_solid=rpv.z_faces,
        area=2 * math.pi * r_riser_o * rpv.dz, pipe=riser, correlation=corr, series=gap)
    rccs = RCCS0Sink(geo["rccs"], rpv, sp, T_sink=s0["rccs"]["boundary_temperature_K"])

    fp = fluid.ports()
    junction = lambda name, a, b: StreamJunction(name, (a, b), mass_flow_epsilon=1e-2, pressure_scale=1e3,  # noqa
                                                  flow_scale=m_in, enthalpy_scale=5e5)
    links = [junction("J04", rv_in.port("to_riser"), riser.ports()["helium_inlet"]),
             junction("J05", riser.ports()["helium_outlet"], upper.port("from_riser")),
             junction("J06_J08", upper.port("to_core"), fp["helium_inlet"]),
             junction("J09_J11", fp["helium_outlet"], lower.port("from_core"))]
    bcs = [FixedSystemPressure(p_out),
           PortBoundary("bc_J03", rv_in.port("from_cold_annulus"), m=m_in, h_in=he.enthalpy(p_out, T_in)),
           PortBoundary("bc_J12", lower.port("to_hot_duct"), p=p_out, h_in=he.enthalpy(p_out, T_out))]
    system = SystemAssembly(bcs + [rv_in, riser, upper, core, lower, rpv], connections=links,
                            thermal_links=[barrel_riser.link(), riser_rpv.link(), rccs.link()])
    parts = dict(he=he, fluid=fluid, solid=solid, power=power, coupling=coupling, riser=riser, rpv=rpv, rccs=rccs,
                 lower=lower, upper=upper, m_in=m_in, T_in=T_in, p_out=p_out)
    return system, parts
