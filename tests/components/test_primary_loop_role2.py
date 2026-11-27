"""Role 2 플레넘·배관·순환기 단독 시험과 Role 2 장치만으로 닫은 1차 루프 시험.

실제 geometry.yaml 노드(하부·상부 플레넘, RV 입구, riser, cold annulus, hot duct)와 노심 유체를 StreamJunction으로
연결하고, 아직 없는 SG(Role 4)와 미확보 노드(110, 340: M10·M11)는 시험 전용 배관으로 대신한다.
SG 대용 배관은 일정 온도 물(193 °C)로 열을 빼는 SG-0 모양 열침원이다. 실제 HTGR 정상상태 결과가 아니다.
"""
from pathlib import Path

import numpy as np
import pytest
import yaml

from interfaces.component import Component
from interfaces.contracts import HeliumInventory
from interfaces.state_registry import StateRegistry
from models._inputs import ModelInputMissing
from models.compressor import Compressor
from models.compressor.performance import QuadraticPerformance
from models.core.fluid import ParallelChannelCoreFluid
from models.piping import Piping
from models.plenum import Plenum
from physics.properties import make_helium_properties
from solvers import DASSL, solve_steady_state
from solvers.initialization import consistent_initialize
from solvers.jacobian import finite_difference_jacobian
from system.assembly import SystemAssembly
from system.connections import StreamJunction
from system.inventory import total_helium_inventory

ROOT = Path(__file__).resolve().parents[2]
HE = make_helium_properties()
GEO = yaml.safe_load((ROOT / "config/geometry.yaml").open(encoding="utf-8"))
OP = yaml.safe_load((ROOT / "config/operation.yaml").open(encoding="utf-8"))
NODES = GEO["flow_network"]["nodes"]
M, P0, T_COLD, T_HOT, T_WATER = 157.1, 6.3373e6, 532.15, 961.0, 466.15


def ini(T, m=M):
    return dict(m_dot=m, T=T, p=P0)


def junction(name, a, b):
    return StreamJunction(name, (a, b), mass_flow_epsilon=1e-2, pressure_scale=1e3, flow_scale=M,
                          enthalpy_scale=5e5)


class SystemPressure(Component):
    """단독 시험 전용: p_sys 이완식 ṗ + (p − p0)/τ = 0 (폐루프에서는 하부 플레넘이 p_sys를 소유)."""

    name = "test_psys"

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=P0, scale=1e5, owner="primary_loop")

    def residual(self, t, y, ydot, context):
        return np.array([float(context.get(ydot, "p_sys")) + (float(context.get(y, "p_sys")) - P0)])

    def residual_ids(self):
        return ("test_psys.relax",)

    def residual_scales(self):
        return np.array([1e5])


class PortBC(Component):
    """단독 시험 전용: 포트 하나에 압력 또는 유량, 그리고 유입 엔탈피를 지정."""

    def __init__(self, name, port, *, p=None, m=None, T_in):
        self.name, self.port, self.p, self.m = name, port, p, m
        self.h = float(HE.enthalpy(P0, T_in))

    def register(self, registry):
        pass

    def residual(self, t, y, ydot, context):
        g = lambda k: float(context.get(y, k))  # noqa: E731
        first = g(self.port.p) - self.p if self.p is not None else g(self.port.m_dot) - self.m
        return np.array([first, g(self.port.h_instream) - self.h])

    def residual_ids(self):
        return (f"{self.name}.a", f"{self.name}.b")

    def residual_scales(self):
        return np.array([1e3 if self.p is not None else M, 1e6])


# ------------------------------------------------------------------ 단독 시험
def test_legacy_constructors_fail_with_missing_input_message():
    for cls in (Piping, Plenum, Compressor):
        with pytest.raises(ModelInputMissing, match="required"):
            cls().register(StateRegistry())
    node = dict(NODES["mc_discharge_plenum"], volume_m3=None)                   # 체적 없으면 키 이름으로 실패
    with pytest.raises(ModelInputMissing, match=r"flow_network\.nodes\.mc_discharge_plenum\.volume_m3"):
        Plenum(node, HE, name="mc_discharge_plenum", initial=ini(T_COLD))
    with pytest.raises(ModelInputMissing, match="I_circ_rotor"):
        Compressor(OP["compressor"], HE, mode="coastdown_inertia", m_rated=M, initial=ini(T_COLD)).register(
            StateRegistry())
    with pytest.raises(ModelInputMissing, match="performance map"):
        Compressor(OP["compressor"], HE, mode="performance", m_rated=M, initial=ini(T_COLD)).register(
            StateRegistry())


def test_piping_steady_pressure_drop_and_heat():
    # riser 12.74 m 상향: 정상 압력강하 = 마찰 + 정수압(ρ g Δz) + 플레넘 경계
    riser = Piping(NODES["riser"], HE, name="riser", initial=ini(T_COLD), inlet=("plenum", 0.5, 1.0),
                   outlet=("plenum", 1.0, 0.5))
    p = riser.ports()
    system = SystemAssembly([SystemPressure(), riser,
                             PortBC("bc_in", p["helium_inlet"], m=M, T_in=T_COLD),
                             PortBC("bc_out", p["helium_outlet"], p=P0, T_in=T_COLD)])
    y, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                              residual_scale=system.residual_scales())
    out = system.outputs(0.0, y)
    rho = float(HE.density(P0, HE.enthalpy(P0, T_COLD)))
    hydro = rho * 9.80665 * 12.74
    assert hydro < out["riser.dp_inlet_to_outlet"] < hydro + 6.0e3   # 마찰 몇 kPa 이내 (PSID riser 설계 4.8 kPa)
    assert isinstance(riser, HeliumInventory)
    assert riser.helium_volumes().sum() == pytest.approx(NODES["riser"]["volume_m3"], rel=1e-3)


def test_plenum_mixing_and_hydrostatic_ports():
    # 3포트 플레넘: 두 유입(300 K, 600 K)이 섞여 한 출구로 나감
    cfg = {"volume_m3": 10.0, "z_bottom_m": 0.0, "z_top_m": 2.0}
    plen = Plenum(cfg, HE, name="mix", initial=ini(450.0, m=10.0), port_elevations=[2.0, 2.0, 0.0],
                  port_labels=("a", "b", "out"))
    system = SystemAssembly([SystemPressure(), plen,
                             PortBC("bc_a", plen.port("a"), m=4.0, T_in=300.0),
                             PortBC("bc_b", plen.port("b"), m=6.0, T_in=600.0),
                             PortBC("bc_out", plen.port("out"), p=P0, T_in=450.0)])
    y, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                              residual_scale=system.residual_scales())
    T_mix = system.outputs(0.0, y)["mix.T_he"]
    assert T_mix == pytest.approx((4 * 300 + 6 * 600) / 10, abs=1e-3)
    reg = system.registry
    assert reg.get(y, plen.port("out").m_dot) == pytest.approx(-10.0)
    rho = float(HE.density(P0, reg.get(y, "mix.h_he")))
    # 아래 포트(z 0)가 위 포트(z 2)보다 ρ g Δz만큼 높은 압력
    assert reg.get(y, plen.port("out").p) - reg.get(y, plen.port("a").p) == pytest.approx(rho * 9.80665 * 2.0)


def test_compressor_prescribed_flow_and_heating():
    comp = Compressor(OP["compressor"], HE, mode="prescribed_flow", m_dot_set=M, m_rated=M, initial=ini(T_COLD))
    p = comp.ports()
    system = SystemAssembly([SystemPressure(), comp,
                             PortBC("bc_in", p["helium_inlet"], p=6.29e6, T_in=528.15),
                             PortBC("bc_out", p["helium_outlet"], p=6.38e6, T_in=T_COLD)])
    y, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                              residual_scale=system.residual_scales())
    out = system.outputs(0.0, y)
    assert out["main_circulator.pressure_rise"] == pytest.approx(9.0e4)
    assert out["main_circulator.Q_to_helium"] == pytest.approx(3.1e6, rel=1e-6)
    T_dis = float(HE.temperature(P0, system.registry.get(y, p["helium_outlet"].h_outflow)))
    assert T_dis - 528.15 == pytest.approx(3.1e6 / (M * HE.cp_const), rel=1e-6)   # 약 3.8 K (PSID)


def test_quadratic_performance_requires_consistent_coefficients():
    perf = QuadraticPerformance(rated_speed=649.26, rated_mass_flow=M, rated_density=5.73, rated_pressure_rise=9.1e4,
                                c0=1.3, c1=0.0, c2=0.3, efficiency=0.8)   # 변수표 예비 계수 (대리값)
    rated = perf.evaluate(speed=649.26, mass_flow=M, inlet_state={"rho": 5.73}, geometry={})
    assert rated["pressure_rise"] == pytest.approx(9.1e4)
    half = perf.evaluate(speed=324.63, mass_flow=M / 2, inlet_state={"rho": 5.73}, geometry={})
    assert half["pressure_rise"] == pytest.approx(9.1e4 / 4)                # 상사법칙 ΔP ∝ ω²
    # 정지 회전자(α = 0)에서는 유동 저항이 된다: ΔP = −c2 ν|ν| ΔP_rated
    locked = perf.evaluate(speed=0.0, mass_flow=M, inlet_state={"rho": 5.73}, geometry={})
    assert locked["pressure_rise"] == pytest.approx(-0.3 * 9.1e4)
    with pytest.raises(ValueError, match="c0 \\+ c1 - c2"):
        QuadraticPerformance(rated_speed=1, rated_mass_flow=1, rated_density=1, rated_pressure_rise=1,
                             c0=1.0, c1=0.5, c2=0.0, efficiency=0.8)


# ------------------------------------------------------------------ 폐루프
class HeatedCore(Component):
    """시험 전용: 노심 유체에 지정 열률 (Role 3 고체·결합 대신)."""

    def __init__(self, fluid, q_of_t):
        self.fluid, self.q_of_t = fluid, q_of_t
        self.name = "core_test"

    def register(self, registry):
        self.fluid.register(registry)

    def residual(self, t, y, ydot, context):
        return self.fluid.residual(t, y, ydot, context, heat_rate=self.q_of_t(t))

    def ports(self):
        return self.fluid.ports()

    def residual_ids(self):
        return self.fluid.residual_ids()

    def residual_scales(self):
        return self.fluid.residual_scales()

    def outputs(self, t, y, context):
        return {**self.fluid.outputs(t, y, context), "core_test.Q": float(np.sum(self.q_of_t(t)))}


def build_loop(*, steady, core_power=lambda t: 1.0):
    UA = 1.7e6                                               # [시험 가정] SG-0 모양 열침원 UA [W/K]

    def sink(t, y, context, pipe):
        p_sys = float(context.get(y, "p_sys"))
        T = np.asarray(HE.temperature(p_sys, context.get(y, pipe.key("h_he"))), dtype=float)
        return -UA / pipe.n * (T - T_WATER)

    synth = lambda n, A, D, L, z0, z1: dict(cells=n, flow_area_m2=A, hydraulic_diameter_m=D, length_m=L,  # noqa
                                            z_inlet_m=z0, z_outlet_m=z1)
    lower = Plenum(NODES["lower_plenum"], HE, name="lower_plenum", initial=ini(T_HOT), owns_p_sys=True,
                   steady_pressure_closure=P0 if steady else None, port_labels=("from_core", "to_hot_duct"))
    upper = Plenum(NODES["upper_plenum"], HE, name="upper_plenum", initial=ini(T_COLD),
                   port_labels=("from_riser", "to_core"))
    rv_in = Plenum(NODES["rv_inlet_mcss"], HE, name="rv_inlet_mcss", initial=ini(T_COLD),
                   port_labels=("from_cold_annulus", "to_riser"))
    hot = Piping(NODES["hot_duct"], HE, name="hot_duct", initial=ini(T_HOT), inlet=("plenum", 0.5, 1.0))
    sg = Piping(synth(4, 9.23, 0.05, 14.08, 1.58, -12.5), HE, name="sg_standin", initial=ini(700.0),
                heat_source=sink)                            # [시험 가정] SG He 측 대용 (Role 4)
    ret = Piping(synth(3, 2.0, 0.2, 20.7, -12.5, 8.2), HE, name="annulus_standin", initial=ini(T_COLD))  # 340 대용
    comp = Compressor(OP["compressor"], HE, mode="prescribed_flow", m_dot_set=M, m_rated=M, initial=ini(T_COLD))
    dis = Piping(synth(2, 1.0, 0.5, 6.62, 8.2, 1.58), HE, name="discharge_standin", initial=ini(T_COLD))  # 110 대용
    cold = Piping(NODES["cold_annulus"], HE, name="cold_annulus", initial=ini(T_COLD), outlet=("plenum", 1.0, 0.5))
    riser = Piping(NODES["riser"], HE, name="riser", initial=ini(T_COLD), inlet=("plenum", 0.5, 1.0),
                   outlet=("plenum", 1.0, 0.5))
    fluid = ParallelChannelCoreFluid(GEO, HE, channels=("average", "hot"),
                                     initial=dict(m_dot=M, T_inlet=T_COLD, T_outlet=T_HOT, p_inlet=P0))
    dist = OP["reactor"]["power_distribution"]
    axial = np.asarray(dist["axial_fraction"], dtype=float)
    shares = np.array([dist["channel_fraction"]["average"], dist["channel_fraction"]["hot"]])
    q0 = np.concatenate([3.5e8 * s / shares.sum() * axial for s in shares])
    core = HeatedCore(fluid, lambda t: q0 * core_power(t))

    fp = fluid.ports()
    links = [junction("J_core_out", fp["helium_outlet"], lower.port("from_core")),
             junction("J12", lower.port("to_hot_duct"), hot.ports()["helium_inlet"]),
             junction("J13", hot.ports()["helium_outlet"], sg.ports()["helium_inlet"]),
             junction("J17", sg.ports()["helium_outlet"], ret.ports()["helium_inlet"]),
             junction("J19", ret.ports()["helium_outlet"], comp.ports()["helium_inlet"]),
             junction("J01", comp.ports()["helium_outlet"], dis.ports()["helium_inlet"]),
             junction("J02", dis.ports()["helium_outlet"], cold.ports()["helium_inlet"]),
             junction("J03", cold.ports()["helium_outlet"], rv_in.port("from_cold_annulus")),
             junction("J04", rv_in.port("to_riser"), riser.ports()["helium_inlet"]),
             junction("J05", riser.ports()["helium_outlet"], upper.port("from_riser")),
             junction("J06", upper.port("to_core"), fp["helium_inlet"])]
    parts = [lower, upper, rv_in, hot, sg, ret, comp, dis, cold, riser, core]
    system = SystemAssembly(parts, connections=links)
    inventory = [lower, upper, rv_in, hot, sg, ret, dis, cold, riser, fluid]
    return system, inventory, dict(sg=sg, comp=comp, core=core, fluid=fluid, lower=lower)


@pytest.fixture(scope="module")
def loop_steady():
    system, inventory, parts = build_loop(steady=True)
    y, ydot = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                                 residual_scale=system.residual_scales())
    return system, inventory, parts, y, ydot


def test_closed_loop_has_single_p_sys_and_is_square(loop_steady):
    system, inventory, parts, y, ydot = loop_steady
    reg = system.registry
    assert reg.spec("p_sys").owner == "primary_loop" and reg.spec("p_sys").kind == "differential"
    assert sum(name == "p_sys" for name in reg.names) == 1
    assert len(system.residual_ids()) == reg.size
    assert "lower_plenum.pressure_reference" in system.residual_ids()
    assert "lower_plenum.steady_pressure_closure" in system.residual_ids()


def test_closed_loop_steady_state_energy_and_pressure(loop_steady):
    system, inventory, parts, y, ydot = loop_steady
    # Newton 목표는 초기 잔차에 상대적이다 (rtol 1e-9 × 초기 노름). 척도 잔차 1e-6 이하면 수렴으로 본다.
    assert np.linalg.norm(system.residual(0.0, y, ydot) / system.residual_scales(), ord=np.inf) < 1e-6
    out = system.outputs(0.0, y)
    reg = system.registry
    assert reg.get(y, "p_sys") == pytest.approx(P0)
    assert reg.get(y, "lower_plenum.p_fluid_cell") == pytest.approx(P0)          # 국부 기준식
    q_core, q_comp, q_sink = out["core_test.Q"], out["main_circulator.Q_to_helium"], out["sg_standin.Q_wall"]
    assert q_core + q_comp + q_sink == pytest.approx(0.0, abs=1e-6 * q_core)   # 정상 열수지
    assert out["main_circulator.m_dot"] == pytest.approx(M)
    # 순환기 승압 = 루프 비가역 손실 (+ 정수압 합 ≈ 0). 손실계수 미보정이라 PSID 91 kPa보다 작다.
    dp = out["main_circulator.pressure_rise"]
    assert 1.5e4 < dp < 9.1e4
    # 노심 입구는 순환기 출구 엔탈피 (승온 약 3.8 K)
    T_core_in = float(HE.temperature(P0, reg.get(y, parts["fluid"].ports()["helium_inlet"].h_instream)))
    T_sink_out = out["sg_standin.T_he"][-1]
    assert T_core_in - T_sink_out == pytest.approx(3.1e6 / (M * HE.cp_const), rel=1e-3)
    # 저마하수 재고량
    total = total_helium_inventory(inventory, 0.0, y, system.context)
    assert 400.0 < total < 2585.0


def test_closed_loop_transient_conserves_mass_and_energy(loop_steady):
    system_s, _, _, y_s, _ = loop_steady
    system, inventory, parts = build_loop(steady=False,
                                          core_power=lambda t: 1.0 + 0.1 * min(max(t, 0.0) / 0.5, 1.0))
    mask = system.registry.differential_mask()
    y0, yp0 = consistent_initialize(system.residual, 0.0, y_s, mask, residual_scale=system.residual_scales())
    # 희소 Jacobian 패턴: 정상해에서 ∂F/∂y와 ∂F/∂ẏ의 0이 아닌 원소 (Role 5 builder 참고)
    Jy = finite_difference_jacobian(lambda z: system.residual(0.0, z, yp0), y0)
    Jd = finite_difference_jacobian(lambda z: system.residual(0.0, y0, z), yp0)
    pattern = (Jy != 0) | (Jd != 0) | np.eye(y0.size, dtype=bool)
    result = DASSL(dt=0.05, rtol=1e-6, atol=1e-6, max_step=0.5, residual_scale=system.residual_scales()).solve(
        system.residual, t_span=(0.0, 2.0), y0=y0, ydot0=yp0, differential_mask=mask, jacobian_sparsity=pattern)
    assert result.t[-1] == pytest.approx(2.0)
    reg = system.registry
    M_he = np.array([total_helium_inventory(inventory, t, s, system.context) for t, s in zip(result.t, result.y)])
    assert np.max(np.abs(M_he - M_he[0])) < 1e-7 * M_he[0]                    # 닫힌 루프 질량 보존

    def internal_energy(t, s):
        p_sys = float(reg.get(s, "p_sys"))
        U = 0.0
        for comp in inventory:
            V = comp.helium_volumes()
            mass = comp.helium_inventory(t, s, system.context)
            keys = ([comp._ck(ch, "h_he") for ch in comp.channels] if hasattr(comp, "channels")
                    else [comp.key("h_he")])
            h = np.concatenate([np.atleast_1d(reg.get(s, k)) for k in keys])
            U += float(np.sum(mass * h - V * p_sys))
        return U

    U = np.array([internal_energy(t, s) for t, s in zip(result.t, result.y)])
    outs = [system.outputs(t, s) for t, s in zip(result.t, result.y)]
    q_net = np.array([o["core_test.Q"] + o["main_circulator.Q_to_helium"] + o["sg_standin.Q_wall"] for o in outs])
    integral = np.concatenate(([0.0], np.cumsum(0.5 * (q_net[1:] + q_net[:-1]) * np.diff(result.t))))
    assert np.max(np.abs((U - U[0]) - integral)) < 0.02 * np.max(np.abs(integral))   # dU/dt = Σ Q (사다리꼴)
    assert reg.get(result.y[-1], "p_sys") > reg.get(result.y[0], "p_sys")          # 가열 → 계통압 상승
