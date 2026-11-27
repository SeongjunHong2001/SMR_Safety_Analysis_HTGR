"""Role 2 노심 유체 (models/core/fluid.py) 단독 시험.

실제 geometry.yaml과 R1 헬륨 물성으로 병렬 유로 노심 유체를 조립하고, 시험 전용 경계(포트 조건·p_sys)와
지정 열률로 정상해·과도·역류를 푼다. 고체·결합(Role 3)과 플레넘(Role 2 후속)은 아직 없으므로
경계 Component는 이 시험에서만 쓴다. 실제 HTGR 정상상태 결과가 아니다.
"""
from pathlib import Path

import numpy as np
import pytest
import yaml

from interfaces.component import Component
from interfaces.contracts import HeliumInventory
from models.core import Core
from models.core.fluid import ParallelChannelCoreFluid, initial_from_operation
from physics.properties import make_helium_properties
from solvers import DASSL, solve_steady_state
from solvers.initialization import consistent_initialize
from system.assembly import SystemAssembly
from system.inventory import total_helium_inventory

ROOT = Path(__file__).resolve().parents[2]
HE = make_helium_properties()
GEO = yaml.safe_load((ROOT / "config/geometry.yaml").open(encoding="utf-8"))
OP = yaml.safe_load((ROOT / "config/operation.yaml").open(encoding="utf-8"))
# [시험 가정] 우회 유로 수력직경은 geometry에서 미정(null)이라 블록 간극 후보 4 mm 평행평판으로 둔다.
BYPASS_OVERRIDE = {"bypass": {"hydraulic_diameter_m": 0.004, "shape": "parallel_plates"}}
P_OUT = 6.3373e6


class Boundaries(Component):
    """시험 전용: p_sys와 플레넘 쪽 포트 조건 (입구 유량·유입 엔탈피, 출구 압력·역류 유입 엔탈피)."""

    def __init__(self, fluid_name, *, m_in, T_in, p_out, T_back, p_sys=P_OUT):
        self.name = "test_bc"
        self.f = fluid_name
        self.m_in, self.p_out, self.p_sys = m_in, p_out, p_sys
        self.h_in, self.h_back = float(HE.enthalpy(p_sys, T_in)), float(HE.enthalpy(p_sys, T_back))

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=self.p_sys, scale=1e5,
                          owner="primary_loop")

    def residual(self, t, y, ydot, context):
        g = lambda k: float(context.get(y, k))  # noqa: E731
        m_in = self.m_in(t) if callable(self.m_in) else self.m_in
        # p_sys는 미분 상태: 정상에서는 p_sys = p0, 과도에서는 ṗ_sys가 정해지도록 이완식 (τ = 1 s)
        return np.array([float(context.get(ydot, "p_sys")) + (g("p_sys") - self.p_sys) / 1.0,
                         g(f"{self.f}.helium_inlet.m_dot") - m_in,
                         g(f"{self.f}.helium_inlet.h_instream") - self.h_in,
                         g(f"{self.f}.helium_outlet.p") - self.p_out,
                         g(f"{self.f}.helium_outlet.h_instream") - self.h_back])

    def residual_ids(self):
        return ("test_bc.p_sys", "test_bc.m_in", "test_bc.h_in", "test_bc.p_out", "test_bc.h_back")

    def residual_scales(self):
        return np.array([1e5, 100.0, 1e6, 1e5, 1e6])


class HeatedFluid(Component):
    """시험 전용: 지정 열률 q(t) [W, cell_order]를 노심 유체에 넘긴다 (Role 3 결합 대신)."""

    def __init__(self, fluid, q):
        self.fluid, self.q = fluid, q
        self.name = fluid.name + "_test"

    def register(self, registry):
        self.fluid.register(registry)

    def residual(self, t, y, ydot, context):
        return self.fluid.residual(t, y, ydot, context, heat_rate=self.q(t))

    def ports(self):
        return self.fluid.ports()

    def residual_ids(self):
        return self.fluid.residual_ids()

    def residual_scales(self):
        return self.fluid.residual_scales()

    def outputs(self, t, y, context):
        return self.fluid.outputs(t, y, context)


def uniform_power(fluid, total_W):
    """operation.yaml 균일 출력분포: 연료 셀(축 분율 0.1)에 유로 분율대로."""
    dist = OP["reactor"]["power_distribution"]
    axial = np.asarray(dist["axial_fraction"], dtype=float)
    q = []
    for ch in fluid.channels:
        q.append(total_W * dist["channel_fraction"][ch.label] * axial)
    return np.concatenate(q)


def build(fluid, q, **bc):
    system = SystemAssembly([HeatedFluid(fluid, q), Boundaries(fluid.name, **bc)])
    return system


def steady(system):
    y, ydot = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                                 residual_scale=system.residual_scales())
    return y, ydot


@pytest.fixture(scope="module")
def rated():
    init = initial_from_operation(OP)
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channel_overrides=BYPASS_OVERRIDE)
    q = uniform_power(fluid, 3.5e8)
    system = build(fluid, lambda t: q, m_in=init["m_dot"], T_in=init["T_inlet"], p_out=P_OUT,
                   T_back=init["T_outlet"])
    y, ydot = steady(system)
    return fluid, system, y, ydot, q


def test_structure_is_square_and_ports_are_valid(rated):
    fluid, system, y, ydot, _ = rated
    n = system.registry.size
    assert n == 3 * (17 + 17 + 18) + 8 + 1 + 1          # 유로 3개 + 포트 8 + 우회 보정 1 + p_sys 1
    assert len(system.residual_ids()) == n and np.all(system.residual_scales() > 0)
    assert set(fluid.ports()) == {"helium_inlet", "helium_outlet"}
    assert len(fluid.cell_order) == 51 and fluid.cell_order[17] == ("hot", 0)
    assert system.diagnose_jacobian(0.0, y, ydot, cj=0.0).rank == n
    assert system.diagnose_jacobian(0.0, y, ydot, cj=1.0).rank == n


def test_rated_steady_state_heat_balance_and_bypass(rated):
    fluid, system, y, ydot, q = rated
    assert np.linalg.norm(system.residual(0.0, y, ydot) / system.residual_scales(), ord=np.inf) < 1e-8
    out = system.outputs(0.0, y)
    m = 157.1
    T_mix = out["core.fluid.T_outlet_mixed"]
    assert T_mix == pytest.approx(532.15 + q.sum() / (m * HE.cp_const), abs=1e-6)
    assert out["core.fluid.bypass_fraction"] == pytest.approx(0.11, abs=1e-10)
    reg = system.registry
    assert reg.get(y, "core.fluid.helium_inlet.m_dot") == pytest.approx(-reg.get(y, "core.fluid.helium_outlet.m_dot"))
    for label in ("average", "hot", "bypass"):
        faces = reg.get(y, f"core.fluid.{label}.m_dot_fluid_face")
        assert np.allclose(faces, faces[0], rtol=1e-10)    # 정상: 유로 안 유량 일정
    # 평균·고온 유로는 같은 Δp에서 기둥당 유량이 거의 같다 (같은 출력·유로)
    T_avg = out["core.fluid.average.T_he"][-1]
    # 보고서 [2] C §5.1: 11 % 우회가 단열이고 350 MW 전량이 연료 유로로 가면 유효 채널 출구 1014 K
    assert T_avg == pytest.approx(1014.0, abs=1.0)
    T_hot = out["core.fluid.hot.T_he"][-1]
    assert abs(T_hot - T_avg) < 30.0
    assert out["core.fluid.bypass.T_he"][-1] == pytest.approx(532.15, abs=1e-6)   # 우회 유로는 가열 없음


def test_rated_pressure_drop_is_in_benchmark_range(rated):
    fluid, system, y, ydot, _ = rated
    out = system.outputs(0.0, y)
    dp = out["core.fluid.dp_plenum_to_plenum"]
    # R3 (INL) Ex. 2 노심 압력강하 19.95 kPa. 하부 블록 수력직경을 연료 냉각공 값으로 둔 현재 geometry에서는
    # 하부 블록 마찰이 과대평가된다 (handoff 참고). 범위만 확인한다.
    assert 1.0e4 < dp < 4.0e4
    p_avg = system.registry.get(y, "core.fluid.average.p_fluid_cell")
    assert p_avg[0] < 6.3373e6 + dp and p_avg[-1] > P_OUT - 1e3


def test_helium_inventory_contract(rated):
    fluid, system, y, _, _ = rated
    assert isinstance(fluid, HeliumInventory)
    V = fluid.helium_volumes()
    assert V.shape == (51,) and np.all(V > 0)
    vol = GEO["core"]["channels"]
    assert V[:17].sum() == pytest.approx(vol["average"]["volume_m3"], rel=1e-3)
    assert V[17:34].sum() == pytest.approx(vol["hot"]["volume_m3"], rel=1e-3)
    assert V[34:].sum() == pytest.approx(vol["bypass"]["volume_m3"], rel=1e-2)
    M = total_helium_inventory([fluid], 0.0, y, system.context)
    assert 50.0 < M < 120.0          # 노심 유로 약 21.5 m³ × 3–6 kg/m³


def test_coupling_state_is_ready_for_role3(rated):
    fluid, system, y, _, _ = rated
    s = fluid.coupling_state(0.0, y, system.context)
    for key in ("T_bulk", "rho", "mu", "k", "cp", "m_dot", "D_h", "A_flow", "dz"):
        assert s[key].shape == (51,), key
    assert s["shape"][40] == "parallel_plates" and s["shape"][0] == "circular"
    assert s["m_dot"][:17].sum() / 17 == pytest.approx(system.outputs(0.0, y)["core.fluid.average.m_dot_inlet"])


def test_reverse_flow_steady_state_takes_lower_plenum_enthalpy():
    init = dict(m_dot=-5.0, T_inlet=900.0, T_outlet=700.0, p_inlet=P_OUT)
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channel_overrides=BYPASS_OVERRIDE)
    zero = np.zeros(51)
    system = build(fluid, lambda t: zero, m_in=-5.0, T_in=600.0, p_out=P_OUT, T_back=900.0)
    y, ydot = steady(system)
    reg = system.registry
    assert reg.get(y, "core.fluid.helium_inlet.m_dot") == pytest.approx(-5.0)
    # 가열 없는 역류: 하부 플레넘에서 들어온 900 K가 위로 올라가 입구 포트로 나간다
    T_out_top = float(HE.temperature(P_OUT, reg.get(y, "core.fluid.helium_inlet.h_outflow")))
    assert T_out_top == pytest.approx(900.0, abs=1e-6)
    # 역류에서 하부(출구 포트)가 높은 압력
    assert reg.get(y, "core.fluid.helium_inlet.p") < P_OUT


def test_transient_power_step_conserves_mass_and_heats_outlet():
    # 빠른 시험을 위해 고온 유로 1개만 (1기둥 유량), 우회 없음
    init = initial_from_operation(OP)
    init["m_dot"] = 2.1
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channels=("hot",))
    q0 = uniform_power(fluid, 3.5e8) / OP["reactor"]["power_distribution"]["channel_fraction"]["hot"] / 66.0
    q_of_t = lambda t: q0 * (1.0 + 0.1 * min(max(t, 0.0) / 0.5, 1.0))  # noqa: E731  0.5 s에 10 % 증가
    system = build(fluid, q_of_t, m_in=init["m_dot"], T_in=532.15, p_out=P_OUT, T_back=960.0)
    y0, _ = steady(system)
    mask = system.registry.differential_mask()
    y0, yp0 = consistent_initialize(system.residual, 0.0, y0, mask, residual_scale=system.residual_scales())
    result = DASSL(dt=0.05, rtol=1e-6, atol=1e-6, max_step=0.5,
                   residual_scale=system.residual_scales()).solve(
        system.residual, t_span=(0.0, 2.0), y0=y0, ydot0=yp0, differential_mask=mask)
    assert result.t[-1] == pytest.approx(2.0)
    reg = system.registry
    T_mix = [system.outputs(t, s)["core.fluid.T_outlet_mixed"] for t, s in zip(result.t, result.y)]
    assert T_mix[-1] > T_mix[0] + 10.0
    # 질량 보존: M(t) − M(0) = ∫ (포트 유입 합) dt  (p_sys 일정, 사다리꼴 적분)
    M = np.array([total_helium_inventory([fluid], t, s, system.context) for t, s in zip(result.t, result.y)])
    net = np.array([reg.get(s, "core.fluid.helium_inlet.m_dot") + reg.get(s, "core.fluid.helium_outlet.m_dot")
                    for s in result.y])
    integral = np.concatenate(([0.0], np.cumsum(0.5 * (net[1:] + net[:-1]) * np.diff(result.t))))
    assert np.max(np.abs((M - M[0]) - integral)) < 2e-3 * abs(M[-1] - M[0]) + 1e-6
    assert M[-1] < M[0]               # 가열로 밀도 감소 → 헬륨이 빠져나감


def test_inside_core_composite_with_minimal_solid():
    """Core 상위 조립기에 주입해 포트 이름·열률 경로가 맞는지 확인 (고체·결합은 시험 전용 최소 모델)."""
    init = initial_from_operation(OP)
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channel_overrides=BYPASS_OVERRIDE)
    q_shape = uniform_power(fluid, 1.0)            # 열률 분배 모양 (합 1)

    class LumpedSolid(Component):
        name = "test_solid"
        def register(self, registry):
            registry.register("test_solid.T", kind="differential", unit="K", initial=1000.0, scale=1000.0)
        def residual(self, t, y, ydot, context, *, heat_rate):
            return np.array([1.0e8 * context.get(ydot, "test_solid.T") - heat_rate[0]])
        def residual_ids(self):
            return ("test_solid.energy",)
        def residual_scales(self):
            return np.array([1.0e6])
        def ports(self):
            return {}
        def outputs(self, t, y, context):
            return {}

    class Power:
        def evaluate(self, t, y, context):
            return np.array([3.5e8])

    class Coupling:
        def evaluate(self, t, y, context):
            T_s = context.get(y, "test_solid.T")
            T_f = fluid.coupling_state(t, y, context)["T_bulk"][:17].mean()
            Q = 2.0e6 * (T_s - T_f)
            return Q * q_shape, np.array([Q])

    core = Core(fluid=fluid, solid=LumpedSolid(), power=Power(), coupling=Coupling())
    assert set(core.ports()) == {"fluid.helium_inlet", "fluid.helium_outlet"}
    system = SystemAssembly([core, Boundaries(fluid.name, m_in=init["m_dot"], T_in=532.15, p_out=P_OUT,
                                              T_back=960.0)])
    y, ydot = steady(system)
    out = system.outputs(0.0, y)
    assert out["core.fluid.T_outlet_mixed"] == pytest.approx(532.15 + 3.5e8 / (157.1 * HE.cp_const), abs=1e-5)
    assert system.registry.get(y, "test_solid.T") > out["core.fluid.T_outlet_mixed"]


def test_input_errors_name_the_missing_item():
    init = initial_from_operation(OP)
    with pytest.raises(ValueError, match=r"core\.channels\.bypass\.hydraulic_diameter_m"):
        ParallelChannelCoreFluid(GEO, HE, initial=init)              # 우회 D_h가 geometry에서 null
    with pytest.raises(ValueError, match="bypass_K_cal"):
        ParallelChannelCoreFluid(GEO, HE, initial=init, channel_overrides=BYPASS_OVERRIDE,
                                 mode_bypass="calibrated_parallel")
    with pytest.raises(ValueError, match="initial"):
        ParallelChannelCoreFluid(GEO, HE, channel_overrides=BYPASS_OVERRIDE)
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channel_overrides=BYPASS_OVERRIDE)
    system = build(fluid, lambda t: np.zeros(10), m_in=157.1, T_in=532.15, p_out=P_OUT, T_back=960.0)
    with pytest.raises(ValueError, match="heat_rate must have shape"):
        system.residual(0.0, system.registry.initial_state(), np.zeros(system.registry.size))
