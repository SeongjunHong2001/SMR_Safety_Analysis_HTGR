"""Role 3 노심 고체·출력·결합 시험과 Core 상위 조립 (유체 + 고체 + 출력 + 결합) 시험.

실제 geometry.yaml·operation.yaml(균일 출력분포)·R1 물성으로 조립한다. 플레넘 대신 시험 전용 경계
(tests/components/test_core_fluid.py의 Boundaries)를 쓴다. 실제 HTGR 정상상태 검증 결과가 아니다.
"""
import math

import numpy as np
import pytest

from interfaces.contracts import SimulationContext
from models._inputs import ModelInputMissing
from models.core import Core
from models.core.coupling import UnitCellCoreCoupling
from models.core.fluid import ParallelChannelCoreFluid, initial_from_operation
from models.core.power import SpecifiedCorePower
from models.core.solid import UnitCellCoreSolid
from physics.heat_transfer.correlations import make_convection_correlation
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN
from physics.heat_transfer.solid_properties import R1SolidProperties
from solvers import DASSL, solve_steady_state
from solvers.initialization import consistent_initialize
from solvers.jacobian import finite_difference_jacobian
from system.assembly import SystemAssembly
from tests.components.test_core_fluid import BYPASS_OVERRIDE, GEO, HE, OP, P_OUT, Boundaries

SP = R1SolidProperties()
DIST = OP["reactor"]["power_distribution"]


def build_core(*, channels=("average", "hot", "bypass"), rings=True, power=3.5e8, m_dot=None):
    init = initial_from_operation(OP)
    if m_dot is not None:
        init["m_dot"] = m_dot
    overrides = BYPASS_OVERRIDE if "bypass" in channels else None
    fluid = ParallelChannelCoreFluid(GEO, HE, initial=init, channels=channels, channel_overrides=overrides)
    heated = tuple(c for c in channels if c != "bypass")
    solid = UnitCellCoreSolid(GEO, SP, helium=HE, channels=heated, rings=rings,
                              initial=dict(T=np.linspace(init["T_inlet"], init["T_outlet"], 17)))
    dist = DIST
    if heated == ("hot",):                                     # 고온 유로 단독 시험: 분율을 그 유로에 몰아줌
        dist = dict(DIST, channel_fraction={"hot": 1.0})
    pw = SpecifiedCorePower(solid, dist, power)
    coupling = UnitCellCoreCoupling(fluid, solid, make_convection_correlation("mceligot"),
                                    bypass_regions=GEO["core"]["channels"]["bypass"]["type_I_fraction"])
    core = Core(fluid=fluid, solid=solid, power=pw, coupling=coupling)
    system = SystemAssembly([core, Boundaries(fluid.name, m_in=init["m_dot"], T_in=init["T_inlet"], p_out=P_OUT,
                                              T_back=init["T_outlet"])])
    return system, dict(fluid=fluid, solid=solid, power=pw, coupling=coupling, init=init)


@pytest.fixture(scope="module")
def rated():
    system, parts = build_core()
    y, ydot = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                                 residual_scale=system.residual_scales())
    return system, parts, y, ydot


# ------------------------------------------------------------------ 물성·출력·형상
def test_solid_properties_match_r1_and_report3():
    assert SP.conductivity("fuel_compact") == 20.0 and SP.conductivity("fuel_block") == 37.0
    assert SP.conductivity("graphite_H451") == 66.0 and SP.conductivity("graphite_2020") == 35.0
    assert SP.conductivity("alloy_800H") == 17.8
    assert SP.cp("graphite_H451", 1000.0) == pytest.approx(1759.0, abs=2.0)      # 보고서 [3]: 1000 K ≈ 1759
    assert SP.cp("alloy_800H", 700.0) == 515.0
    assert SP.density("graphite_H451") == 1850.0 and SP.emissivity("graphite_2020") == 0.85
    with pytest.raises(ValueError, match="ex2b_fixed"):
        R1SolidProperties(mode="ex2c_fluence")
    with pytest.raises(ValueError, match="unknown solid material"):
        SP.conductivity("concrete")


def test_unit_cell_geometry_and_conductances():
    solid = UnitCellCoreSolid(GEO, SP, helium=HE, initial=dict(T=800.0))
    avg, hot = solid.column("average"), solid.column("hot")
    assert (avg.N_f, avg.N_c) == (53 * 210 + 12 * 186, 53 * 108 + 12 * 95)      # 층당 연료공·냉각공
    assert (hot.N_f, hot.N_c) == (210, 108)
    assert avg.n_fuel == 10 and solid.grid_size == 2 * (10 + 17) + 4 * 17
    # 두 평행 원통 형상계수: p 18.796 mm, r_f 6.35 mm, r_c 7.94 mm
    arg = (0.018796 ** 2 - 0.00635 ** 2 - 0.00794 ** 2) / (2 * 0.00635 * 0.00794)
    assert solid.S_pair == pytest.approx(2 * math.pi / math.acosh(arg), rel=1e-12)
    # 컴팩트–블록 컨덕턴스를 손으로 계산한 값과 비교 (800 K, 6.39 MPa)
    T, p = 800.0, 6.39e6
    h_gap = float(HE.conductivity(p, T)) / 1e-4 + 4 * STEFAN_BOLTZMANN * solid.eps_gap * T ** 3
    ratio = hot.N_c / hot.N_f
    R = ratio / (8 * math.pi * 20.0) + ratio / (2 * math.pi * 0.006225 * h_gap) + 0.5 / (6 * 66.0 * solid.S_pair)
    g = solid.fuel_block_conductance(hot, np.array([T]), np.array([T]), p)
    assert g[0] == pytest.approx(hot.N_c * hot.dz[hot.fuel][0] / R, rel=1e-12)
    # 연료 영역 흑연 + 컴팩트 체적 ≈ 기둥 체적 − 냉각공 (66기둥 × 육각 0.112237 m²)
    fuel_cells = avg.fuel
    total = avg.V_block[fuel_cells] + avg.V_fuel + avg.N_f * math.pi * (0.00635 ** 2 - 0.006225 ** 2) * avg.dz[fuel_cells]
    assert np.allclose(total, avg.A_axial * avg.dz[fuel_cells], rtol=1e-12)


def test_internal_conduction_conserves_energy():
    solid = UnitCellCoreSolid(GEO, SP, helium=HE, initial=dict(T=700.0))
    T = np.random.default_rng(1).uniform(550.0, 1100.0, solid.grid_size)
    assert abs(solid.conduction_net(T, 6.39e6).sum()) < 1e-6 * np.abs(solid.conduction_net(T, 6.39e6)).max()
    assert np.allclose(solid.conduction_net(np.full(solid.grid_size, 900.0), 6.39e6), 0.0)


def test_specified_power_distribution():
    solid = UnitCellCoreSolid(GEO, SP, helium=HE, initial=dict(T=800.0))
    pw = SpecifiedCorePower(solid, DIST, lambda t: 3.5e8 * (1 + t))
    q = pw.evaluate(1.0, None, None)
    assert q.sum() == pytest.approx(7.0e8)
    assert np.all(q[solid.slices[("average", "block")]] == 0) and np.all(q[solid.slices[("core_barrel",)]] == 0)
    assert q[solid.slices[("hot", "fuel")]].sum() == pytest.approx(7.0e8 * DIST["channel_fraction"]["hot"])
    with pytest.raises(ModelInputMissing, match="D7"):
        SpecifiedCorePower(solid, DIST, 3.5e8, mode="fission_plus_decay")
    with pytest.raises(ValueError, match="no fuel"):
        SpecifiedCorePower(solid, dict(DIST, channel_fraction={"average": 0.5, "hot": 0.0, "bypass": 0.5}), 3.5e8)


def test_outer_link_heat_rate_goes_to_barrel():
    solid = UnitCellCoreSolid(GEO, SP, helium=HE, initial=dict(T=800.0), outer_link="core_rpv")
    system, _ = build_core()
    ctx = SimulationContext(system.registry, heat_rates={"core_rpv.to_a": np.full(17, -1.0e3)})
    q = solid.external_rate(ctx)
    assert q.sum() == pytest.approx(-1.7e4) and np.all(q[solid.slices[("core_barrel",)]] == -1.0e3)
    with pytest.raises(ModelInputMissing, match="core_rpv.to_a"):
        solid.external_rate(SimulationContext(system.registry))
    assert set(solid.ports()) == {f"outer_surface[{i}]" for i in range(17)}


# ------------------------------------------------------------------ Core 정격 정상해
def test_core_composite_is_square_and_regular(rated):
    system, parts, y, ydot = rated
    assert system.registry.size == 166 + 122                       # 유체 166 (p_sys 포함) + 고체 격자 122
    assert system.diagnose_jacobian(0.0, y, ydot, cj=1.0).rank == system.registry.size


def test_core_rated_steady_state(rated):
    system, parts, y, ydot = rated
    assert np.linalg.norm(system.residual(0.0, y, ydot) / system.residual_scales(), ord=np.inf) < 1e-6
    out = system.outputs(0.0, y)
    q_fluid, q_solid = parts["coupling"].evaluate(0.0, y, system.context)
    assert q_fluid.sum() == pytest.approx(3.5e8, rel=1e-9)        # 외면 단열: 출력 전부 헬륨으로
    assert q_fluid.sum() == pytest.approx(q_solid.sum(), rel=1e-12)
    assert out["core.fluid.T_outlet_mixed"] == pytest.approx(532.15 + 3.5e8 / (157.1 * HE.cp_const), abs=1e-5)
    # 보고서 [2] C §5.1: 우회가 39.1 MW를 흡수하면 유효 채널 출구 687 °C (960 K). 이 모델: 반사체 전도로 우회가 가열됨.
    # 우회 D_h = 4A/P_w(#31)로 젖은 면적이 4 mm 가정의 절반이라 우회 흡수 ≈ 34.5 MW, 평균 유로 출구 ≈ 966 K
    q_bypass = q_fluid[34:].sum()
    assert 3.0e7 < q_bypass < 4.5e7
    for ch, sl in (("bypass", slice(34, None)), ("average", slice(0, 17))):
        m_ch = out[f"core.fluid.{ch}.m_dot_inlet"]                 # 유로별 에너지 수지: ṁ c_p ΔT = 흡수 열
        assert m_ch * HE.cp_const * (out[f"core.fluid.{ch}.T_he"][-1] - 532.15) == pytest.approx(q_fluid[sl].sum(),
                                                                                                rel=1e-6)
    assert 955.0 < out["core.fluid.average.T_he"][-1] < 975.0
    # 연료 온도: 균일 출력이라 R3(Ex. 2, 실제 분포) 최고 1261 K보다 낮다
    T_max = out["core.solid.T_fuel_max"]
    assert 1050.0 < T_max < 1300.0
    centers = out["core.solid.average.T_compact_center"]
    assert np.all(np.diff(centers) > 0)                            # 하향 유동: 아래로 갈수록 뜨거움
    assert np.all(out["core.solid.average.T_fuel"] > out["core.solid.average.T_block"][parts["solid"].column("average").fuel])


# ------------------------------------------------------------------ 과도
def test_hot_column_transient_energy_balance():
    system, parts = build_core(channels=("hot",), rings=False, power=3.5e8 / 66, m_dot=157.1 / 66 * 0.89)
    y_s, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                                residual_scale=system.residual_scales())
    pw = parts["power"]
    pw.total_power = lambda t: 3.5e8 / 66 * (1.0 + 0.2 * min(max(t, 0.0) / 0.2, 1.0))
    mask = system.registry.differential_mask()
    y0, yp0 = consistent_initialize(system.residual, 0.0, y_s, mask, residual_scale=system.residual_scales())
    Jy = finite_difference_jacobian(lambda z: system.residual(0.0, z, yp0), y0)
    Jd = finite_difference_jacobian(lambda z: system.residual(0.0, y0, z), yp0)
    pattern = (Jy != 0) | (Jd != 0) | np.eye(y0.size, dtype=bool)
    result = DASSL(dt=0.05, rtol=1e-6, atol=1e-6, max_step=1.0, residual_scale=system.residual_scales()).solve(
        system.residual, t_span=(0.0, 5.0), y0=y0, ydot0=yp0, differential_mask=mask, jacobian_sparsity=pattern)
    assert result.t[-1] == pytest.approx(5.0)
    solid, coupling = parts["solid"], parts["coupling"]
    for t, s, sd in zip(result.t[::5], result.y[::5], result.ydot[::5]):
        T = solid.temperatures(s, system.context)
        stored = float(np.sum(solid.heat_capacity(T) * solid._ydot(sd, system.context)))
        q_fluid, _ = coupling.evaluate(t, s, system.context)
        assert stored == pytest.approx(pw.total(t) - q_fluid.sum(), rel=1e-6, abs=1e-3 * pw.total(t))
    T_first = system.outputs(0.0, result.y[0])["core.solid.T_fuel_max"]
    T_last = system.outputs(5.0, result.y[-1])["core.solid.T_fuel_max"]
    assert T_last > T_first + 5.0                                  # 출력 20 % 증가 → 연료 온도 상승
