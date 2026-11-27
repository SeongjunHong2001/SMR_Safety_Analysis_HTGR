"""S0 core_benchmark 첫 실제 계산 (R1 Ex. 2b 경계, 균일 출력분포 [가정]).

``tests/integration/s0_core_benchmark.py``의 조립 원형으로 정상해를 풀고 R1·R3 참고값, 열수지, 판정 기준과 비교한다.
비교값은 판정 기준이 아니며 R3는 다른 출력분포·경계를 썼다 (steady_state.yaml comparison 주석).
"""
import numpy as np
import pytest

from solvers import solve_steady_state
from tests.integration.s0_core_benchmark import build


@pytest.fixture(scope="module")
def s0():
    system, parts = build()
    y, ydot = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                                 residual_scale=system.residual_scales())
    rates = system._trial_context(0.0, y).heat_rates
    return system, parts, y, ydot, system.outputs(0.0, y), rates


def test_s0_is_square_and_converged(s0):
    system, parts, y, ydot, out, rates = s0
    assert system.registry.size == len(system.residual_ids())
    assert np.linalg.norm(system.residual(0.0, y, ydot) / system.residual_scales(), ord=np.inf) < 1e-6
    assert system.registry.get(y, "p_sys") == pytest.approx(6.39e6)
    assert out["core.fluid.bypass_fraction"] == pytest.approx(0.11, abs=1e-9)


def test_s0_heat_balance_and_links_conserve(s0):
    system, parts, y, ydot, out, rates = s0
    for name in ("barrel_riser", "riser_rpv"):
        assert rates[f"{name}.to_a"].sum() == pytest.approx(-rates[f"{name}.to_b"].sum(), rel=1e-10)
    q_rccs = -rates["rpv_rccs.to_a"].sum()
    he = parts["he"]
    h_in = float(he.enthalpy(6.39e6, 532.15))
    h_out = float(system.registry.get(y, parts["lower"].port("from_core").h_instream))
    # 정상: 출력 = 헬륨 엔탈피 증가 + RCCS 손실 (상·하부 헤드 단열)
    assert 157.1 * (h_out - h_in) == pytest.approx(3.5e8 - q_rccs, rel=1e-7)
    assert 1.0e5 < q_rccs < 2.0e6


def test_s0_results_against_benchmark_references(s0):
    system, parts, y, ydot, out, rates = s0
    # R1 Table I.8 결과: 노심 출구 ≈ 687 °C (960.15 K)
    assert out["core.fluid.T_outlet_mixed"] == pytest.approx(960.15, abs=2.0)
    # 판정 기준: RPV 최고 370 °C (643.15 K) 미만 (R9, 보고서 D9)
    assert out["rpv.T_max"] < 643.15
    # R3 (INL, 다른 출력분포·경계): 배럴 평균 595 K, 노심 Δp 19.95 kPa, 연료 최고 1261 K
    assert out["core.solid.core_barrel.T"].mean() == pytest.approx(595.15, abs=60.0)
    assert 1.5e4 < out["core.fluid.dp_plenum_to_plenum"] < 3.5e4
    assert 1050.0 < out["core.solid.T_fuel_max"] < 1300.0
    # riser는 배럴에서 열을 받아 노심 입구 온도가 입구 경계보다 약간 높다
    assert 532.15 < out["riser.T_he"][-1] < 545.0
