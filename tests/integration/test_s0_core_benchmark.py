"""S0 core_benchmark 실제 계산 (R1 Ex. 2b 경계, 균일 출력분포 [가정]) — Role 5 builder·실행기 경로.

``system.simulation.run_scenario``로 조립 → 정상해 → 과도 잔차 재확인 → 무외란 DASSL 유지 → 출력까지 돌리고
R1·R3 참고값, 열수지, 판정 기준과 비교한다. 비교값은 판정 기준이 아니며 R3는 다른 출력분포·경계를 썼다
(steady_state.yaml comparison 주석). 가정 파일 없이 config의 원문·계산 입력만으로 돈다 (#31: 우회 D_h 4A/P_w,
하부 블록 셀별 D_h를 R1·PSID 치수로 채움).
"""
import csv
import json

import numpy as np
import pytest

from system.simulation import flatten_names, run_scenario

@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("s0_core")
    report = run_scenario("steady_state", "core_benchmark", output_dir=out, log=lambda *_: None)
    built, y = report.transient.built, report.transient.y
    system = built.system
    return report, built, y, system.outputs(0.0, y), system._trial_context(0.0, y).heat_rates


def test_s0_is_square_converged_and_verified(run):
    report, built, y, out, rates = run
    system = built.system
    assert system.registry.size == len(system.residual_ids()) == 347
    assert report.steady.residual_norm < 1e-6 and report.transient.residual_norm < 1e-6
    assert system.registry.get(y, "p_sys") == pytest.approx(6.39e6)
    assert out["core.fluid.bypass_fraction"] == pytest.approx(0.11, abs=1e-9)


def test_s0_heat_balance_and_links_conserve(run):
    report, built, y, out, rates = run
    for name in ("barrel_riser", "riser_rpv"):
        assert rates[f"{name}.to_a"].sum() == pytest.approx(-rates[f"{name}.to_b"].sum(), rel=1e-10)
    q_rccs = -rates["rpv_rccs.to_a"].sum()
    he = built.parts["he"]
    h_in = float(he.enthalpy(6.39e6, 532.15))
    h_out = float(built.system.registry.get(y, built.parts["lower_plenum"].port("inlet").h_instream))
    # 정상: 출력 = 헬륨 엔탈피 증가 + RCCS 손실 (상·하부 헤드 단열)
    assert 157.1 * (h_out - h_in) == pytest.approx(3.5e8 - q_rccs, rel=1e-7)
    assert 1.0e5 < q_rccs < 2.0e6
    assert report.energy_balance["relative_imbalance"] < 1e-9


def test_s0_results_against_benchmark_references(run):
    report, built, y, out, rates = run
    s, checks = report.summary, report.checks
    # R1 Table I.8 결과: 노심 출구 ≈ 687 °C (960.15 K)
    assert s["core_outlet_temperature_K"] == pytest.approx(960.15, abs=2.0)
    assert checks["comparison"]["core_outlet_temperature_K"]["reference"] == 960.15
    # 판정 기준: RPV 최고 370 °C (643.15 K) 미만 (R9, 보고서 D9)
    assert checks["acceptance"]["rpv_temperature_max_K"]["pass"]
    # R3 (INL, 다른 출력분포·경계): 배럴 평균 595 K, 노심 Δp 19.95 kPa, 연료 최고 1261 K
    assert s["core_barrel_temperature_avg_K"] == pytest.approx(595.15, abs=60.0)
    assert s["core_pressure_drop_Pa"] == pytest.approx(19950.0, rel=0.10)   # 하부 블록 D_h 반영 후 20.7 kPa
    assert 1050.0 < s["fuel_temperature_max_K"] < 1300.0
    assert sum(v for k, v in s.items() if k.startswith("channel_mass_flow")) == pytest.approx(157.1, rel=1e-9)
    # riser는 배럴에서 열을 받아 노심 입구 온도가 입구 경계보다 약간 높다
    assert 532.15 < out["riser.T_he"][-1] < 545.0


def test_s0_hold_transient_and_outputs(run):
    report, built, y, out, rates = run
    hold = report.hold
    assert hold.result.t[-1] == pytest.approx(10.0) and hold.result.rejected_steps == 0
    assert hold.ydot0_differential_max < 1e-8 and hold.drift_scaled_max < 1e-8
    assert hold.balance["mass_error_rel_max"] < 1e-10 and hold.balance["energy_error_rel_max"] < 1e-8
    meta = json.loads(report.paths["meta"].read_text(encoding="utf-8"))
    assert meta["scenario"] == "steady_state" and meta["condition"] == "core_benchmark"
    assert meta["solver"]["hold_transient"]["accepted_steps"] == len(hold.result.orders)
    assert meta["assumptions"] == [] and meta["calibration"] == []
    assert meta["units"]["p_sys"] == "Pa" and meta["code_version"]["commit"]
    with report.paths["states"].open(encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["time_s", *flatten_names(built.system.registry)[0]]
    assert len(rows) == len(hold.result.t) + 1
    assert np.allclose([float(v) for v in rows[1][1:]], hold.result.y[0])
