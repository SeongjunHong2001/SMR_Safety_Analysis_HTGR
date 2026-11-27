"""S0 system_loop 실제 계산 — Role 5 builder·실행기 경로.

폐루프 정상 초기화(기준압 폐합, SG-0 U 보정, 순환기 유량 지정) → 과도 잔차 재확인 → 무외란 DASSL 유지,
그리고 보정 U로 만든 과도 계의 출력 증가 과도에서 질량 보존을 본다.

PSID 정격값은 비교용이며 판정 기준이 아니다. K 보정 전이라 압력은 PSID와 직접 비교하지 않는다 (보정은 test_calibration).
SG 용기 내부 체적은 R19 치수로 채웠고(#31), 상부 영역을 110·350에 나누는 비율만 config/assumptions/provisional_s0.yaml [가정]이다.
"""
import json

import numpy as np
import pytest

from solvers import DASSL
from solvers.initialization import consistent_initialize
from system.builder import build_system
from system.simulation import jacobian_pattern, run_scenario, transfer_state
from validation import conservation

ASSUME = ["config/assumptions/provisional_s0.yaml"]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("s0_loop")
    report = run_scenario("steady_state", "system_loop", assumptions=ASSUME, output_dir=out, log=lambda *_: None)
    return report


def test_loop_ledger_and_closure(run):
    steady, ver = run.steady, run.transient
    reg, ids = steady.built.system.registry, steady.built.system.residual_ids()
    assert reg.size == len(ids) == 481 and ver.built.system.registry.size == 480
    assert sum(n == "p_sys" for n in reg.names) == 1 and reg.spec("p_sys").owner == "primary_loop"
    assert "lower_plenum.steady_pressure_closure" in ids and "sg.secondary.U_calibration" in ids
    assert steady.residual_norm < 1e-6 and ver.residual_norm < 1e-6     # 정상해가 교체 전 과도 잔차도 만족
    assert reg.get(steady.y, "p_sys") == pytest.approx(6.3373e6)


def test_loop_heat_balance_and_sg_calibration(run):
    s, bal = run.summary, run.energy_balance
    assert s["sg_heat_removal_W"] == pytest.approx(3.5e8 + s["circulator_heat_W"] - s["rccs_heat_removal_W"],
                                                   rel=1e-8)                    # 정상 폐루프 열수지
    assert bal["relative_imbalance"] < 1e-8
    assert s["sg_heat_removal_W"] == pytest.approx(3.522e8, abs=1.0e6)          # PSID SG 순 열부하 352.2 MW
    assert s["sg_helium_outlet_temperature_K"] == pytest.approx(528.15, abs=1e-6)   # U 보정 목표 255 °C
    assert 350.0 < s["sg_overall_U_W_m2K"] < 600.0
    assert s["core_outlet_temperature_K"] == pytest.approx(960.15, abs=3.0)
    assert run.checks["acceptance"]["rpv_temperature_max_K"]["pass"]
    assert 0.0 < s["circulator_pressure_rise_Pa"] < 9.1e4                       # K 미보정: PSID 91 kPa보다 작음
    # He 재고량: 유로 체적 + 정체 공간(배럴–RPV 간극, UPTPS 위, 저온 He) ≈ R22 2,585 kg (RV + SG 용기 + cross duct)
    geo = run.steady.built.inputs.geometry
    he = run.steady.built.parts["he"]
    V_stag = sum(v["volume_m3"] for v in geo["flow_network"]["stagnant_volumes"].values())
    m_stag = V_stag * float(he.density(s["p_sys_Pa"], he.enthalpy(s["p_sys_Pa"], 532.15)))
    assert s["helium_inventory_kg"] + m_stag == pytest.approx(2585.0, rel=0.03)
    assert set(run.checks["comparison"]) >= {"circulator_pressure_rise_Pa", "sg_heat_removal_W", "helium_inventory_kg"}


def test_loop_hold_transient_and_metadata(run):
    hold = run.hold
    assert hold.result.t[-1] == pytest.approx(10.0)
    assert hold.drift_scaled_max < 1e-6
    assert hold.balance["mass_error_rel_max"] < 1e-10 and hold.balance["energy_error_rel_max"] < 1e-8
    meta = json.loads(run.paths["meta"].read_text(encoding="utf-8"))
    assert meta["p_sys_owner"] == "lower_plenum" and meta["replaced_residual_id"] == "lower_plenum.mass"
    assert len(meta["assumptions"]) == 2 and meta["assumptions_outside_domain"] == []
    assert meta["solver"]["verification"]["U_sg_fixed"] == pytest.approx(run.summary["sg_overall_U_W_m2K"])


def test_loop_power_step_conserves_mass(run):
    ver = run.transient
    power = lambda t: 3.5e8 * (1 + 0.1 * min(max(t, 0.0) / 0.5, 1.0))  # noqa: E731
    built = build_system(ver.built.inputs, mode="transient", U_sg=ver.U_sg, power=power)
    system, reg = built.system, built.system.registry
    y0 = transfer_state(ver.built, ver.y, built)
    mask = reg.differential_mask()
    pattern = jacobian_pattern(system, 0.0, y0, np.zeros_like(y0))
    z_order = np.concatenate([np.flatnonzero(~mask), np.flatnonzero(mask)])
    y0, yp0 = consistent_initialize(system.residual, 0.0, y0, mask, residual_scale=system.residual_scales(),
                                    jacobian_sparsity=pattern[:, z_order])
    result = DASSL(dt=0.05, rtol=1e-5, atol=1e-5, max_step=1.0, residual_scale=system.residual_scales()).solve(
        system.residual, t_span=(0.0, 3.0), y0=y0, ydot0=yp0, differential_mask=mask, jacobian_sparsity=pattern)
    assert result.t[-1] == pytest.approx(3.0)
    bal = conservation.transient_balance(built, result.t, result.y)
    assert bal["mass_error_rel_max"] < 1e-8                                    # 닫힌 루프 질량 보존
    assert bal["energy_error_rel_max"] < 2e-3                                  # 램프 꺾임의 사다리꼴 오차 (≈0.8e-3)
    assert bal["energy_error_rel_end"] < 2e-4                                  # 끝값 ≈ 0.08e-3
    assert reg.get(result.y[-1], "p_sys") > reg.get(result.y[0], "p_sys")      # 가열 → 계통압 상승
    out0, out1 = system.outputs(0.0, result.y[0]), system.outputs(3.0, result.y[-1])
    assert out1["core.solid.T_fuel_max"] > out0["core.solid.T_fuel_max"]
