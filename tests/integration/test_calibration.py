"""Role 5 압력손실 K 보정 시험: 연결부 K_form → 면 적용 규칙, 보정 구간, 보정 파일 규칙, PSID 설계손실 도달.

보정 자체는 S0 system_loop 정상해 3번(약 20 s). 저장된 config/calibration/s0_pressure_loss_K.yaml이 현재 입력으로
다시 보정한 값과 같은지도 본다 — 입력이 바뀌면 이 시험이 실패하고 파일을 다시 만들어야 한다 (python main.py --calibrate).
"""
import dataclasses

import pytest
import yaml

from system.builder import ROOT, build_system, load_inputs
from system.calibration import calibrate_pressure_losses, segments, with_junction_k
from system.simulation import run_hold, settings_from_scenario, verify_steady

CAL_FILE = "config/calibration/s0_pressure_loss_K.yaml"


@pytest.fixture(scope="module")
def loop_inputs():
    return load_inputs("steady_state", "system_loop")


def test_junction_k_goes_to_the_adjacent_flow_face(loop_inputs):
    K = {"J12": 4.0, "J06": 10.0, "J16": 100.0, "J02": -0.1, "J05": 1.5}
    built = build_system(with_junction_k(loop_inputs, K))
    p = built.parts
    assert p["hot_duct"].inlet == ("plenum", 0.5 + 4.0, 1.0 + 4.0)          # J12: 230 → 240 입구 면
    assert p["cold_annulus"].inlet == pytest.approx(("plenum", 0.4, 0.9))   # J02 음의 증분 (합은 양수)
    assert p["riser"].outlet == ("plenum", 1.0 + 1.5, 0.5 + 1.5)            # J05: riser → 150 출구 면
    assert p["sg_he"].pipe.outlet == ("plenum", 101.0, 100.5)               # J16: SG → 330
    assert p["fluid"].face_K["average"][0] == 10.0 and "hot" not in p["fluid"].face_K
    assert built.k_form == K and built.describe()["junction_K_form"] == K
    with pytest.raises(ValueError, match="no pipe/core/SG face"):
        build_system(with_junction_k(loop_inputs, {"J19": 1.0}))          # MLSV → 순환기
    with pytest.raises(NotImplementedError, match="internal"):
        build_system(with_junction_k(loop_inputs, {"J15": 1.0}))          # SG 내부 FSH → EES
    with pytest.raises(ValueError, match="negative"):
        build_system(with_junction_k(loop_inputs, {"J02": -0.6}))         # 입구 급축소 0.5보다 큰 감소


def test_segments_carry_k_on_the_first_placeable_junction(loop_inputs):
    segs = {s.name: s.k_junctions for s in segments(loop_inputs)}
    assert segs == {"mc_outlet_to_core_inlet": ("J02",), "riser": ("J05",), "core": ("J06", "J07"),
                    "core_outlet_plenum": ("J12",), "hot_duct": ("J13",), "sg_inlet_plenum": ("J14",),
                    "sg_bundle": ("J16",), "sg_outlet": ("J17",), "mc_inlet": ("J18",)}
    # 노심 묶음은 연료 유로(J06 평균, J07 고온)만. Type I 고정 분율 우회(J08)는 분율식이 유량을 정한다


def test_calibration_file_rules(tmp_path):
    with pytest.raises(ValueError, match="calibrated for condition 'system_loop'"):
        load_inputs("steady_state", "core_benchmark", assumptions=[CAL_FILE])
    dep = tmp_path / "cal_dep.yaml"                       # 가정 파일에 기대는 보정은 그 파일도 넘겨야 한다
    dep.write_text(yaml.safe_dump({"meta": {"kind": "calibration", "basis": {"condition": "system_loop"},
                                            "depends_on_assumptions": ["config/assumptions/x.yaml"]},
                                   "geometry": {"junctions": {"J05": {"K_form": 1.0}}}}), encoding="utf-8")
    with pytest.raises(ValueError, match="load those assumption files too"):
        load_inputs("steady_state", "system_loop", assumptions=[dep])
    inputs = load_inputs("steady_state", "system_loop", assumptions=[CAL_FILE])
    assert inputs.files["calibration[0]"] == CAL_FILE
    d = build_system(inputs).describe()
    assert len(d["calibration"]) == 10 and d["assumptions"] == []


@pytest.fixture(scope="module")
def cal(loop_inputs):
    return calibrate_pressure_losses(loop_inputs, log=lambda *_: None)


def test_calibration_reaches_the_psid_allocation(cal):
    for name, target in cal.targets.items():
        assert cal.losses[name]["irreversible"] == pytest.approx(target, abs=1.0)
    assert cal.losses["unassigned"]["irreversible"] == 0.0
    s = cal.summary
    assert s["sum_irreversible_Pa"] == pytest.approx(90800.0, abs=10.0)       # PSID design 합
    assert s["circulator_pressure_rise_Pa"] == pytest.approx(91000.0, abs=1500.0)   # PSID 91 kPa (+가역 항)
    assert 0.0 < s["reversible_net_Pa"] < 2000.0                                # 루프 정수압 + 열가속
    assert cal.iterations <= 4
    # 지점 압력: 노드 위치 정의 차이(정수압·입구 동압) 수준에서 PSID 분포와 맞다
    diffs = {r["station"]: r["difference_Pa"] for r in cal.stations}
    assert diffs[7] == 0.0 and max(abs(v) for v in diffs.values()) < 3.5e3


def test_saved_calibration_file_is_current(cal):
    saved = yaml.safe_load((ROOT / CAL_FILE).read_text(encoding="utf-8"))
    K_saved = {j: v["K_form"] for j, v in saved["geometry"]["junctions"].items()}
    assert set(K_saved) == set(cal.K)
    for j, v in cal.K.items():
        assert K_saved[j] == pytest.approx(v, rel=1e-4, abs=1e-5), f"{j}: rerun python main.py --calibrate"
    assert saved["meta"]["depends_on_assumptions"] == []


def test_calibrated_loop_holds_steady(cal):
    settings = settings_from_scenario(cal.built.inputs)
    ver = verify_steady(cal.steady, settings)
    assert ver.residual_norm < 1e-6
    hold = run_hold(ver, dataclasses.replace(settings.hold, t_end_s=5.0), settings)
    assert hold.drift_scaled_max < 1e-6
    assert hold.balance["mass_error_rel_max"] < 1e-10 and hold.balance["energy_error_rel_max"] < 1e-8
