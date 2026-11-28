"""config/operation.yaml 입력 검증 (Role 1).

기존 스키마 키의 값·단위 범위, 초기값과 정격 분포의 위치 일치,
geometry.yaml 압력손실 배분과의 교차 일관성, 보고서 §7.1 열수지 검산을 확인한다.
입력 정합성 시험이며 실제 HTGR 물리 정확도 시험이 아니다.
"""
import math
from functools import lru_cache
from numbers import Number
from pathlib import Path

import yaml

CONFIG = Path(__file__).resolve().parents[2] / "config"


@lru_cache(maxsize=1)
def _op():
    return yaml.safe_load((CONFIG / "operation.yaml").read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _geo():
    return yaml.safe_load((CONFIG / "geometry.yaml").read_text(encoding="utf-8"))


def _stations():
    return _op()["initial_distribution"]["stations"]


def _station_at(node_id):
    matches = [s for s in _stations() if s["node"] == node_id]
    assert len(matches) == 1, node_id
    return matches[0]


def test_legacy_initial_keys_are_filled_and_in_physical_range():
    op = _op()
    design_p = op["reference"]["psid_R19"]["rpv_design_pressure_Pa"]

    power = op["reactor"]["initial_thermal_power_W"]
    assert isinstance(power, Number) and power > 0

    primary = op["primary"]
    assert 1.0e5 <= primary["initial_pressure_Pa"] <= design_p
    assert 273.15 < primary["initial_temperature_K"] < 2000.0
    assert primary["initial_mass_flow_kg_s"] > 0

    comp = op["compressor"]
    omega = comp["initial_rotation_rate_rad_s"]
    assert math.isclose(omega, comp["rated_speed_rpm"] * 2 * math.pi / 60, abs_tol=0.01)
    assert 0 <= omega <= comp["speed_range_fraction"][1] * comp["rated_speed_rpm"] * 2 * math.pi / 60


def test_initial_values_match_distribution_at_named_locations():
    op, geo = _op(), _geo()
    nodes = geo["flow_network"]["nodes"]
    primary = op["primary"]

    p_node = nodes[primary["initial_pressure_location"]]["node_id"]
    assert math.isclose(primary["initial_pressure_Pa"], _station_at(p_node)["p_Pa"], abs_tol=1.0)

    t_node = nodes[primary["initial_temperature_location"]]["node_id"]
    assert math.isclose(primary["initial_temperature_K"], _station_at(t_node)["T_K"], abs_tol=1e-6)

    assert math.isclose(primary["initial_temperature_K"],
                        op["reference"]["benchmark_R1"]["core_inlet_temperature_K"], abs_tol=1e-6)


def test_distribution_follows_geometry_pressure_loss_segments():
    op, geo = _op(), _geo()
    cal = geo["pressure_loss_calibration"]
    segments = cal["segments"]
    stations = _stations()

    used = [s["loss_segment"] for s in stations[1:]]
    assert sorted(used) == sorted(segments), "모든 압력손실 구간이 정확히 한 번씩 쓰여야 한다"

    for prev, cur in zip(stations, stations[1:]):
        drop = prev["p_Pa"] - cur["p_Pa"]
        assert math.isclose(drop, segments[cur["loss_segment"]]["design_Pa"], abs_tol=60.0), cur["point"]

    rise = op["compressor"]["rated_pressure_rise_Pa"]
    assert rise == cal["circulator_pressure_rise_Pa"]
    assert math.isclose(stations[-1]["p_Pa"] + rise, stations[0]["p_Pa"], abs_tol=500.0)

    psid = op["reference"]["psid_R19"]["pressure_Pa"]
    assert math.isclose(_station_at(300)["p_Pa"], psid["sg_inlet"], abs_tol=5000.0)
    assert math.isclose(_station_at(350)["p_Pa"], psid["circulator_inlet"], abs_tol=5000.0)


def test_distribution_nodes_exist_and_densities_follow_ideal_gas():
    geo = _geo()
    ids = {n["node_id"] for n in geo["flow_network"]["nodes"].values()}
    ids |= {geo["heat_exchanger"]["helium_side"]["fsh_node_id"]}
    gas_constant = _op()["reference"]["benchmark_R1"]["helium_gas_constant_J_kgK"]
    for s in _stations():
        assert s["node"] in ids, s
        assert math.isclose(s["rho_kg_m3"], s["p_Pa"] / (gas_constant * s["T_K"]), rel_tol=2e-3), s["point"]


def test_heat_balances_from_report_section_7_1():
    op = _op()
    bench = op["reference"]["benchmark_R1"]
    psid = op["reference"]["psid_R19"]

    q_core = op["reactor"]["initial_thermal_power_W"]
    dT = bench["core_outlet_temperature_K"] - bench["core_inlet_temperature_K"]
    m_balance = q_core / (bench["helium_cp_J_kgK"] * dT)
    assert math.isclose(m_balance, op["primary"]["initial_mass_flow_kg_s"], rel_tol=5e-3)

    q_sg = q_core + op["compressor"]["rated_heat_to_helium_W"] - psid["rccs_heat_removal_W"]
    assert math.isclose(q_sg, psid["sg_net_heat_load_W"], abs_tol=2.0e5)
    assert op["boundaries"]["sg_heat_sink"]["rated_heat_removal_W"] == psid["sg_net_heat_load_W"]


def test_channel_flows_sum_to_total_and_bypass_matches_geometry():
    op, geo = _op(), _geo()
    flows = op["initial_distribution"]["channel_mass_flow_kg_s"]
    total = op["primary"]["initial_mass_flow_kg_s"]
    assert math.isclose(sum(flows.values()), total, abs_tol=0.1)

    bypass_fraction = sum(geo["core"]["channels"]["bypass"]["type_I_fraction"].values())
    assert math.isclose(flows["bypass"] / total, bypass_fraction, abs_tol=1e-3)


def test_power_distribution_fractions_sum_to_one_and_sit_on_fuel():
    """균일 출력분포 [가정]: 축방향 분율은 연료 셀에만 있고 합 1, 채널 분율은 연료 기둥 수 비례."""
    dist = _op()["reactor"]["power_distribution"]
    core = _geo()["core"]
    axial = dist["axial_fraction"]
    assert len(axial) == core["axial_cells"]
    assert math.isclose(sum(axial), 1.0, abs_tol=1e-9)
    for frac, region in zip(axial, core["cell_region"]):
        assert (frac > 0) == (region == "fuel"), (frac, region)
    if dist["type"] == "uniform":
        fuel = [f for f, r in zip(axial, core["cell_region"]) if r == "fuel"]
        assert max(fuel) == min(fuel)
        assert dist["hot_channel_peaking_factor"] == 1.0

    channel = dist["channel_fraction"]
    assert math.isclose(sum(channel.values()), 1.0, abs_tol=1e-6)
    assert channel["bypass"] == 0.0, "열은 연료에 국부 침적 (R1 §I.8.1)"
    columns = {name: core["channels"][name]["fuel_columns"] for name in ("average", "hot")}
    total = sum(columns.values())
    for name, n in columns.items():
        assert math.isclose(channel[name], n / total, abs_tol=1e-7), name


def test_sg0_water_sink_estimate_reproduces_rated_heat_removal():
    """SG-0: 일정 온도 액체 물 열침원. U 추정값 × 전열면적 × LMTD가 정격 제거열 352.2 MW와 맞는지."""
    sink = _op()["boundaries"]["sg_heat_sink"]
    assert sink["mode"] == "fixed_temperature_liquid_water" and sink["phase_change"] is False
    T_in = sink["rated_helium_inlet_temperature_K"]
    T_out = sink["rated_helium_outlet_temperature_K"]
    T_w = sink["water_temperature_K"]
    assert T_w < T_out < T_in, "물이 He 출구보다 차가워야 열이 빠진다"
    lmtd = ((T_in - T_w) - (T_out - T_w)) / math.log((T_in - T_w) / (T_out - T_w))
    area = _geo()["heat_exchanger"]["heat_transfer_area_total_m2"]
    q = sink["overall_U_estimate_W_m2K"] * area * lmtd
    assert math.isclose(q, sink["rated_heat_removal_W"], rel_tol=2e-3)
    assert sink["water_mass_kg"] > 0


def _walk(node, path=()):
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, path + (str(key),))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _walk(value, path + (str(i),))
    else:
        yield path, node


def test_units_temperatures_absolute_and_pressures_positive():
    """K는 절대온도(차이값 rise/drop/delta 제외), 나머지 단위는 양수.

    PyYAML은 3.5e8처럼 지수 부호가 없는 표기를 문자열로 읽으므로 그런 실수도 잡는다.
    """
    positive_units = ("_Pa", "_W", "_kg_s", "_rad_s", "_kg")
    difference_words = ("rise", "drop", "delta")
    checked = 0
    for path, value in _walk(_op()):
        if value is None:
            continue
        # 단위 접미사는 {inlet: .., outlet: ..}나 목록의 부모 키에 붙는다 → 가장 가까운 단위 키를 쓴다
        unit_key = next((p for p in reversed(path) if p.endswith(("_K",) + positive_units)), None)
        if unit_key is None:
            continue
        if unit_key.endswith("_K") and not any(w in unit_key for w in difference_words):
            assert isinstance(value, Number) and value > 273.15, path
        else:
            assert isinstance(value, Number) and value > 0, path
        checked += 1
    assert checked > 60
