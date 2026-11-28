"""config/scenarios/steady_state.yaml (S0) 입력 검증 (Role 1).

조건 세트 구조, 계산 영역과 연결부의 정합성(경계 누수 없음, 폐루프 닫힘),
operation.yaml·geometry.yaml과의 값 일치, 단위 범위를 확인한다.
입력 정합성 시험이며 실제 HTGR 정상상태 해석 결과 시험이 아니다.
"""
import math
from functools import lru_cache
from numbers import Number
from pathlib import Path

import yaml

CONFIG = Path(__file__).resolve().parents[2] / "config"


def _load(relative):
    return yaml.safe_load((CONFIG / relative).read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _s0():
    return _load("scenarios/steady_state.yaml")


@lru_cache(maxsize=1)
def _op():
    return _load("operation.yaml")


@lru_cache(maxsize=1)
def _geo():
    return _load("geometry.yaml")


def _node_name_to_id():
    return {name: n["node_id"] for name, n in _geo()["flow_network"]["nodes"].items()}


def _active_junctions(condition):
    closed = set(condition["domain"].get("closed_junctions", []))
    return {jid: j for jid, j in _geo()["junctions"].items() if jid not in closed}


def test_structure_follows_schema():
    s0 = _s0()
    assert s0["name"] == "steady_state"
    conditions = s0["conditions"]
    assert isinstance(conditions, dict) and conditions
    assert s0["default_condition"] in conditions
    assert set(conditions) == {"core_benchmark", "system_loop"}
    for name, cond in conditions.items():
        assert cond["reactor"]["thermal_power_W"] > 0, name
        assert cond["outputs"], name


def test_domain_references_exist_in_geometry():
    geo = _geo()
    junctions = geo["junctions"]
    structures = geo["heat_structures"]
    node_ids = {j["from"] for j in junctions.values()} | {j["to"] for j in junctions.values()}
    boundary_ids = {b["node_id"] for b in geo["boundaries"].values()}
    for name, cond in _s0()["conditions"].items():
        domain = cond["domain"]
        nodes = set(domain["nodes"])
        assert nodes <= node_ids, (name, nodes - node_ids)
        for jid in domain.get("closed_junctions", []):
            assert jid in junctions, (name, jid)
        for hs in domain["heat_structures"]:
            assert hs in structures, (name, hs)
            # 열구조체가 닿는 수력 노드는 계산 영역 안이거나 경계 노드(930 등)여야 한다
            for side in ("inner", "outer"):
                refs = structures[hs][side]
                refs = refs if isinstance(refs, list) else [refs]
                for ref in refs:
                    if isinstance(ref, int):
                        assert ref in nodes | boundary_ids, (name, hs, side, ref)


def test_core_benchmark_has_one_inlet_one_outlet_and_no_leaks():
    cond = _s0()["conditions"]["core_benchmark"]
    nodes = set(cond["domain"]["nodes"])
    active = _active_junctions(cond)
    inlet = cond["inlet_boundary"]["junction"]
    outlet = cond["outlet_boundary"]["junction"]

    assert active[inlet]["from"] not in nodes and active[inlet]["to"] in nodes
    assert active[outlet]["from"] in nodes and active[outlet]["to"] not in nodes

    crossing = {jid for jid, j in active.items() if (j["from"] in nodes) != (j["to"] in nodes)}
    assert crossing == {inlet, outlet}, crossing

    outlet_node = _node_name_to_id()[cond["outlet_boundary"]["pressure_location"]]
    assert active[outlet]["from"] == outlet_node


def test_system_loop_is_closed_and_every_node_has_in_and_out_flow():
    cond = _s0()["conditions"]["system_loop"]
    nodes = set(cond["domain"]["nodes"])
    active = _active_junctions(cond)

    crossing = {jid for jid, j in active.items() if (j["from"] in nodes) != (j["to"] in nodes)}
    assert not crossing, f"닫힌 루프인데 계산 영역 밖으로 열린 연결부: {crossing}"

    internal = [j for j in active.values() if j["from"] in nodes and j["to"] in nodes]
    for node in nodes:
        assert any(j["to"] == node for j in internal), f"{node}로 들어오는 연결부 없음"
        assert any(j["from"] == node for j in internal), f"{node}에서 나가는 연결부 없음"

    hydraulic = ({n["node_id"] for n in _geo()["flow_network"]["nodes"].values()}
                 | {c["node_id"] for c in _geo()["core"]["channels"].values()}
                 | {_geo()["heat_exchanger"]["helium_side"]["fsh_node_id"],
                    _geo()["heat_exchanger"]["helium_side"]["ees_node_id"]})
    assert hydraulic <= nodes, f"1차 루프 노드 누락: {hydraulic - nodes}"


def test_core_benchmark_values_match_operation_reference():
    cond = _s0()["conditions"]["core_benchmark"]
    op = _op()
    bench = op["reference"]["benchmark_R1"]
    assert cond["reactor"]["thermal_power_W"] == op["reactor"]["initial_thermal_power_W"] == bench["thermal_power_W"]
    assert cond["inlet_boundary"]["mass_flow_kg_s"] == bench["mass_flow_kg_s"]
    assert cond["inlet_boundary"]["temperature_K"] == bench["core_inlet_temperature_K"]
    assert cond["outlet_boundary"]["pressure_Pa"] == bench["core_outlet_pressure_Pa"]
    assert cond["rccs"]["boundary_temperature_K"] == bench["rccs_boundary_temperature_K"]
    assert cond["comparison"]["core_outlet_temperature_K"] == bench["core_outlet_temperature_K"]


def test_system_loop_values_match_operation():
    cond = _s0()["conditions"]["system_loop"]
    op = _op()
    closure = cond["pressure_closure"]
    assert closure["type"] == op["helium_inventory"]["steady_state_constraint"]
    assert closure["location"] == op["primary"]["initial_pressure_location"]
    assert closure["pressure_Pa"] == op["primary"]["initial_pressure_Pa"]

    circ = cond["circulator"]
    assert circ["mass_flow_kg_s"] == op["primary"]["initial_mass_flow_kg_s"]
    assert circ["heat_to_helium_W"] == op["compressor"]["rated_heat_to_helium_W"]
    assert math.isclose(circ["rotation_rate_rad_s"], op["compressor"]["initial_rotation_rate_rad_s"])

    sink = cond["sg_heat_sink"]
    op_sink = op["boundaries"]["sg_heat_sink"]
    assert sink["model_level"] == op_sink["model_level"]
    assert sink["mode"] == op_sink["mode"] == "fixed_temperature_liquid_water"
    assert sink["water_temperature_K"] == op_sink["water_temperature_K"]
    assert sink["calibrate"] in op_sink and op_sink[sink["calibrate"]] is None, "U는 정상해에서 보정"
    assert sink["calibration_target_helium_outlet_K"] == op_sink["rated_helium_outlet_temperature_K"]
    hx = _geo()["heat_exchanger"]["helium_side"]
    assert set(sink["nodes"]) == {hx["fsh_node_id"], hx["ees_node_id"]}

    comp = cond["comparison"]
    psid = op["reference"]["psid_R19"]
    assert comp["sg_inlet_pressure_Pa"] == psid["pressure_Pa"]["sg_inlet"]
    assert comp["circulator_inlet_pressure_Pa"] == psid["pressure_Pa"]["circulator_inlet"]
    assert comp["circulator_pressure_rise_Pa"] == _geo()["pressure_loss_calibration"]["circulator_pressure_rise_Pa"]
    assert comp["sg_heat_removal_W"] == psid["sg_net_heat_load_W"]
    assert comp["helium_inventory_total_estimate_kg"] == op["helium_inventory"]["reference_kg"]
    assert op["helium_inventory"]["stagnant_volume_temperature_K"] > 273.15


def test_pressure_bases_differ_only_by_documented_offset():
    """두 조건 모두 하부 플레넘 기준. 차이는 보고서 D2의 벤치마크 vs PSID 차이(≈52.7 kPa)."""
    conds = _s0()["conditions"]
    bench = conds["core_benchmark"]["outlet_boundary"]
    loop = conds["system_loop"]["pressure_closure"]
    assert bench["pressure_location"] == loop["location"] == "lower_plenum"
    offset = bench["pressure_Pa"] - loop["pressure_Pa"]
    assert math.isclose(offset, 52700.0, abs_tol=1.0)
    assert offset / loop["pressure_Pa"] < 0.01


def _walk(node, path=()):
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, path + (str(key),))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _walk(value, path + (str(i),))
    else:
        yield path, node


def test_units_and_missing_inputs():
    """K는 절대온도, 압력·열량·유량은 양수. outputs 목록의 이름 문자열은 키가 아니라 검사 대상이 아니다."""
    positive_units = ("_Pa", "_W", "_kg_s", "_rad_s", "_kg")
    checked = 0
    for path, value in _walk(_s0()["conditions"]):
        if value is None:
            continue
        unit_key = next((p for p in reversed(path) if p.endswith(("_K",) + positive_units)), None)
        if unit_key is None:
            continue
        if unit_key.endswith("_K"):
            assert isinstance(value, Number) and value > 273.15, path
        else:
            assert isinstance(value, Number) and value > 0, path
        checked += 1
    assert checked >= 20   # 현재 단위 키 값 22개 — 순회가 값을 놓치지 않았는지 확인

    for item in _s0()["missing_inputs"]:
        assert item["item"] and item["owner"]
        assert set(item["conditions"]) <= set(_s0()["conditions"]), item
