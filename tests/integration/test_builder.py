"""Role 5 builder 단독 시험: 입력 검사·가정 파일·연결부 표 조립·p_sys 소유·정상 폐합·포트 연결·출력 평탄화.

풀이는 하지 않는다 (조립만, 각 0.1 s 안팎). 실제 정상해는 test_s0_core_benchmark.py / test_s0_system_loop.py.
"""
import copy
import dataclasses

import numpy as np
import pytest
import yaml

from interfaces.state_registry import StateRegistry
from models._inputs import ModelInputMissing
from system.builder import MissingInputsError, build_system, check_inputs, load_inputs
from system.simulation import flatten_names, settings_from_scenario

ASSUME = ["config/assumptions/provisional_s0.yaml"]


@pytest.fixture(scope="module")
def loop_inputs():
    return load_inputs("steady_state", "system_loop", assumptions=ASSUME)


def test_missing_inputs_are_listed_with_file_key_and_condition():
    with pytest.raises(MissingInputsError) as err:
        build_system(load_inputs("steady_state", "system_loop"))
    msg = str(err.value)
    assert "scenario 'steady_state' condition 'system_loop'" in msg
    keys = {(m.file, m.key) for m in err.value.items}
    # 남은 M10: SG 용기 상부 영역을 110·350에 나누는 비율 (#31에서 나머지는 R19·R1 치수로 채움)
    assert keys == {("config/geometry.yaml", "flow_network.nodes.mc_discharge_plenum.volume_m3"),
                    ("config/geometry.yaml", "flow_network.nodes.mc_inlet_mlsv.volume_m3")}
    assert isinstance(err.value, ModelInputMissing)
    # 노심 단독 조건은 가정 없이 원문·계산 입력만으로 조립된다
    assert check_inputs(load_inputs("steady_state", "core_benchmark")) == []
    # null을 하나 만들면 파일·키로 보고한다
    inputs = load_inputs("steady_state", "core_benchmark")
    geo = copy.deepcopy(inputs.geometry)
    geo["core"]["channels"]["bypass"]["hydraulic_diameter_m"] = None
    del geo["core"]["channels"]["bypass"]["cell_hydraulic_diameter_m"]
    missing = check_inputs(dataclasses.replace(inputs, geometry=geo))
    assert [(m.file, m.key) for m in missing] == [("config/geometry.yaml", "core.channels.bypass.hydraulic_diameter_m")]


def test_assumptions_fill_only_nulls_and_are_recorded(tmp_path, loop_inputs):
    assert check_inputs(loop_inputs) == []
    keys = {a.key: a.value for a in loop_inputs.assumptions}
    assert keys == {"geometry.flow_network.nodes.mc_discharge_plenum.volume_m3": 24.425,
                    "geometry.flow_network.nodes.mc_inlet_mlsv.volume_m3": 24.425}
    assert loop_inputs.files["assumptions[0]"] == "config/assumptions/provisional_s0.yaml"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump({"geometry": {"flow_network": {"nodes": {"riser": {"flow_area_m2": 9.9}}}}}),
                   encoding="utf-8")
    with pytest.raises(ValueError, match="would replace the sourced value"):
        load_inputs("steady_state", "system_loop", assumptions=[bad])
    odd = tmp_path / "odd.yaml"
    odd.write_text(yaml.safe_dump({"solver": {"dt": 1.0}}), encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported top-level keys"):
        load_inputs("steady_state", "system_loop", assumptions=[odd])


def test_loop_is_assembled_from_the_junction_table(loop_inputs):
    built = build_system(loop_inputs)
    assert built.flow_order == ("main_circulator", "mc_discharge_plenum", "cold_annulus", "rv_inlet_mcss", "riser",
                                "upper_plenum", "core", "lower_plenum", "hot_duct", "sg_inlet_plenum", "sg",
                                "sg_outlet_plenum", "sg_annulus_riser", "mc_inlet_mlsv")
    names = [c[0] for c in built.connections]
    assert "J06_J08" in names and "J09_J11" in names          # 병렬 노심 유로 → 연결 하나
    assert not any("J15" in c[3] for c in built.connections)  # SG 내부 FSH → EES
    assert not any(j in ("J20", "J21", "J22", "J23") for c in built.connections for j in c[3])   # SCS 닫힘
    # Core·HeatExchanger 상위만 Assembly에 등록, 하위 상태 이중 등록 없음
    comps = [c.name for c in built.system.components]
    assert "core" in comps and "sg" in comps and "core.fluid" not in comps
    assert len({id(c) for c in built.inventory}) == len(built.inventory) == 13
    assert built.parts["fluid"] in built.inventory and built.parts["sg_he"] in built.inventory
    assert {a.key for a in built.assumptions_used} == {a.key for a in loop_inputs.assumptions}
    assert any("HS-09" in s for s in built.model_assumptions) and any("350" in s for s in built.model_assumptions)


def test_open_domain_boundaries_and_assumption_scope():
    built = build_system(load_inputs("steady_state", "core_benchmark", assumptions=ASSUME))
    comps = [c.name for c in built.system.components]
    assert comps[:3] == ["bc_p_sys", "bc_J03", "bc_J12"]
    assert built.p_sys_owner == "bc_p_sys" and built.replaced_residual_id is None
    assert built.system.registry.size == 347
    assert built.assumptions_used == ()                     # 110·350 가정은 노심 단독 영역 밖
    assert "geometry.flow_network.nodes.mc_discharge_plenum.volume_m3" in built.describe()["assumptions_outside_domain"]


def test_p_sys_single_owner_and_steady_closure(loop_inputs):
    steady = build_system(loop_inputs, mode="steady")
    ids = steady.system.residual_ids()
    spec = steady.system.registry.spec("p_sys")
    assert spec.kind == "differential" and spec.owner == "primary_loop"
    assert steady.p_sys_owner == "lower_plenum" and steady.replaced_residual_id == "lower_plenum.mass"
    assert "lower_plenum.steady_pressure_closure" in ids and "lower_plenum.mass" not in ids
    assert "lower_plenum.pressure_reference" in ids                     # 국부 기준식은 정상 폐합과 별개
    with pytest.raises(ValueError, match="U_sg"):
        build_system(loop_inputs, mode="transient")
    transient = build_system(loop_inputs, mode="transient", U_sg=470.0)
    t_ids = transient.system.residual_ids()
    assert "lower_plenum.mass" in t_ids and not any(i.endswith("steady_pressure_closure") for i in t_ids)
    # 과도 계 = 정상 계 − SG-0 U 보정 미지수
    assert set(steady.system.registry.names) - set(transient.system.registry.names) == {"sg.secondary.U_overall"}
    # 시나리오가 교체 ID를 다르게 적으면 실패
    cond = copy.deepcopy(loop_inputs.condition)
    cond["pressure_closure"]["replaced_residual_id"] = "upper_plenum.mass"
    with pytest.raises(ValueError, match="replaced_residual_id"):
        build_system(dataclasses.replace(loop_inputs, condition=cond))


def test_domain_errors_are_reported(loop_inputs):
    cond = copy.deepcopy(loop_inputs.condition)
    cond["domain"]["nodes"].remove(110)                                 # J01·J02가 영역을 자름 → 경계 아님
    with pytest.raises(ValueError, match="neither closed nor a boundary"):
        build_system(dataclasses.replace(loop_inputs, condition=cond))
    cond = copy.deepcopy(loop_inputs.condition)
    cond["domain"]["nodes"].remove(210)                                 # 노심 유로 일부만
    with pytest.raises(ValueError, match="all core nodes"):
        build_system(dataclasses.replace(loop_inputs, condition=cond))
    cond = copy.deepcopy(loop_inputs.condition)
    cond["sg_heat_sink"]["model_level"] = "SG-1"
    with pytest.raises(NotImplementedError, match="SG-0"):
        build_system(dataclasses.replace(loop_inputs, condition=cond))
    bench = load_inputs("steady_state", "core_benchmark", assumptions=ASSUME)
    cond = copy.deepcopy(bench.condition)
    cond["default_subcase"] = "2c"
    with pytest.raises(NotImplementedError, match="fluence"):
        build_system(dataclasses.replace(bench, condition=cond))


def test_run_settings_come_from_the_scenario(loop_inputs):
    s = settings_from_scenario(loop_inputs)
    assert s.sparse_jacobian and s.steady_residual_max == 1e-6 and s.hold.t_end_s == 10.0
    scen = copy.deepcopy(loop_inputs.scenario)
    del scen["solution"]["newton"]
    with pytest.raises(ModelInputMissing, match=r"steady_state\.yaml: missing solution\.newton"):
        settings_from_scenario(dataclasses.replace(loop_inputs, scenario=scen))


def test_flatten_names_uses_registry_shape_and_detects_duplicates():
    reg = StateRegistry()
    reg.register("a", kind="algebraic", initial=np.zeros(2), unit="K", shape=(2,))
    reg.register("b", kind="differential", initial=1.0, unit="Pa")
    reg.freeze()
    assert flatten_names(reg) == (["a[0]", "a[1]", "b"], ["K", "K", "Pa"])
    reg = StateRegistry()
    reg.register("a", kind="algebraic", initial=np.zeros(1), unit="K", shape=(1,))
    reg.register("a[0]", kind="algebraic", initial=0.0, unit="K")
    reg.freeze()
    with pytest.raises(ValueError, match="duplicate flattened"):
        flatten_names(reg)


def test_solid_sensible_energy_integrates_temperature_dependent_capacity():
    from validation.conservation import SOLID_T_REF_K, _sensible_energy
    C = lambda T: 2.0 + 0.01 * np.asarray(T)            # noqa: E731  C(T) 선형 → ∫ = 2ΔT + 0.005(T² − T_ref²)
    T = np.array([500.0, 1200.0])
    exact = np.sum(2.0 * (T - SOLID_T_REF_K) + 0.005 * (T ** 2 - SOLID_T_REF_K ** 2))
    assert _sensible_energy(C, T) == pytest.approx(exact, rel=1e-12)
