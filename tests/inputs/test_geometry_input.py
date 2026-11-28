"""config/geometry.yaml 입력 검증 (Role 1).

YAML 구문, 기존 스키마 키, 셀 (N,)/면 (N+1,) 배열 길이, 양의 기하,
노드·연결부 참조, 보고서 계산값의 자기일관성을 확인한다.
형상 입력의 정합성 시험이며 실제 HTGR 물리 정확도 시험이 아니다.
"""
import math
from functools import lru_cache
from numbers import Number
from pathlib import Path

import numpy as np
import yaml

GEOMETRY = Path(__file__).resolve().parents[2] / "config" / "geometry.yaml"


@lru_cache(maxsize=1)
def _geo():
    return yaml.safe_load(GEOMETRY.read_text(encoding="utf-8"))


def _node_ids(geo):
    ids = {n["node_id"] for n in geo["flow_network"]["nodes"].values()}
    ids |= {c["node_id"] for c in geo["core"]["channels"].values()}
    hx = geo["heat_exchanger"]
    ids |= {hx["helium_side"]["fsh_node_id"], hx["helium_side"]["ees_node_id"],
            hx["secondary_side"]["node_id"]}
    ids |= {n["node_id"] for n in geo["scs"]["nodes"].values()}
    ids |= {b["node_id"] for b in geo["boundaries"].values()}
    return ids


def test_legacy_schema_keys_present_and_typed():
    geo = _geo()
    core = geo["core"]
    n = core["axial_cells"]
    assert isinstance(n, int) and n > 0
    assert core["channel_length_m"] > 0
    assert core["flow_area_m2"] > 0
    assert len(core["elevation_profile_m"]) == n + 1

    hx = geo["heat_exchanger"]
    m = hx["axial_cells"]
    assert isinstance(m, int) and m > 0
    assert len(hx["heat_transfer_area_m2"]) == m
    assert all(a > 0 for a in hx["heat_transfer_area_m2"])

    location = geo["compressor"]["location"]
    assert isinstance(location, str)
    assert location in geo["flow_network"]["nodes"]


def test_core_elevation_profile_matches_channel_and_plena():
    geo = _geo()
    core = geo["core"]
    z = np.asarray(core["elevation_profile_m"], dtype=float)
    assert np.all(np.diff(z) < 0), "노심은 하향류이므로 유동 순서 face 고도가 감소해야 한다"
    assert math.isclose(z[0] - z[-1], core["channel_length_m"], abs_tol=1e-4)

    nodes = geo["flow_network"]["nodes"]
    assert math.isclose(z[0], nodes["upper_plenum"]["z_bottom_m"], abs_tol=1e-4)
    assert math.isclose(z[-1], nodes["lower_plenum"]["z_top_m"], abs_tol=1e-4)

    regions = core["cell_region"]
    assert len(regions) == core["axial_cells"]
    fuel = [i for i, r in enumerate(regions) if r == "fuel"]
    assert len(fuel) == core["layout"]["blocks_per_column"]
    assert fuel == list(range(fuel[0], fuel[-1] + 1)), "연료 셀은 연속해야 한다"

    active = geo["rpv_internals"]["axial_regions_m"]["active_core"]
    assert math.isclose(z[fuel[0]], active[1], abs_tol=1e-4)
    assert math.isclose(z[fuel[-1] + 1], active[0], abs_tol=1e-4)
    dz = -np.diff(z)
    assert np.allclose(dz[fuel], core["layout"]["block_height_m"], atol=2e-4)


def test_core_channel_areas_volumes_and_block_holes():
    geo = _geo()
    core = geo["core"]
    ch = core["channels"]
    layout = core["layout"]
    block = core["fuel_block"]

    n_cols = layout["fuel_columns_standard"] + layout["fuel_columns_rsc"]
    assert ch["average"]["fuel_columns"] + ch["hot"]["fuel_columns"] == n_cols

    def hole_area(kind):
        large = block["coolant_holes_large"][kind] * math.pi * block["coolant_hole_large_radius_m"] ** 2
        small = block["coolant_holes_small"][kind] * math.pi * block["coolant_hole_small_radius_m"] ** 2
        return large + small

    per_block = block["coolant_flow_area_per_block_m2"]
    for kind in ("standard", "rsc"):
        assert math.isclose(hole_area(kind), per_block[kind], rel_tol=1e-3)

    total = (layout["fuel_columns_standard"] * per_block["standard"]
             + layout["fuel_columns_rsc"] * per_block["rsc"])
    assert math.isclose(core["flow_area_m2"], total, rel_tol=1e-3)
    assert math.isclose(ch["average"]["flow_area_m2"] + ch["hot"]["flow_area_m2"],
                        core["flow_area_m2"], rel_tol=1e-3)
    assert math.isclose(ch["hot"]["flow_area_m2"], per_block["standard"], rel_tol=1e-3)

    length = core["channel_length_m"]
    for name in ("average", "hot"):
        c = ch[name]
        assert math.isclose(c["volume_m3"], c["flow_area_m2"] * length, rel_tol=0.01)


def test_core_cell_flow_areas_follow_r1_porosity_and_reporting_layers():
    """하부반사체 셀 면적 = R1 Table I.9 공극률 × 육각 블록 면적, 연료 셀은 채널 면적. 보고 층은 Fig. II.7."""
    core = _geo()["core"]
    n = core["axial_cells"]
    regions = core["cell_region"]
    layout = core["layout"]
    porosity = core["block_porosity"]
    a_hex = core["hex_block_area_m2"]
    assert math.isclose(a_hex, math.sqrt(3) / 2 * layout["block_flat_to_flat_m"] ** 2, rel_tol=1e-4)

    lower = regions.index("lower_reflector")
    std_avg = layout["fuel_columns_standard"] - core["channels"]["hot"]["fuel_columns"]
    expected = {
        "average": (std_avg * porosity["bottom_reflector_standard"]
                    + layout["fuel_columns_rsc"] * porosity["bottom_reflector_rsc"]) * a_hex,
        "hot": porosity["bottom_reflector_standard"] * a_hex,
    }
    dz = -np.diff(np.asarray(core["elevation_profile_m"], dtype=float))
    for name, lower_area in expected.items():
        ch = core["channels"][name]
        area = np.asarray(ch["cell_flow_area_m2"], dtype=float)
        assert len(area) == n and np.all(area > 0), name
        assert all(area[i] == ch["flow_area_m2"] for i, r in enumerate(regions) if r == "fuel"), name
        assert math.isclose(area[lower], lower_area, rel_tol=1e-4), name
        assert math.isclose(float(area @ dz), ch["volume_m3"], rel_tol=0.01), name

    # 하부 블록 (PSID): 블록당 유로면적 = 유로 수 × π/4 × d². 평균채널은 65기둥 (RSC도 같게)
    channels = core["lower_block_channels"]
    per_block = {
        "lower_transition_reflector": channels["transition_block"]["count"]
        * math.pi / 4 * channels["transition_block"]["diameter_m"] ** 2,
        "flow_distribution_block": channels["flow_distribution_block"]["count"]
        * math.pi / 4 * channels["flow_distribution_block"]["diameter_m"] ** 2,
        "support_post_block": channels["post_block"]["vertex_circle_count_per_block"]
        * math.pi / 4 * channels["post_block"]["diameter_m"] ** 2,
    }
    porosity_key = {"lower_transition_reflector": "transition_block",
                    "flow_distribution_block": "flow_distribution_block",
                    "support_post_block": "post_block"}
    n_avg = core["channels"]["average"]["fuel_columns"]
    for region, a_block in per_block.items():
        i = regions.index(region)
        assert math.isclose(core["channels"]["hot"]["cell_flow_area_m2"][i], a_block, rel_tol=1e-3), region
        assert math.isclose(core["channels"]["average"]["cell_flow_area_m2"][i], n_avg * a_block, rel_tol=1e-3), region
        assert math.isclose(porosity[porosity_key[region]], a_block / a_hex, rel_tol=1e-3), region
    # 116 mm = 67 mm × √3 → 꼭짓점 원 2개분 = 67 mm 원 6개 (설계가 유로면적을 유지)
    assert math.isclose(per_block["flow_distribution_block"], per_block["support_post_block"], rel_tol=2e-3)
    # 상부반사체 = 연료 냉각공 배열 (R1 §I.5.1.4)
    for i, r in enumerate(regions):
        if r.startswith("upper_reflector"):
            assert core["channels"]["average"]["cell_flow_area_m2"][i] == core["channels"]["average"]["flow_area_m2"]

    layers = core["r1_reporting_layer"]
    assert len(layers) == n and layers[0] is None, "CRE는 R1 보고 메시 밖"
    assert all(a >= b for a, b in zip(layers[1:], layers[2:]))
    assert [layers[i] for i, r in enumerate(regions) if r == "fuel"] == list(range(12, 2, -1))
    # R1 Fig. II.7: layer 13은 ½블록(39.65 cm), layer 14는 1블록(79.30 cm)
    assert math.isclose(dz[layers.index(13)], 0.3965, abs_tol=2e-4)
    assert math.isclose(dz[layers.index(14)], 0.7930, abs_tol=2e-4)


def test_bypass_type_I_fractions_sum_to_eleven_percent():
    fractions = _geo()["core"]["channels"]["bypass"]["type_I_fraction"]
    assert len(fractions) == 7
    assert all(f > 0 for f in fractions.values())
    assert math.isclose(sum(fractions.values()), 0.11, abs_tol=1e-9)


def test_sg_area_cells_and_tubes_are_consistent():
    hx = _geo()["heat_exchanger"]
    n = hx["axial_cells"]
    tubes = hx["tubes"]
    area = np.asarray(hx["heat_transfer_area_m2"], dtype=float)
    total = hx["heat_transfer_area_total_m2"]

    assert math.isclose(area.sum(), total, rel_tol=1e-5)
    outer = tubes["count"] * math.pi * tubes["outer_diameter_m"] * tubes["effective_length_m"]
    assert math.isclose(outer, total, rel_tol=1e-3)

    height = np.asarray(hx["helium_side"]["cell_height_m"], dtype=float)
    assert len(height) == n
    bundle = hx["bundle"]
    assert math.isclose(height.sum(), bundle["ees_height_m"] + bundle["fsh_height_m"], rel_tol=1e-6)

    sec = hx["secondary_side"]
    mapping = sec["cell_to_helium_cell"]
    assert sorted(mapping) == list(range(n))
    assert mapping == list(range(n - 1, -1, -1)), "SG는 대향류: 2차측 셀 순서가 He 셀 순서의 역순"
    tube_len = np.asarray(sec["tube_length_per_cell_m"], dtype=float)
    assert len(tube_len) == n
    assert math.isclose(tube_len.sum(), tubes["effective_length_m"], rel_tol=1e-4)
    per_tube_circumference = tubes["count"] * math.pi * tubes["outer_diameter_m"]
    assert np.allclose(tube_len, area[mapping] / per_tube_circumference, rtol=1e-3)

    d_i = tubes["outer_diameter_m"] - 2 * tubes["wall_thickness_m"]
    assert math.isclose(sec["hydraulic_diameter_m"], d_i, rel_tol=1e-4)
    assert math.isclose(sec["flow_area_m2"], tubes["count"] * math.pi * d_i ** 2 / 4, rel_tol=1e-3)
    assert math.isclose(sec["volume_m3"], sec["flow_area_m2"] * sec["length_m"], rel_tol=1e-3)

    # 관 내경·내면적·벽 체적도 같은 PSID 관 세트(350관 × 164.5 m)에서 계산
    n_t, L_t, d_o = tubes["count"], tubes["effective_length_m"], tubes["outer_diameter_m"]
    assert n_t == 350, "SG 관 수는 PSID 기준 (팀 결정 2026-11-26)"
    assert math.isclose(tubes["inner_diameter_m"], d_i, rel_tol=1e-4)
    assert math.isclose(tubes["inner_surface_area_total_m2"], n_t * math.pi * d_i * L_t, rel_tol=1e-4)
    assert math.isclose(tubes["wall_volume_m3"], n_t * math.pi / 4 * (d_o ** 2 - d_i ** 2) * L_t, rel_tol=1e-4)


def test_radial_structure_is_strictly_outward():
    radii = list(_geo()["core"]["radial_m"].values())
    assert all(r > 0 for r in radii)
    assert all(a < b for a, b in zip(radii, radii[1:]))


def test_rpv_axial_regions_are_contiguous_from_outlet_plenum_up():
    regions = _geo()["rpv_internals"]["axial_regions_m"]
    names = list(regions)
    stack = names[names.index("outlet_plenum"):]
    for name in names:
        bottom, top = regions[name]
        assert bottom < top, name
    for lower, upper in zip(stack, stack[1:]):
        assert math.isclose(regions[lower][1], regions[upper][0], abs_tol=1e-4), (lower, upper)


def test_junctions_and_heat_structures_reference_defined_nodes():
    geo = _geo()
    ids = _node_ids(geo)
    junctions = geo["junctions"]
    for jid, j in junctions.items():
        assert j["from"] in ids, (jid, j["from"])
        assert j["to"] in ids, (jid, j["to"])
        assert j["from"] != j["to"], jid

    connected = {j["from"] for j in junctions.values()} | {j["to"] for j in junctions.values()}
    hydraulic = ({n["node_id"] for n in geo["flow_network"]["nodes"].values()}
                 | {c["node_id"] for c in geo["core"]["channels"].values()})
    assert hydraulic <= connected, hydraulic - connected

    for seg_name, seg in geo["pressure_loss_calibration"]["segments"].items():
        for jid in seg["junctions"]:
            assert jid in junctions, (seg_name, jid)

    structures = geo["heat_structures"]
    for hs_id, hs in structures.items():
        for side in ("inner", "outer"):
            refs = hs[side]
            refs = refs if isinstance(refs, list) else [refs]
            for ref in refs:
                if isinstance(ref, int):
                    assert ref in ids, (hs_id, side, ref)
                elif isinstance(ref, str) and ref.startswith("HS-"):
                    assert ref in structures, (hs_id, side, ref)


PRIMARY_LOOP_JUNCTIONS = [f"J{i:02d}" for i in range(1, 20)]   # J01–J19: 순환기 → 노심 → SG → 순환기


def _node_elevations(geo):
    """노드 번호 → (z_inlet, z_outlet). 1차 루프 수력 노드."""
    z = {n["node_id"]: (n["z_inlet_m"], n["z_outlet_m"]) for n in geo["flow_network"]["nodes"].values()}
    core_z = geo["core"]["elevation_profile_m"]
    for channel in geo["core"]["channels"].values():
        z[channel["node_id"]] = (core_z[0], core_z[-1])
    hs = geo["heat_exchanger"]["helium_side"]
    sg_z = hs["elevation_profile_m"]
    z[hs["fsh_node_id"]] = (sg_z[0], sg_z[1])
    z[hs["ees_node_id"]] = (sg_z[1], sg_z[-1])
    return z


def test_primary_loop_elevations_close_at_every_junction():
    """M08: 연결부마다 앞 노드 출구 고도 = 뒤 노드 입구 고도 → 루프 Σ Δz = 0."""
    geo = _geo()
    z = _node_elevations(geo)
    for jid in PRIMARY_LOOP_JUNCTIONS:
        j = geo["junctions"][jid]
        z_out, z_in = z[j["from"]][1], z[j["to"]][0]
        assert z_out is not None and z_in is not None, jid
        assert math.isclose(z_out, z_in, abs_tol=1e-3), (jid, z_out, z_in)

    for n in geo["flow_network"]["nodes"].values():
        if n.get("dz_m") is not None:
            assert math.isclose(n["dz_m"], n["z_outlet_m"] - n["z_inlet_m"], abs_tol=0.02), n["node_id"]
        if n.get("z_bottom_m") is not None and n.get("z_top_m") is not None:
            assert n["z_bottom_m"] < n["z_top_m"], n["node_id"]


def test_sg_elevations_match_bundle_layout_and_vessel():
    geo = _geo()
    hx = geo["heat_exchanger"]
    hs = hx["helium_side"]
    z = np.asarray(hs["elevation_profile_m"], dtype=float)
    assert len(z) == hx["axial_cells"] + 1
    assert np.all(np.diff(z) < 0), "SG He 측은 하향류"

    expected = np.asarray(hs["cell_height_m"], dtype=float)
    expected[hs["fsh_ees_gap_attached_to_cell"]] += hx["bundle"]["fsh_ees_gap_height_m"]
    assert np.allclose(-np.diff(z), expected, atol=2e-3)

    vessel = hx["vessel"]
    assert math.isclose(vessel["z_top_m"] - vessel["z_bottom_m"], vessel["length_m"], abs_tol=0.05)

    pe = geo["plant_elevations"]
    nodes = geo["flow_network"]["nodes"]
    stack = [vessel["z_bottom_m"], pe["sg_subassembly_bottom_z_m"], pe["sg_shroud_bottom_z_m"],
             z[-1], z[0], pe["reference_plane_P_z_m"], pe["sg_subassembly_top_z_m"],
             vessel["z_top_m"], pe["circulator_top_z_m"]]
    assert all(a < b for a, b in zip(stack, stack[1:])), stack
    assert nodes["sg_inlet_plenum"]["z_top_m"] == pe["sg_subassembly_top_z_m"]
    assert nodes["sg_outlet_plenum"]["z_bottom_m"] == pe["sg_subassembly_bottom_z_m"]
    assert nodes["sg_annulus_riser"]["z_inlet_m"] == pe["sg_shroud_bottom_z_m"]
    assert nodes["main_circulator"]["z_inlet_m"] == vessel["z_top_m"]
    for name in ("cold_annulus", "hot_duct"):
        assert nodes[name]["z_inlet_m"] == nodes[name]["z_outlet_m"] == pe["reference_plane_P_z_m"]


def test_thermal_centers_follow_report_equation_6b():
    geo = _geo()
    tc = geo["plant_elevations"]["thermal_centers"]
    regions = geo["rpv_internals"]["axial_regions_m"]
    assert math.isclose(tc["core_z_m"], sum(regions["active_core"]) / 2, abs_tol=0.01)

    hs = geo["heat_exchanger"]["helium_side"]
    z = hs["elevation_profile_m"]
    fsh_bottom = z[0] - hs["cell_height_m"][0]
    assert math.isclose(tc["fsh_center_z_m"], (z[0] + fsh_bottom) / 2, abs_tol=0.01)
    assert math.isclose(tc["ees_center_z_m"], (z[1] + z[-1]) / 2, abs_tol=0.01)

    f = tc["fsh_duty_fraction"]
    z_sg = f * tc["fsh_center_z_m"] + (1 - f) * tc["ees_center_z_m"]
    assert math.isclose(tc["sg_z_m"], z_sg, abs_tol=0.01)
    assert math.isclose(tc["delta_z_m"], tc["core_z_m"] - z_sg, abs_tol=0.1)
    # SG 번들 전체가 노심보다 낮다 (R19 §5.3.3.4.1): FSH 상단 < 출구 플레넘 바닥
    assert z[0] < regions["outlet_plenum"][0]


def test_pressure_loss_targets_close_the_loop():
    cal = _geo()["pressure_loss_calibration"]
    design = sum(s["design_Pa"] for s in cal["segments"].values())
    initial = sum(s["initial_Pa"] for s in cal["segments"].values())
    margin = sum(s["margin_Pa"] for s in cal["segments"].values())
    assert math.isclose(design, cal["circulator_pressure_rise_Pa"], abs_tol=500.0)
    assert math.isclose(initial + margin, design, abs_tol=500.0)


def _numbers(value):
    if isinstance(value, bool):
        return
    if isinstance(value, Number):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from _numbers(item)


def _walk(node, path=()):
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, path + (str(key),))
    else:
        yield path, node


def test_lengths_areas_volumes_are_positive_or_null():
    """null은 미정이며 0으로 채우지 않는다. 고도(z, dz, elevation; 음수 허용)와 원점 포함 반경만 예외."""
    exempt = {"r_m", "surface_roughness_m"}
    checked = 0
    for path, value in _walk(_geo()):
        key = path[-1]
        if not key.endswith(("_m", "_m2", "_m3")) or key in exempt:
            continue
        if key.startswith(("z", "dz")) or key.endswith("_z_m") or "elevation" in key:
            continue
        if value is None:
            continue
        numbers = list(_numbers(value))
        assert numbers, path
        assert all(x > 0 for x in numbers), path
        checked += 1
    assert checked > 50
