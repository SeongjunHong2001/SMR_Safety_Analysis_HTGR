"""Role 5 압력손실 K 보정 (S0 system_loop 정상해, 보고서 식 7).

목표 (geometry.yaml ``pressure_loss_calibration.segments``, PSID R19 Table 5.3-5 설계손실)
    구간마다 비가역 손실 Δp_irr = Σ(마찰 + 형상손실) = design_Pa. 정수압 ρgΔz와 가역 동압 변화(플레넘 베르누이,
    열가속)는 목표에 넣지 않고 노드 고도·상태로 따로 계산한다 (geometry.yaml 주석의 [가정]).
보정 변수
    구간마다 K 하나. 구간 연결부 목록에서 면을 실을 수 있는 첫 연결부에 둔다 (builder.junction_face).
    병렬 노심 유로 묶음(J06–J08)은 같은 K를 연료 유로 연결부 모두에 주고, Type I 고정 분율 우회 유로는 제외한다
    (우회 유량은 분율식과 dp_bypass_cal이 정한다).
갱신 (식 7)
    K_j ← K_j + (Δp_target − Δp_irr) / q_j,   q_j = ṁ|ṁ| / (2 ρ A²) (K 면, 대표 유로)
    저마하수 모델에서 밀도는 p_sys와 엔탈피로만 정해지고 순환기 유량이 지정이라, 비노심 구간은 한 번에 맞는다.
    노심은 평균·고온 유로 분배가 조금 바뀌어 두세 번 반복한다.
손실 배정
    유로 면의 마찰과 내부 면 형상손실 → 그 노드의 구간 (목록에 없는 노드는 출구 연결부의 구간),
    경계 면 형상손실 → 그 면의 연결부 구간. 노심은 유량이 가장 큰 유로(평균 유로)를 경로로 본다.
결과
    보정 K는 별도 파일(config/calibration/*.yaml)에 쓴다. builder는 이 파일을 넘길 때만 쓴다 (``--calibration``).
"""
from __future__ import annotations

import copy
import dataclasses
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from models._inputs import ModelInputMissing
from system.builder import (ROOT, ScenarioInputs, _get, _units, _connections, build_system, core_channel_nodes,
                            junction_face)
from system.builder import fingerprint_of
from system.simulation import _code_version, settings_from_scenario, solve_steady


class CalibrationError(RuntimeError):
    pass


@dataclass(frozen=True)
class Segment:
    name: str
    design_Pa: float
    junctions: tuple
    nodes: tuple
    k_junctions: tuple          # 보정 K를 받는 연결부 (같은 값)


def segments(inputs: ScenarioInputs, units=None, unit_of=None):
    """보정 구간과 각 구간의 K 연결부."""
    cal = _get(inputs.geometry, "pressure_loss_calibration.segments")
    if not cal:
        raise ModelInputMissing(f"{inputs.files['geometry']}: pressure_loss_calibration.segments")
    if units is None:
        units, unit_of = _units(inputs)
    conns = _connections(inputs, units, unit_of)
    group_of = {j: c.junctions for c in conns for j in c.junctions}
    bypass = (_get(inputs.condition, "bypass_model") == "type_I_fixed_fraction")
    out = []
    for name, seg in cal.items():
        design = seg.get("design_Pa")
        if design is None:
            raise ModelInputMissing(f"{inputs.files['geometry']}: pressure_loss_calibration.segments.{name}.design_Pa")
        k_j = None
        for j in seg.get("junctions") or ():
            try:
                face = junction_face(inputs, units, unit_of, j)
            except (NotImplementedError, ValueError):
                continue                                   # 내부 연결(J15)·면 없는 연결(J19)
            if face is None:
                continue
            group = group_of.get(j, (j,))
            k_j = tuple(g for g in group
                        if not (bypass and junction_face(inputs, units, unit_of, g)[2] == "bypass"))
            break
        if not k_j:
            raise CalibrationError(f"segment {name}: no junction with a pipe/core/SG face to carry K")
        out.append(Segment(name, float(design), tuple(seg.get("junctions") or ()), tuple(seg.get("nodes") or ()),
                           k_j))
    return out


def _sg_nodes(inputs):
    hs = inputs.geometry["heat_exchanger"]["helium_side"]
    return [int(hs["fsh_node_id"]), int(hs["ees_node_id"])]


def segment_losses(built, y, segs=None):
    """구간별 비가역 손실(마찰+형상)과 유로 면의 정수압·가속 항 [Pa], 그리고 K 면 동압 q [Pa].

    gravity·acceleration은 배관·노심·SG 면 항만 더한 참고값이다 (플레넘 포트의 정수압은 빠진다).
    """
    inputs, ctx = built.inputs, built.system.context
    segs = segs or segments(inputs)
    seg_of_j = {j: s.name for s in segs for j in s.junctions}
    seg_of_n = {int(n): s.name for s in segs for n in s.nodes}
    res = {s.name: {"irreversible": 0.0, "gravity": 0.0, "acceleration": 0.0} for s in segs}
    res["unassigned"] = {"irreversible": 0.0, "gravity": 0.0, "acceleration": 0.0}
    q_face = {}
    core_ch_node = {v: k for k, v in core_channel_nodes(inputs.geometry).items()}
    juncs = inputs.geometry["junctions"]

    def add(seg, kind, value):
        res[seg if seg in res else "unassigned"][kind] += float(value)

    def assign(terms, node_seg, j_in, j_out, rep_q_key):
        n = terms["fric"].size
        for f in range(n):
            add(node_seg, "irreversible", terms["fric"][f])
            add(node_seg, "gravity", terms["grav"][f])
            if f == 0:
                add(seg_of_j.get(j_in), "irreversible", terms["form"][f])
                add(seg_of_j.get(j_in), "acceleration", terms["acc"][f])
            elif f == n - 1:
                add(seg_of_j.get(j_out), "irreversible", terms["form"][f])
                add(seg_of_j.get(j_out), "acceleration", terms["acc"][f])
            else:
                add(node_seg, "irreversible", terms["form"][f])
                add(node_seg, "acceleration", terms["acc"][f])
        for j, f in ((j_in, 0), (j_out, -1)):
            m, rho, A = terms["m"][f], terms["rho"][f], terms["A"][f]
            q_face[(rep_q_key, j)] = float(m * abs(m) / (2.0 * rho * A * A))

    for key in built.flow_order:
        if key == "core":
            fluid = built.parts["fluid"]
            terms_all = fluid.face_pressure_terms(0.0, y, ctx)
            for label, terms in terms_all.items():
                node = core_ch_node[label]
                j_in = next(j for j, d in juncs.items() if int(d["to"]) == node)
                j_out = next(j for j, d in juncs.items() if int(d["from"]) == node)
                if label == max(terms_all, key=lambda k: terms_all[k]["m"][0]):
                    assign(terms, seg_of_n.get(node), j_in, j_out, "core")
                else:                                     # 병렬 유로: K 면 q만 기록
                    for j, f in ((j_in, 0), (j_out, -1)):
                        m, rho, A = terms["m"][f], terms["rho"][f], terms["A"][f]
                        q_face[("core", j)] = float(m * abs(m) / (2.0 * rho * A * A))
            continue
        comp = built.parts["sg_he"] if key == "sg" else built.parts[key]
        if not hasattr(comp, "face_pressure_terms"):
            continue
        conn_in = next(c for c in built.connections if c[2] == key)
        conn_out = next(c for c in built.connections if c[1] == key)
        j_in, j_out = conn_in[3][-1], conn_out[3][0]
        node_ids = _sg_nodes(inputs) if key == "sg" else [int(inputs.geometry["flow_network"]["nodes"][key]["node_id"])]
        node_seg = next((seg_of_n[n] for n in node_ids if n in seg_of_n), seg_of_j.get(j_out))
        assign(comp.face_pressure_terms(0.0, y, ctx), node_seg, j_in, j_out, key)
    return res, q_face


@dataclass
class Calibration:
    K: dict                                   # 연결부 → K_form
    losses: dict                              # 구간 → {irreversible, gravity, acceleration} [Pa]
    targets: dict                             # 구간 → design_Pa
    iterations: int
    history: list = field(default_factory=list)
    built: object = None
    steady: object = None
    stations: list = field(default_factory=list)
    summary: dict = field(default_factory=dict)
    segments: list = field(default_factory=list)


def _base_forward_k(built, unit, junction, K_current):
    """K 면의 기본 정방향 경계 K (연결부 증분 제외)."""
    face = junction_face(built.inputs, *_units(built.inputs), junction)
    side = face[1]
    if unit == "core":
        fluid = built.parts["fluid"]
        return fluid.inlet_K[0] if side == "inlet" else fluid.outlet_K[0]
    comp = built.parts["sg_he"].pipe if unit == "sg" else built.parts[unit]
    return (comp.inlet if side == "inlet" else comp.outlet)[1] - K_current


def with_junction_k(inputs: ScenarioInputs, K: dict) -> ScenarioInputs:
    geo = copy.deepcopy(inputs.geometry)
    for j, v in K.items():
        geo["junctions"][j]["K_form"] = float(v)
    return dataclasses.replace(inputs, geometry=geo)


def station_pressures(built, y):
    """operation.initial_distribution 지점별 모델 압력 vs PSID (플레넘=셀, 유로=첫 셀, 순환기=토출 포트)."""
    reg = built.system.registry
    nodes = built.inputs.geometry["flow_network"]["nodes"]
    by_id = {int(c["node_id"]): k for k, c in nodes.items() if isinstance(c, dict) and "node_id" in c}
    rows = []
    for st in _get(built.inputs.operation, "initial_distribution.stations") or ():
        nid = int(st["node"])
        if nid in _sg_nodes(built.inputs):
            key, p = "sg", float(np.asarray(reg.get(y, "sg.helium.p_fluid_cell"))[0])
        elif nid in by_id and by_id[nid] in built.parts:
            key = by_id[nid]
            comp = built.parts[key]
            if key == "main_circulator":
                p = float(reg.get(y, comp.ports()["helium_outlet"].p))
            else:
                p = float(np.atleast_1d(np.asarray(reg.get(y, f"{key}.p_fluid_cell")))[0])
        else:
            continue
        rows.append({"station": st["station"], "node": nid, "point": st["point"], "unit": key,
                     "psid_Pa": float(st["p_Pa"]), "model_Pa": p, "difference_Pa": p - float(st["p_Pa"])})
    return rows


def calibrate_pressure_losses(inputs: ScenarioInputs, *, tol_Pa=1.0, max_iter=6, log=print) -> Calibration:
    """S0 system_loop 정상해에서 구간별 K를 식 7로 반복 보정한다."""
    if not inputs.closed_loop:
        raise CalibrationError("K calibration needs the closed-loop condition (system_loop)")
    settings = settings_from_scenario(inputs)
    segs = segments(inputs)
    K = {j: float(inputs.geometry["junctions"][j].get("K_form") or 0.0) for s in segs for j in s.k_junctions}
    targets = {s.name: s.design_Pa for s in segs}
    history = []
    for it in range(1, max_iter + 1):
        built = build_system(with_junction_k(inputs, K), mode="steady")
        steady = solve_steady(built, settings)
        losses, q = segment_losses(built, steady.y, segs)
        err = {s.name: s.design_Pa - losses[s.name]["irreversible"] for s in segs}
        history.append({"iteration": it, "max_error_Pa": max(abs(e) for e in err.values()), "K": dict(K)})
        log(f"  K iteration {it}: max |Δp_target − Δp_irr| = {history[-1]['max_error_Pa']:.3g} Pa "
            f"({steady.wall_s:.1f} s)")
        if history[-1]["max_error_Pa"] <= tol_Pa:
            break
        for s in segs:
            rep = s.k_junctions[0]
            unit = junction_face(inputs, *_units(inputs), rep)[0]
            qk = q[(unit, rep)]
            if not qk > 0:
                raise CalibrationError(f"segment {s.name}: no forward flow at the K face {rep}")
            new = K[rep] + err[s.name] / qk
            base = _base_forward_k(built, unit, rep, K[rep])
            if base + new < 0:
                raise CalibrationError(f"segment {s.name}: target {s.design_Pa:.0f} Pa is below the friction "
                                       f"already in the model ({losses[s.name]['irreversible']:.0f} Pa with "
                                       f"K_{rep} = {K[rep]:.3g}); the face K would be negative")
            for j in s.k_junctions:
                K[j] = new
    else:
        raise CalibrationError(f"K calibration did not converge in {max_iter} iterations: {history[-1]}")
    cal = Calibration(K=K, losses=losses, targets=targets, iterations=it, history=history, built=built, steady=steady,
                      segments=segs)
    cal.stations = station_pressures(built, steady.y)
    out = built.system.outputs(0.0, steady.y)
    rise = float(out["main_circulator.pressure_rise"])
    irr = float(sum(v["irreversible"] for v in losses.values()))
    cal.summary = {"circulator_pressure_rise_Pa": rise,
                   "sum_design_Pa": float(sum(targets.values())),
                   "sum_irreversible_Pa": irr,
                   # 승압 − 비가역 = 루프 정수압(플레넘 포함) + 열가속·단면변화 (가역 항 순합)
                   "reversible_net_Pa": rise - irr,
                   "unassigned_irreversible_Pa": float(losses["unassigned"]["irreversible"])}
    return cal


def write_calibration(cal: Calibration, path, *, log=print) -> Path:
    """보정 K를 builder가 읽는 덮어쓰기 파일로 쓴다 (geometry.junctions.<J>.K_form)."""
    path = Path(path)
    path = path if path.is_absolute() else ROOT / path
    inputs = cal.built.inputs
    assumption_files = sorted({a.source for a in inputs.assumptions
                               if not a.key.startswith("geometry.junctions.")})
    doc = {
        "meta": {
            "kind": "calibration",
            "quantity": "연결부 형상손실계수 K_form (정·역 대칭 [가정]), 무차원, K 면의 동압 기준",
            "method": "보고서 식 7: K ← K + (Δp_target − Δp_irr)/q, 구간 비가역 손실 = 마찰 + 형상손실",
            "basis": {"scenario": inputs.scenario_name, "condition": inputs.condition_name,
                      "targets": "geometry.yaml pressure_loss_calibration.segments.*.design_Pa (PSID R19 Table 5.3-5)"},
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "code_version": _code_version(),
            "depends_on_assumptions": assumption_files,
            # 유효 입력값 지문 (주석·개행 무관, K_form 제외). builder.load_inputs가 K 적용 전에 대조한다
            "input_fingerprint": fingerprint_of(inputs),
            "iterations": cal.iterations,
            "segments_Pa": {k: {"target": cal.targets[k], "irreversible": round(v["irreversible"], 3),
                                "gravity": round(v["gravity"], 3), "acceleration": round(v["acceleration"], 3)}
                            for k, v in cal.losses.items() if k in cal.targets},
            "circulator_pressure_rise_Pa": round(cal.summary["circulator_pressure_rise_Pa"], 1),
            "reversible_net_Pa": round(cal.summary["reversible_net_Pa"], 1),
            "note": "input_fingerprint와 같은 입력(형상·운전점·시나리오 조건·참조 물성·적용한 가정)에서만 유효하다. "
                    "다르면 builder가 적용을 거부하므로 python main.py --calibrate로 다시 만든다.",
        },
        "geometry": {"junctions": {j: {"K_form": round(float(v), 6)} for j, v in sorted(cal.K.items())}},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("# 자동 생성: python main.py --calibrate (system/calibration.py). 손으로 고치지 말고 다시 보정한다.\n"
              "# builder는 --calibration 으로 넘길 때만 이 K를 쓴다.\n")
    path.write_text(header + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")
    log(f"  wrote calibration: {path}")
    return path
