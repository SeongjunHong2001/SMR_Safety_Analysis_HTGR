"""Role 5 생산용 builder: config YAML → SystemAssembly.

입력 (모두 명시적. 코드 안에 숨은 기본 수치를 두지 않는다)
    config/geometry.yaml, config/operation.yaml, config/scenarios/<scenario>.yaml 의 conditions.<condition>
    선택: config/assumptions/*.yaml — geometry/operation의 null(미확보) 입력을 [가정]값으로 채운다.
          넘긴 경우에만 쓰고, null인 키만 채울 수 있다(원문 값 덮어쓰기는 오류). 적용 항목은 BuiltSystem에 남는다.

조립 규칙
    계산 영역      conditions.<c>.domain.nodes (노드 ID)
    노드 → 장치    core.channels.*.node_id (200·210·220) → Core 하나 (fluid·solid·power·coupling 주입)
                   heat_exchanger.helium_side.{fsh,ees}_node_id (310·320) → HeatExchanger 하나 (SG-0)
                   flow_network.nodes: pipe → Piping, branch → Plenum, pump → Compressor,
                   valve → Plenum (열린 무손실 밸브 [가정], mlsv: open일 때만)
    유체 연결      geometry.junctions 중 양 끝이 영역 안이고 closed_junctions가 아닌 것.
                   같은 장치 쌍 사이의 병렬 연결(J06–J08)은 StreamJunction 하나, 장치 내부 연결(J15)은 생략.
    개방 영역 경계 inlet_boundary.junction(유량·온도), outlet_boundary.junction(압력) → PortBoundary,
                   p_sys는 FixedSystemPressure (하부 플레넘이 소유하지 않음)
    폐루프 기준압  pressure_closure.location 플레넘이 p_sys를 단일 소유. 정상(mode="steady")에서는 그 플레넘의
                   질량식(``<location>.mass``)을 p_sys − p0 = 0으로 바꾼다 (interface_spec §5, 보고서 D2).
    열경계         노심 배럴(HS-06) → riser, riser → RPV(HS-07, He 간극), RPV → RCCS-0 (ThermalLink 3개)
                   나머지 영역 열구조체(HS-09·10·11·13)는 미구현이라 단열 [가정]으로 기록한다.
    재고량         He 체적을 가진 장치(Core는 fluid, SG는 He 측)를 한 번씩 → system.inventory.total_helium_inventory

시험용 조립(tests/fixtures)과 분리된 생산 진입점이다. 실행 흐름은 system/simulation.py가 묶는다.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from interfaces.ports import FluidStreamPort
from models._inputs import ModelInputMissing
from models.compressor import Compressor
from models.core import Core
from models.core.coupling import UnitCellCoreCoupling
from models.core.fluid import ParallelChannelCoreFluid
from models.core.power import SpecifiedCorePower
from models.core.solid import UnitCellCoreSolid
from models.heat_exchanger import HeatExchanger
from models.heat_exchanger.coupling import SG0Coupling
from models.heat_exchanger.helium_side import SGHeliumSide
from models.heat_exchanger.secondary_side import FixedTemperatureWaterSink
from models.heat_exchanger.wall import NoStorageWall
from models.piping import Piping
from models.plenum import Plenum
from models.rccs.rccs import RCCS0Sink
from models.vessel import Vessel
from models.vessel.vessel import SurfaceFluidConvection, helium_gap_conductance
from physics.heat_transfer.correlations import make_convection_correlation
from physics.heat_transfer.solid_properties import R1SolidProperties
from physics.properties import make_helium_properties
from system.assembly import SystemAssembly
from system.boundaries import FixedSystemPressure, PortBoundary
from system.connections import StreamJunction
from system.inventory import total_helium_inventory

ROOT = Path(__file__).resolve().parents[1]
MODES = ("steady", "transient")
MODELED_HEAT_STRUCTURES = {"HS-01", "HS-02", "HS-03", "HS-04", "HS-05", "HS-06", "HS-08",   # Core solid
                           "HS-07",                                                         # Vessel
                           "HS-12"}                                                         # SG-0 (축열 없는 벽)
PLENUM_INLET = ("plenum", 0.5, 1.0)     # 플레넘 → 유로: 급축소 K 0.5 (역류 급확대 1.0) + 베르누이
PLENUM_OUTLET = ("plenum", 1.0, 0.5)    # 유로 → 플레넘: 급확대 K 1.0 (역류 급축소 0.5) + 베르누이
PIPE_END = ("pipe", 0.0, 0.0)           # 유로 ↔ 유로/장치: 형상손실은 연결부 K (미보정 = 0)


# =====================================================================
# 입력 읽기와 검사
# =====================================================================
@dataclass(frozen=True)
class Assumption:
    """assumptions 파일이 채운 null 입력 한 개."""

    source: str     # assumptions 파일 (ROOT 기준 상대 경로)
    key: str        # 예: geometry.core.channels.bypass.hydraulic_diameter_m
    value: object


@dataclass(frozen=True)
class MissingInput:
    file: str
    key: str
    note: str = ""


class MissingInputsError(ModelInputMissing):
    """실행 전 입력 검사에서 찾은 미정(null) 입력 목록. 파일·키·시나리오·조건을 모두 적는다."""

    def __init__(self, scenario, condition, items):
        self.scenario, self.condition, self.items = scenario, condition, tuple(items)
        lines = [f"scenario '{scenario}' condition '{condition}': {len(self.items)} required input(s) missing"]
        lines += [f"  - {m.file}: {m.key}" + (f"  ({m.note})" if m.note else "") for m in self.items]
        lines.append("  원문 값을 config에 넣거나, [가정]값을 config/assumptions/*.yaml로 명시해 --assumptions로 넘긴다.")
        super().__init__("\n".join(lines))


@dataclass(frozen=True)
class ScenarioInputs:
    scenario_name: str
    condition_name: str
    geometry: dict
    operation: dict
    scenario: dict
    condition: dict
    files: dict                       # 역할 → ROOT 기준 상대 경로 문자열
    assumptions: tuple = ()           # (Assumption, ...)

    @property
    def closed_loop(self) -> bool:
        return "pressure_closure" in self.condition


def _rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _load_yaml(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"config file not found: {path}")
    with path.open(encoding="utf-8") as f:
        return yaml.safe_load(f)


def _fill_nulls(base: dict, override: dict, prefix: str, source: str, applied: list):
    for k, v in override.items():
        path = f"{prefix}.{k}"
        current = base.get(k)
        if isinstance(v, dict):
            if current is None:
                base[k] = current = {}
            elif not isinstance(current, dict):
                raise ValueError(f"{source}: {path} is a mapping but the config value is {current!r}")
            _fill_nulls(current, v, path, source, applied)
        elif current is None:
            base[k] = v
            applied.append(Assumption(source, path, v))
        elif current != v:
            raise ValueError(f"{source}: {path} = {v!r} would replace the sourced value {current!r}; "
                             "assumptions may only fill null inputs")


def load_inputs(scenario="steady_state", condition=None, *, config_dir=None, assumptions=()) -> ScenarioInputs:
    """시나리오·조건과 config를 읽고, 넘긴 assumptions 파일로 null 입력만 채운다."""
    config_dir = Path(config_dir) if config_dir is not None else ROOT / "config"
    paths = {"geometry": config_dir / "geometry.yaml", "operation": config_dir / "operation.yaml",
             "scenario": config_dir / "scenarios" / f"{scenario}.yaml"}
    geometry = copy.deepcopy(_load_yaml(paths["geometry"]))
    operation = copy.deepcopy(_load_yaml(paths["operation"]))
    scen = _load_yaml(paths["scenario"])
    conditions = scen.get("conditions") or {}
    condition = condition or scen.get("default_condition")
    if condition not in conditions:
        raise KeyError(f"{_rel(paths['scenario'])}: condition {condition!r} not in {sorted(conditions)}")
    applied: list[Assumption] = []
    files = {k: _rel(p) for k, p in paths.items()}
    calibrations, fingerprints = [], []
    for i, a in enumerate(assumptions or ()):
        a = Path(a)
        a = a if a.is_absolute() else ROOT / a
        data = _load_yaml(a) or {}
        unknown = set(data) - {"meta", "geometry", "operation"}
        if unknown:
            raise ValueError(f"{_rel(a)}: unsupported top-level keys {sorted(unknown)} (meta, geometry, operation)")
        source = _rel(a)
        _fill_nulls(geometry, data.get("geometry") or {}, "geometry", source, applied)
        _fill_nulls(operation, data.get("operation") or {}, "operation", source, applied)
        meta = data.get("meta") or {}
        if meta.get("kind") == "calibration":           # 보정 결과 파일 (system/calibration.py가 씀)
            extra = set(data.get("geometry") or {}) - {"junctions"}
            extra |= {f"junctions.{j}.{k}" for j, d in ((data.get("geometry") or {}).get("junctions") or {}).items()
                      for k in (d or {}) if k != "K_form"}
            if extra or data.get("operation"):
                raise ValueError(f"{source}: a calibration file may only set geometry.junctions.<J>.K_form "
                                 f"(found {sorted(extra) or ['operation']})")
            files[f"calibration[{i}]"] = source
            calibrations.append((source, tuple(meta.get("depends_on_assumptions") or ())))
            fingerprints.append((source, meta.get("input_fingerprint")))
            basis = (meta.get("basis") or {}).get("condition")
            if basis is not None and basis != condition:
                raise ValueError(f"{source} was calibrated for condition '{basis}', not '{condition}'")
        else:
            files[f"assumptions[{i}]"] = source
    loaded = {v for k, v in files.items() if k.startswith("assumptions[")}
    for source, deps in calibrations:
        missing = [d for d in deps if d not in loaded]
        if missing:
            raise ValueError(f"{source} was calibrated with {missing}; load those assumption files too")
    if calibrations:
        # 보정 K는 그 K를 만든 입력에서만 유효하다. 적용 전에 유효 입력값 지문을 대조한다
        current = input_fingerprint(geometry, operation, conditions[condition])
        for source, expected in fingerprints:
            _check_calibration_fingerprint(source, expected, current)
    return ScenarioInputs(scenario, condition, geometry, operation, scen, copy.deepcopy(conditions[condition]),
                          files, tuple(applied))


# =====================================================================
# 유효 입력값 지문 (K 보정 파일 검증)
# =====================================================================
FINGERPRINT_VERSION = 1
_CONDITION_DOC_KEYS = ("description", "outputs", "comparison", "acceptance")   # 풀이에 쓰지 않는 설명·비교 항목


class CalibrationInputMismatch(ValueError):
    """보정 파일을 만든 입력과 현재 입력이 다르다 (다시 보정해야 함)."""


def _normalize(obj):
    """YAML 파싱값을 지문용으로 정규화: 수는 float, 키는 문자열, 그 밖은 그대로 (bool·None·문자열)."""
    if isinstance(obj, dict):
        return {str(k): _normalize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_normalize(v) for v in obj]
    if isinstance(obj, bool) or obj is None or isinstance(obj, str):
        return obj
    if isinstance(obj, (int, float)):
        return float(obj)
    return str(obj)                      # 날짜 등


def _digest(obj) -> str:
    text = json.dumps(_normalize(obj), sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def input_fingerprint(geometry, operation, condition) -> dict:
    """주석·개행·수 표기에 영향받지 않는 유효 입력값 지문.

    부분: geometry(연결부 K_form과 meta 제외), operation(meta·reference 제외), condition(설명·비교 항목 제외),
    helium_properties(참조 파일의 helium 절), solid_properties(참조 파일 전체). 가정 파일로 채운 값은 이미
    geometry·operation에 들어 있으므로 함께 들어간다. 부분별 SHA-256과 전체 SHA-256을 돌려준다.
    """
    geo = copy.deepcopy(geometry)
    geo.pop("meta", None)
    for j in (geo.get("junctions") or {}).values():
        if isinstance(j, dict):
            j.pop("K_form", None)
    op = {k: v for k, v in operation.items() if k not in ("meta", "reference")}
    # 물성 파일은 경로가 아니라 내용으로 지문에 넣는다 (아래 helium/solid_properties)
    cond = {k: v for k, v in condition.items()
            if k not in _CONDITION_DOC_KEYS and k not in ("helium_properties", "solid_properties")}
    he_ref = str(condition.get("helium_properties") or "")
    he_file, _, section = he_ref.partition("#")
    he_data = _load_yaml(ROOT / he_file) if he_file else {}
    he_part = (he_data or {}).get(section or "helium") if isinstance(he_data, dict) else None
    solid_file = condition.get("solid_properties") or he_file
    solid_part = _load_yaml(ROOT / solid_file) if solid_file else None
    parts = {"geometry": _digest(geo), "operation": _digest(op), "condition": _digest(cond),
             "helium_properties": _digest(he_part), "solid_properties": _digest(solid_part)}
    return {"version": FINGERPRINT_VERSION, "sha256": _digest(parts), "parts": parts}


def _check_calibration_fingerprint(source, expected, current):
    hint = "입력이 바뀌었으면 `python main.py --calibrate`로 보정 파일을 다시 만든다."
    if not isinstance(expected, dict) or "sha256" not in expected:
        raise CalibrationInputMismatch(f"{source}: meta.input_fingerprint가 없다 (이전 형식). {hint}")
    if expected.get("version") != current["version"]:
        raise CalibrationInputMismatch(f"{source}: 지문 형식 {expected.get('version')} ≠ {current['version']}. {hint}")
    if expected["sha256"] != current["sha256"]:
        differ = sorted(k for k, v in current["parts"].items() if (expected.get("parts") or {}).get(k) != v)
        raise CalibrationInputMismatch(f"{source}: 보정 K를 만든 입력과 현재 입력이 다르다 (다른 부분: "
                                       f"{', '.join(differ) or '알 수 없음'}). {hint}")


def fingerprint_of(inputs: "ScenarioInputs") -> dict:
    return input_fingerprint(inputs.geometry, inputs.operation, inputs.condition)


def _get(mapping, dotted):
    node = mapping
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


# =====================================================================
# 영역 → 장치 단위 (Unit)
# =====================================================================
@dataclass
class _Unit:
    key: str                    # Assembly Component 이름
    kind: str                   # pipe | branch | pump | valve | core | sg
    node_ids: tuple
    config: dict | None = None
    inlet: str | None = None    # 들어오는 연결 이름 (StreamJunction 또는 경계)
    outlet: str | None = None


@dataclass(frozen=True)
class _Connection:
    name: str
    upstream: str | None        # unit key (영역 밖이면 None)
    downstream: str | None
    junctions: tuple


def _units(inputs: ScenarioInputs):
    geo, cond = inputs.geometry, inputs.condition
    gfile = inputs.files["geometry"]
    domain = [int(n) for n in _get(cond, "domain.nodes") or ()]
    if not domain:
        raise ModelInputMissing(f"{inputs.files['scenario']}: conditions.{inputs.condition_name}.domain.nodes")
    core_ids = {int(ch["node_id"]): name for name, ch in (_get(geo, "core.channels") or {}).items()
                if isinstance(ch, dict) and ch.get("node_id") is not None}
    hs = _get(geo, "heat_exchanger.helium_side") or {}
    sg_ids = {int(hs[k]) for k in ("fsh_node_id", "ees_node_id") if hs.get(k) is not None}
    node_cfg = {int(c["node_id"]): (name, c) for name, c in (_get(geo, "flow_network.nodes") or {}).items()
                if isinstance(c, dict) and c.get("node_id") is not None}
    unit_of, units = {}, {}
    for nid in domain:
        if nid in core_ids:
            key, kind, cfg = "core", "core", None
        elif nid in sg_ids:
            key, kind, cfg = "sg", "sg", None
        elif nid in node_cfg:
            key, cfg = node_cfg[nid]
            kind = cfg.get("type")
            if kind not in ("pipe", "branch", "pump", "valve"):
                raise ValueError(f"{gfile}: flow_network.nodes.{key}.type {kind!r} is not supported by the builder")
        else:
            raise ValueError(f"{inputs.files['scenario']}: domain node {nid} is not defined in {gfile}")
        unit_of[nid] = key
        if key in units:
            units[key].node_ids += (nid,)
        else:
            units[key] = _Unit(key, kind, (nid,), cfg)
    for group, ids in (("core", set(core_ids)), ("sg", sg_ids)):
        present = ids & set(domain)
        if present and present != ids:
            raise ValueError(f"{inputs.files['scenario']}: domain must contain all {group} nodes {sorted(ids)} "
                             f"or none (got {sorted(present)})")
    return units, unit_of


def _connections(inputs: ScenarioInputs, units, unit_of):
    """영역 안 연결부를 장치 쌍으로 묶고, 잘린 연결부는 경계로 표시한다."""
    cond = inputs.condition
    closed = set(_get(cond, "domain.closed_junctions") or ())
    boundary_j = {}
    for side in ("inlet_boundary", "outlet_boundary"):
        j = _get(cond, f"{side}.junction")
        if j is not None:
            boundary_j[j] = side
    grouped: dict[tuple, list] = {}
    for jname, j in (inputs.geometry.get("junctions") or {}).items():
        if jname in closed:
            continue
        a, b = unit_of.get(int(j["from"])), unit_of.get(int(j["to"]))
        if a is None and b is None:
            continue
        if a is not None and a == b:
            continue                                      # 장치 내부 (예: J15 FSH → EES)
        if a is None or b is None:
            if jname not in boundary_j:
                raise ValueError(f"{inputs.files['scenario']}: junction {jname} ({j['from']} → {j['to']}) "
                                 f"crosses the domain of '{inputs.condition_name}' but is neither closed nor a boundary")
        grouped.setdefault((a, b), []).append(jname)
    conns = []
    for (a, b), js in grouped.items():
        name = js[0] if len(js) == 1 else f"{js[0]}_{js[-1]}"
        conns.append(_Connection(name, a, b, tuple(js)))
        for key, attr in ((a, "outlet"), (b, "inlet")):
            if key is None:
                continue
            unit = units[key]
            if getattr(unit, attr) is not None:
                raise NotImplementedError(f"node group '{key}' has more than one {attr} connection "
                                          f"({getattr(unit, attr)}, {name}); multi-port branches need port elevations")
            setattr(unit, attr, name)
    for unit in units.values():
        if unit.inlet is None or unit.outlet is None:
            raise ValueError(f"{inputs.files['scenario']}: node group '{unit.key}' in '{inputs.condition_name}' "
                             f"has no {'inlet' if unit.inlet is None else 'outlet'} connection")
    missing_bc = set(boundary_j) - {c.junctions[0] for c in conns if c.upstream is None or c.downstream is None}
    if missing_bc:
        raise ValueError(f"{inputs.files['scenario']}: boundary junction(s) {sorted(missing_bc)} do not cut the domain")
    return conns


FLOW_KINDS = ("pipe", "core", "sg")      # 형상손실 K를 실을 면이 있는 장치


def core_channel_nodes(geometry) -> dict:
    """노심 유로 노드 ID → 유로 이름 (average/hot/bypass)."""
    return {int(ch["node_id"]): name for name, ch in (_get(geometry, "core.channels") or {}).items()
            if isinstance(ch, dict) and ch.get("node_id") is not None}


def junction_face(inputs, units, unit_of, jname):
    """연결부 J의 형상손실 K를 실을 면: (장치, "inlet"|"outlet", 노심 유로|None) 또는 None(영역 밖).

    K_form은 그 면의 기본 경계 손실(플레넘 경계 급축소 0.5·급확대 1.0, 그 밖 0)에 더하는 증분이다.
    음수도 되지만(기본 날카로운 모서리보다 매끄러운 입구) 면의 정·역 K 합은 0 이상이어야 한다.

    하류 장치가 유로(pipe/core/sg)면 그 입구 면, 아니면 상류 장치의 출구 면. 양쪽 다 유로가 아니면(플레넘–순환기) 오류.
    장치 내부 연결(J15)은 지원하지 않는다.
    """
    j = inputs.geometry["junctions"][jname]
    a, b = unit_of.get(int(j["from"])), unit_of.get(int(j["to"]))
    if a is not None and a == b:
        raise NotImplementedError(f"junction {jname} is internal to '{a}'; internal face K is not supported")
    core_nodes = core_channel_nodes(inputs.geometry)
    if b is not None and units[b].kind in FLOW_KINDS:
        return (b, "inlet", core_nodes.get(int(j["to"])) if b == "core" else None)
    if a is not None and units[a].kind in FLOW_KINDS:
        return (a, "outlet", core_nodes.get(int(j["from"])) if a == "core" else None)
    if a is None or b is None:
        return None                                    # 면이 영역 밖 장치에 있음 → 이 조건에서는 쓰이지 않음
    raise ValueError(f"{inputs.files['geometry']}: junction {jname} ({j['from']} → {j['to']}) has no pipe/core/SG "
                     "face to carry K_form")


def _junction_k(inputs, units, unit_of):
    """영역 안 연결부의 K_form 합을 면별로 모은다. {(장치, 면, 유로): K}, 적용한 연결부 목록."""
    closed = set(_get(inputs.condition, "domain.closed_junctions") or ())
    faces, used = {}, []
    for jname, j in (inputs.geometry.get("junctions") or {}).items():
        K = j.get("K_form") if isinstance(j, dict) else None
        if K is None or jname in closed:
            continue
        K = float(K)
        if not np.isfinite(K):
            raise ValueError(f"junctions.{jname}.K_form must be finite (got {K})")
        face = junction_face(inputs, units, unit_of, jname)
        if face is None:
            continue
        faces[face] = faces.get(face, 0.0) + K
        used.append(jname)
    return faces, tuple(used)


def _with_k(boundary, K):
    """경계 면 (종류, K_정, K_역)에 연결부 K_form 증분을 정·역 대칭으로 더한다. 합이 음수면 오류."""
    kind, kf, kr = boundary
    if kf + K < 0 or kr + K < 0:
        raise ValueError(f"junction K_form {K:+.4g} makes the face loss negative (base {kf:g}/{kr:g})")
    return (kind, kf + K, kr + K)


def _flow_order(units, conns):
    """기준 장치(폐루프: p_sys 소유 플레넘, 개방: 입구 경계)부터 유동 방향 순서."""
    by_inlet = {c.name: c for c in conns}
    start = next((c.downstream for c in conns if c.upstream is None), None)
    order, key, seen = [], start, set()
    if key is None:
        key = next(iter(units))
    while key is not None and key not in seen:
        seen.add(key)
        order.append(key)
        key = by_inlet[units[key].outlet].downstream
    if len(order) != len(units):
        raise ValueError(f"flow path does not visit every node group: {sorted(set(units) - set(order))}")
    return order


# =====================================================================
# 필수 입력 검사
# =====================================================================
def check_inputs(inputs: ScenarioInputs) -> list:
    """조건 실행에 필요한 입력 중 null/누락인 것을 모두 모은다 (생성 전에 한 번에 보고)."""
    geo, op, cond = inputs.geometry, inputs.operation, inputs.condition
    gf, of, sf = inputs.files["geometry"], inputs.files["operation"], inputs.files["scenario"]
    cpath = f"conditions.{inputs.condition_name}"
    missing = []

    def need(mapping, dotted, file, prefix="", note=""):
        if _get(mapping, dotted) is None:
            missing.append(MissingInput(file, f"{prefix}{dotted}", note))

    units, _ = _units(inputs)
    for unit in units.values():
        base = f"flow_network.nodes.{unit.key}."
        if unit.kind == "pipe":
            for k in ("flow_area_m2", "hydraulic_diameter_m", "length_m", "z_inlet_m", "z_outlet_m"):
                need(unit.config, k, gf, base, "M11" if "area" in k or "diameter" in k else "")
        elif unit.kind in ("branch", "valve"):
            for k in ("volume_m3", "z_inlet_m", "z_outlet_m"):
                need(unit.config, k, gf, base, "M10" if k == "volume_m3" else "")
        elif unit.kind == "pump":
            need(unit.config, "z_inlet_m", gf, base)
        elif unit.kind == "core":
            for name, ch in (_get(geo, "core.channels") or {}).items():
                if isinstance(ch, dict) and "node_id" in ch:
                    need(ch, "flow_area_m2", gf, f"core.channels.{name}.")
                    if ch.get("cell_hydraulic_diameter_m") is None:
                        need(ch, "hydraulic_diameter_m", gf, f"core.channels.{name}.")
            need(geo, "core.elevation_profile_m", gf)
            need(op, "reactor.power_distribution", of)
        elif unit.kind == "sg":
            for k in ("elevation_profile_m", "hydraulic_diameter_m"):
                need(geo, f"heat_exchanger.helium_side.{k}", gf, note="M11" if "diameter" in k else "")
            need(geo, "heat_exchanger.heat_transfer_area_m2", gf)
            need(op, "boundaries.sg_heat_sink.water_temperature_K", of)
            need(cond, "sg_heat_sink.calibration_target_helium_outlet_K", sf, f"{cpath}.")
    need(cond, "reactor.thermal_power_W", sf, f"{cpath}.")
    need(cond, "rccs.boundary_temperature_K", sf, f"{cpath}.")
    need(op, "helium_inventory.stagnant_volume_temperature_K", of, note="정체 공간 He 재고 추정")
    for k in ("rccs.boundary.radius_m", "rccs.boundary.emissivity", "core.radial_m.core_barrel_outer",
              "core.radial_m.core_barrel_inner", "core.radial_m.riser_outer", "core.radial_m.rpv_inner"):
        need(geo, k, gf)
    if inputs.closed_loop:
        need(cond, "pressure_closure.pressure_Pa", sf, f"{cpath}.")
        need(cond, "circulator.mass_flow_kg_s", sf, f"{cpath}.")
    else:
        for k in ("inlet_boundary.mass_flow_kg_s", "inlet_boundary.temperature_K", "outlet_boundary.pressure_Pa"):
            need(cond, k, sf, f"{cpath}.")
    return missing


# =====================================================================
# 조립 결과
# =====================================================================
@dataclass
class BuiltSystem:
    system: SystemAssembly
    inputs: ScenarioInputs
    mode: str
    parts: dict                                  # 이름 → Component/하위 모델
    inventory: tuple                             # He 체적 장치 (한 번씩)
    flow_order: tuple                            # 유동 순서 장치 이름
    connections: tuple                           # (이름, 상류, 하류, 원래 연결부들)
    model_assumptions: tuple = ()                # builder 수준 모델 가정 (문장)
    p_sys_owner: str = ""
    replaced_residual_id: str | None = None      # 정상 폐합에서 바뀐 질량식 ID
    assumptions_used: tuple = ()                 # 이 영역에 해당하는 Assumption (나머지는 영역 밖이라 미사용)
    k_form: dict = field(default_factory=dict)   # 적용한 연결부 K_form 증분

    def total_inventory(self, t, y) -> float:
        return total_helium_inventory(self.inventory, t, y, self.system.context)

    def describe(self) -> dict:
        cal_files = {v for k, v in self.inputs.files.items() if k.startswith("calibration[")}
        return {"scenario": self.inputs.scenario_name, "condition": self.inputs.condition_name, "mode": self.mode,
                "unknowns": self.system.registry.size, "flow_order": list(self.flow_order),
                "connections": [dict(name=c[0], upstream=c[1], downstream=c[2], junctions=list(c[3]))
                                for c in self.connections],
                "p_sys_owner": self.p_sys_owner, "replaced_residual_id": self.replaced_residual_id,
                "input_files": dict(self.inputs.files),
                "assumptions": [dict(source=a.source, key=a.key, value=a.value) for a in self.assumptions_used
                                if a.source not in cal_files],
                "calibration": [dict(source=a.source, key=a.key, value=a.value) for a in self.assumptions_used
                                if a.source in cal_files],
                "assumptions_outside_domain": [a.key for a in self.inputs.assumptions if a not in self.assumptions_used],
                "junction_K_form": dict(self.k_form),
                "model_assumptions": list(self.model_assumptions)}


def _assumption_in_domain(a: Assumption, units) -> bool:
    """가정 키가 이 계산 영역의 장치에 해당하는지 (영역 밖 노드·SG 값은 쓰이지 않는다)."""
    parts = a.key.split(".")
    if parts[:3] == ["geometry", "flow_network", "nodes"] and len(parts) > 3:
        return parts[3] in units
    if parts[:2] == ["geometry", "heat_exchanger"] or parts[:3] == ["operation", "boundaries", "sg_heat_sink"]:
        return "sg" in units
    if parts[:2] == ["geometry", "core"]:
        return "core" in units
    return True


def _properties(inputs):
    cond, sf = inputs.condition, inputs.files["scenario"]
    he_ref = cond.get("helium_properties")
    if not he_ref:
        raise ModelInputMissing(f"{sf}: conditions.{inputs.condition_name}.helium_properties")
    he_file, _, section = str(he_ref).partition("#")
    if section not in ("", "helium"):
        raise ValueError(f"{sf}: helium_properties section must be '#helium' (got {he_ref!r})")
    he = make_helium_properties("ideal_gas", path=ROOT / he_file)
    solid_file = cond.get("solid_properties") or he_file
    subcase = cond.get("default_subcase")
    if subcase is not None:
        props = _get(cond, f"benchmark_subcases.{subcase}.properties")
        if props != "fixed":
            raise NotImplementedError(f"{sf}: subcase {subcase} needs {props!r} solid properties; only R1 Ex. 2b "
                                      "fixed properties are implemented (fluence map missing)")
    return he, R1SolidProperties(ROOT / solid_file, mode="ex2b_fixed")


def _station_temperatures(op):
    return {int(s["node"]): float(s["T_K"]) for s in (_get(op, "initial_distribution.stations") or ())
            if s.get("node") is not None and s.get("T_K") is not None}


def _require_mode_values(cond, key, expected, sf, cpath):
    value = _get(cond, key)
    if value not in expected:
        raise NotImplementedError(f"{sf}: {cpath}.{key} = {value!r}; builder supports {expected}")


def build_system(inputs: ScenarioInputs, *, mode="steady", U_sg=None, power=None) -> BuiltSystem:
    """``inputs``의 조건을 Component로 조립한다.

    mode="steady"     정상 초기화용: 폐루프 기준압 폐합(질량식 1개 교체), SG-0 U 보정 미지수
    mode="transient"  과도용: 원래 질량식, 고정 U (``U_sg``: 정상해의 보정값)
    power             총 출력 [W] 상수 또는 t → W 함수. None이면 조건의 reactor.thermal_power_W
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    missing = check_inputs(inputs)
    if missing:
        raise MissingInputsError(inputs.scenario_name, inputs.condition_name, missing)
    geo, op, cond = inputs.geometry, inputs.operation, inputs.condition
    sf, cpath = inputs.files["scenario"], f"conditions.{inputs.condition_name}"
    steady = mode == "steady"
    closed = inputs.closed_loop
    units, unit_of = _units(inputs)
    conns = _connections(inputs, units, unit_of)
    he, sp = _properties(inputs)
    corr = make_convection_correlation("mceligot")
    cp = float(he.cp_const)
    model_assumptions = []

    # ---- 운전점과 초기 추정값 (operation.initial_distribution) ------------------------------
    Q = float(cond["reactor"]["thermal_power_W"]) if power is None else power
    Q0 = float(Q(0.0)) if callable(Q) else float(Q)
    if closed:
        _require_mode_values(cond, "pressure_closure.type", ("pressure",), sf, cpath)
        _require_mode_values(cond, "circulator.boundary_type", ("prescribed_mass_flow",), sf, cpath)
        m = float(cond["circulator"]["mass_flow_kg_s"])
        p0 = float(cond["pressure_closure"]["pressure_Pa"])
    else:
        m = float(cond["inlet_boundary"]["mass_flow_kg_s"])
        p0 = float(cond["outlet_boundary"]["pressure_Pa"])
    order = _flow_order(units, conns)
    stations = _station_temperatures(op)
    if not closed:
        stations[units[order[0]].node_ids[0]] = float(cond["inlet_boundary"]["temperature_K"])
    T_guess, T_prev = {}, None
    for key in order:
        T = stations.get(units[key].node_ids[0], T_prev)
        if T is None:
            raise ModelInputMissing(f"{inputs.files['operation']}: initial_distribution.stations has no temperature "
                                    f"for the first node group '{key}' of {cpath}")
        if units[key].kind == "core":
            T = (T, T + Q0 / (m * cp))
            T_prev = T[1]
        else:
            T_prev = T
        T_guess[key] = T
    for i, key in enumerate(order):                     # SG: 입·출구 추정 평균
        if units[key].kind == "sg":
            T_out = T_guess[order[(i + 1) % len(order)]]
            T_guess[key] = 0.5 * (T_guess[key] + (T_out[0] if isinstance(T_out, tuple) else T_out))
    T_cold = min(t[0] if isinstance(t, tuple) else t for t in T_guess.values())
    ini = lambda key: dict(m_dot=m, T=T_guess[key], p=p0)  # noqa: E731

    k_faces, k_junctions = _junction_k(inputs, units, unit_of)

    def end(key, side):
        """배관 끝 경계: 이웃이 플레넘이면 플레넘 경계식(급축소 0.5/급확대 1.0), 아니면 0. 연결부 K_form을 더한다."""
        conn = next(c for c in conns if c.name == (units[key].inlet if side == "inlet" else units[key].outlet))
        other = conn.upstream if side == "inlet" else conn.downstream
        if other is not None and units[other].kind in ("branch", "valve"):
            base = PLENUM_INLET if side == "inlet" else PLENUM_OUTLET
        else:
            base = PIPE_END
        return _with_k(base, k_faces.get((key, side, None), 0.0))

    # ---- 장치 생성 ---------------------------------------------------------------------
    owner = cond["pressure_closure"]["location"] if closed else None
    replaced = None
    if closed:
        if owner not in units or units[owner].kind != "branch":
            raise ValueError(f"{sf}: {cpath}.pressure_closure.location {owner!r} must be a plenum in the domain")
        if steady:
            replaced = f"{owner}.mass"
            declared = _get(cond, "pressure_closure.replaced_residual_id")
            if declared is not None and declared != replaced:
                raise ValueError(f"{sf}: {cpath}.pressure_closure.replaced_residual_id {declared!r} "
                                 f"but the p_sys-owning plenum replaces {replaced!r}")
    riser_node = _get(geo, "heat_structures.HS-06.outer")
    riser_key = unit_of.get(int(riser_node)) if riser_node is not None else None
    if riser_key is None or units[riser_key].kind != "pipe":
        raise ValueError(f"{inputs.files['geometry']}: heat_structures.HS-06.outer must be a pipe node in the domain")
    comps, parts, inventory = {}, {"he": he, "sp": sp}, []
    U_value = None
    for key in order:
        unit = units[key]
        if unit.kind in ("branch", "valve"):
            if unit.kind == "valve":
                _require_mode_values(cond, "mlsv", ("open",), sf, cpath)
                model_assumptions.append(f"{key} (밸브 {unit.node_ids[0]}): 열린 무손실 플레넘 — 체크밸브 모델 미구현")
            comp = Plenum(unit.config, he, name=key, initial=ini(key), owns_p_sys=key == owner,
                          steady_pressure_closure=p0 if (key == owner and steady) else None,
                          port_labels=("inlet", "outlet"))
            inventory.append(comp)
        elif unit.kind == "pipe":
            heat = None
            if key == riser_key:
                def heat(t, y, context, pipe):
                    return context.heat_rates["barrel_riser.to_b"] + context.heat_rates["riser_rpv.to_b"]
            comp = Piping(unit.config, he, name=key, initial=ini(key), inlet=end(key, "inlet"),
                          outlet=end(key, "outlet"), heat_source=heat)
            inventory.append(comp)
        elif unit.kind == "pump":
            comp = Compressor(op["compressor"], he, name=key, mode="prescribed_flow", m_dot_set=m, m_rated=m,
                              initial=ini(key))
        elif unit.kind == "core":
            T_in, T_out = T_guess[key]
            power_cfg = op["reactor"]["power_distribution"]
            declared = _get(cond, "reactor.power_distribution")
            if declared is not None and declared != power_cfg.get("type"):
                raise ValueError(f"{sf}: {cpath}.reactor.power_distribution {declared!r} differs from "
                                 f"{inputs.files['operation']} reactor.power_distribution.type {power_cfg.get('type')!r}")
            n_z = len(geo["core"]["elevation_profile_m"]) - 1
            core_face_K = {}
            for (ukey, side, label), K in k_faces.items():
                if ukey == "core":
                    arr = core_face_K.setdefault(label, np.zeros(n_z + 1))
                    arr[0 if side == "inlet" else -1] += K
            fluid = ParallelChannelCoreFluid(geo, he, initial=dict(m_dot=m, T_inlet=T_in, T_outlet=T_out, p_inlet=p0),
                                             mode_bypass=cond["bypass_model"], face_K=core_face_K or None)
            solid = UnitCellCoreSolid(geo, sp, helium=he, initial=dict(T=np.linspace(T_in, T_out, n_z)),
                                      outer_link="barrel_riser")
            pw = SpecifiedCorePower(solid, power_cfg, Q)
            coupling = UnitCellCoreCoupling(fluid, solid, corr,
                                            bypass_regions=geo["core"]["channels"]["bypass"]["type_I_fraction"])
            comp = Core(fluid=fluid, solid=solid, power=pw, coupling=coupling)
            parts.update(fluid=fluid, solid=solid, power=pw, coupling=coupling)
            inventory.append(fluid)
        elif unit.kind == "sg":
            _require_mode_values(cond, "sg_heat_sink.model_level", ("SG-0",), sf, cpath)
            sink_cfg = op["boundaries"]["sg_heat_sink"]
            T_w = _get(cond, "sg_heat_sink.water_temperature_K")
            if T_w is not None and float(T_w) != float(sink_cfg["water_temperature_K"]):
                raise ValueError(f"{sf}: {cpath}.sg_heat_sink.water_temperature_K differs from "
                                 f"{inputs.files['operation']} boundaries.sg_heat_sink")
            sg_he = SGHeliumSide(geo, he, initial=ini(key), name="sg.helium", inlet=end(key, "inlet"),
                                 outlet=end(key, "outlet"))
            if steady:
                sink = FixedTemperatureWaterSink(
                    sink_cfg, helium_side=sg_he,
                    calibrate_outlet_K=float(cond["sg_heat_sink"]["calibration_target_helium_outlet_K"]))
            else:
                if U_sg is None:
                    raise ValueError("mode='transient' with SG-0 needs U_sg (the steady-state calibrated U)")
                U_value = float(U_sg)
                sink = FixedTemperatureWaterSink(sink_cfg, helium_side=sg_he, U=U_value)
            comp = HeatExchanger(helium_side=sg_he, wall=NoStorageWall(n_cells=sg_he.n), secondary_side=sink,
                                 coupling=SG0Coupling(sg_he, sink), name="sg")
            parts.update(sg_he=sg_he, sink=sink)
            inventory.append(sg_he)
        comps[key] = comp
    parts.update(comps)

    # ---- 열경계: 배럴 → riser → RPV → RCCS-0 --------------------------------------------
    _require_mode_values(cond, "rccs.model_level", ("RCCS-0",), sf, cpath)
    riser, solid = comps[riser_key], parts["solid"]
    rpv = Vessel(geo["rpv"], sp, name="rpv", z_faces=riser.z_faces, initial=dict(T=T_cold - 50.0),
                 links=(("riser_rpv", "to_a"), ("rpv_rccs", "to_a")))
    rad = geo["core"]["radial_m"]
    r_bo, r_bi = float(rad["core_barrel_outer"]), float(rad["core_barrel_inner"])
    r_ro, r_vi = float(rad["riser_outer"]), float(rad["rpv_inner"])
    z_core = np.asarray(geo["core"]["elevation_profile_m"], dtype=float)
    dz_core = np.abs(np.diff(z_core))
    barrel_half = 2 * math.pi * sp.conductivity("alloy_800H") * dz_core / math.log(r_bo / (0.5 * (r_bi + r_bo)))
    barrel_key = solid._state_key(("core_barrel",))
    rccs = RCCS0Sink(geo["rccs"], rpv, sp, T_sink=float(cond["rccs"]["boundary_temperature_K"]))
    links = [
        SurfaceFluidConvection("barrel_riser", solid_T=lambda y, c: c.get(y, barrel_key), z_solid=z_core,
                               area=2 * math.pi * r_bo * dz_core, pipe=riser, correlation=corr,
                               series=lambda Ts, Tf, p: barrel_half).link(),
        SurfaceFluidConvection("riser_rpv", solid_T=lambda y, c: c.get(y, rpv.T_key), z_solid=rpv.z_faces,
                               area=2 * math.pi * r_ro * rpv.dz, pipe=riser, correlation=corr,
                               series=helium_gap_conductance(he, r_ro, r_vi, rpv.dz, sp.emissivity("alloy_800H"),
                                                             sp.emissivity("pressure_vessel"))).link(),
        rccs.link(),
    ]
    parts.update(rpv=rpv, rccs=rccs)
    unmodeled = [hs for hs in _get(cond, "domain.heat_structures") or () if hs not in MODELED_HEAT_STRUCTURES]
    if unmodeled:
        model_assumptions.append(f"열구조체 {', '.join(unmodeled)}: 미구현 → 단열")

    # ---- 유체 연결과 경계 -----------------------------------------------------------------
    def port(key, side):
        comp = comps[key]
        if units[key].kind in ("branch", "valve"):
            return comp.port(side)
        if units[key].kind == "core":
            return parts["fluid"].ports()[f"helium_{side}"]
        return comp.ports()[f"helium_{side}"]

    junctions, boundaries = [], []
    if not closed:
        boundaries.append(FixedSystemPressure(p0))
        model_assumptions.append("개방 영역 p_sys: 출구 압력으로 이완 고정 (τ = 1 s)")
    for c in conns:
        if c.upstream is not None and c.downstream is not None:
            junctions.append(StreamJunction(c.name, (port(c.upstream, "outlet"), port(c.downstream, "inlet")),
                                            mass_flow_epsilon=1e-2, pressure_scale=1e3, flow_scale=m,
                                            enthalpy_scale=5e5))
        elif c.upstream is None:                        # 입구 경계: 유량 + 유입 온도
            T_in = float(cond["inlet_boundary"]["temperature_K"])
            boundaries.append(PortBoundary(f"bc_{c.name}", port(c.downstream, "inlet"), m_dot=m,
                                           h_in=he.enthalpy(p0, T_in)))
        else:                                           # 출구 경계: 압력 (+ 역류 시 유입 온도 = 추정값)
            T_up = T_guess[c.upstream]
            T_up = T_up[1] if isinstance(T_up, tuple) else T_up
            boundaries.append(PortBoundary(f"bc_{c.name}", port(c.upstream, "outlet"), p=p0,
                                           h_in=he.enthalpy(p0, T_up)))
    if k_junctions:
        sources = sorted({a.source for a in inputs.assumptions if a.key.startswith("geometry.junctions.")})
        model_assumptions.append(f"연결부 K_form {len(k_junctions)}개 적용 ({', '.join(k_junctions)}; 출처 "
                                 f"{', '.join(sources) or inputs.files['geometry']}) — 정·역 대칭 [가정], "
                                 "플레넘–유로 경계 급축소 0.5 / 급확대 1.0에 더함")
    else:
        model_assumptions.append("연결부 K_form 미보정 (0). 플레넘–유로 경계는 급축소 0.5 / 급확대 1.0")

    system = SystemAssembly(boundaries + [comps[k] for k in order] + [rpv], connections=junctions,
                            thermal_links=links)
    built = BuiltSystem(system=system, inputs=inputs, mode=mode, parts=parts, inventory=tuple(inventory),
                        flow_order=tuple(order),
                        connections=tuple((c.name, c.upstream, c.downstream, c.junctions) for c in conns),
                        model_assumptions=tuple(model_assumptions),
                        p_sys_owner=owner if closed else boundaries[0].name, replaced_residual_id=replaced,
                        assumptions_used=tuple(a for a in inputs.assumptions
                                               if (a.key.split(".")[2] in k_junctions
                                                   if a.key.startswith("geometry.junctions.")
                                                   else _assumption_in_domain(a, units))),
                        k_form={j: float(geo["junctions"][j]["K_form"]) for j in k_junctions})
    if U_value is not None:
        parts["U_sg"] = U_value
    validate_built(built, junctions + boundaries)
    return built


# =====================================================================
# 조립 후 검사
# =====================================================================
def validate_built(built: BuiltSystem, connection_objects=()):
    """p_sys 단일 소유, 정상 폐합 ID, 미연결 유체 포트를 검사한다 (SystemAssembly 검사에 더해)."""
    system = built.system
    reg = system.registry
    p_names = [n for n in reg.names if n == "p_sys"]
    if len(p_names) != 1:
        raise ValueError(f"p_sys must be registered exactly once (got {len(p_names)})")
    spec = reg.spec("p_sys")
    if spec.kind != "differential" or spec.owner != "primary_loop":
        raise ValueError(f"p_sys must be differential with owner 'primary_loop' (got {spec.kind}, {spec.owner})")
    ids = system.residual_ids()
    closures = [i for i in ids if i.endswith(".steady_pressure_closure")]
    expected = 1 if built.replaced_residual_id else 0
    if len(closures) != expected:
        raise ValueError(f"expected {expected} steady pressure closure residual(s), got {closures}")
    if built.replaced_residual_id and built.replaced_residual_id in ids:
        raise ValueError(f"{built.replaced_residual_id} must be replaced in the steady system")
    if built.p_sys_owner in built.parts and not built.replaced_residual_id and f"{built.p_sys_owner}.mass" not in ids:
        raise ValueError(f"transient system must keep {built.p_sys_owner}.mass")
    used = {}
    for obj in connection_objects:
        if hasattr(obj, "connected_ports"):
            ports = obj.connected_ports()
        elif isinstance(obj, StreamJunction):
            ports = obj.ports
        else:
            continue
        for p in ports:
            k = reg.resolve(p.m_dot)
            if k in used:
                raise ValueError(f"fluid port {p.m_dot} is connected twice")
            used[k] = obj.name
    for comp in system.components:
        for pname, p in comp.ports().items():
            if isinstance(p, FluidStreamPort) and reg.resolve(p.m_dot) not in used:
                raise ValueError(f"fluid port {comp.name}.{pname} ({p.m_dot}) is not connected")
