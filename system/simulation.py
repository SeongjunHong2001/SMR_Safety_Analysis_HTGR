"""Role 5 실행기: 입력 검사 → 정상해 → 과도 잔차 재확인 → 일치 초기화 → DASSL 무외란 유지 → 보존 진단 → 출력.

척도 구분
    상태 척도   registry.scale_vector()   — 유한차분 Jacobian 단계, 상태 변화량 비교
    잔차 척도   system.residual_scales()  — Newton 수렴 판정과 잔차 노름 (scaled ∞-노름)
풀이 설정은 시나리오 파일의 ``solution`` 절에서 읽는다. config/solver.yaml은 3변수 데모 전용이라 쓰지 않는다.

정상 초기화 (steady_state.yaml solution, interface_spec §5)
    1. builder mode="steady": 폐루프는 p_sys 소유 플레넘의 질량식 1개를 p_sys − p0 = 0으로, SG-0 U는 보정 미지수
    2. ydot = 0 Newton (희소 유한차분 Jacobian: F_y와 F_ydot 비영 패턴의 합집합 + 대각)
    3. 검증: builder mode="transient"(원래 질량식, 보정 U 고정)에 정상해를 이름별로 옮겨 ydot = 0 잔차를 다시 계산
    4. 과도 계에서 일치 초기화 후 DASSL로 무외란 유지 적분 — 정상해가 유지되는지와 전역 질량·에너지 수지 확인
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from models._inputs import ModelInputMissing
from solvers import DASSL
from solvers.initialization import consistent_initialize
from solvers.jacobian import finite_difference_jacobian
from solvers.newton import NewtonOptions, newton_solve
from system.builder import ROOT, BuiltSystem, build_system, load_inputs
from validation import conservation

SG_U_KEY = "sg.secondary.U_overall"


# =====================================================================
# 설정
# =====================================================================
@dataclass(frozen=True)
class HoldSettings:
    t_end_s: float
    dt_initial_s: float
    rtol: float
    atol: float
    max_step_s: float
    max_order: int = 5


@dataclass(frozen=True)
class RunSettings:
    newton: NewtonOptions
    sparse_jacobian: bool
    steady_residual_max: float          # 정상해·검증 잔차 허용 (scaled ∞-노름)
    hold: HoldSettings | None


def settings_from_scenario(inputs) -> RunSettings:
    """시나리오 ``solution`` 절 → RunSettings. 필요한 키가 없으면 파일·키를 적어 실패한다."""
    sol = inputs.scenario.get("solution") or {}
    sf = inputs.files["scenario"]

    def req(mapping, key, where):
        if not isinstance(mapping, dict) or mapping.get(key) is None:
            raise ModelInputMissing(f"{sf}: missing solution.{where}{key} (scenario '{inputs.scenario_name}')")
        return mapping[key]

    if req(sol, "method", "") != "solve_steady_state":
        raise NotImplementedError(f"{sf}: solution.method {sol['method']!r} (only solve_steady_state)")
    nw = req(sol, "newton", "")
    newton = NewtonOptions(max_iter=int(req(nw, "max_iter", "newton.")), rtol=float(req(nw, "rtol", "newton.")),
                           atol=float(req(nw, "atol", "newton.")))
    acc = req(sol, "acceptance", "")
    hold = None
    h = sol.get("hold_transient")
    if h is not None:
        hold = HoldSettings(*(float(req(h, k, "hold_transient.")) for k in
                              ("t_end_s", "dt_initial_s", "rtol", "atol", "max_step_s")),
                            max_order=int(h.get("max_order", 5)))
    return RunSettings(newton=newton, sparse_jacobian=req(sol, "jacobian", "") == "sparse_finite_difference",
                       steady_residual_max=float(req(acc, "scaled_residual_max", "acceptance.")), hold=hold)


# =====================================================================
# 공통 도구
# =====================================================================
def jacobian_pattern(system, t, y, ydot):
    """F_y와 F_ydot 비영 패턴의 합집합 + 대각 (bool 밀집 행렬).

    y와 ydot을 같은 열에서 서로 다른 배율(1, 0.37)로 함께 섭동해 유한차분 한 번으로 합집합을 얻는다.
    """
    y, ydot = np.asarray(y, dtype=float), np.asarray(ydot, dtype=float)
    scale = system.registry.scale_vector()

    def g(z):
        d = (z - y / scale) * scale
        return system.residual(t, y + d, ydot + 0.37 * d)

    J = finite_difference_jacobian(g, y / scale, step_scales=1.0)
    return (J != 0) | np.eye(y.size, dtype=bool)


def scaled_norm(system, t, y, ydot):
    r, ids, scales = system.residual_details(t, y, ydot)
    z = np.abs(r / scales)
    i = int(np.argmax(z))
    return float(z[i]), ids[i]


def transfer_state(src: BuiltSystem, y_src, dst: BuiltSystem):
    """이름별로 상태를 옮긴다. dst의 모든 상태가 src에 있어야 한다 (src에만 있는 보정 미지수는 버린다)."""
    rs, rd = src.system.registry, dst.system.registry
    missing = [n for n in rd.names if n not in rs.names]
    if missing:
        raise ValueError(f"states missing in the source system: {missing[:5]}")
    y = rd.initial_state()
    for name in rd.names:
        y[rd.spec(name).indices] = np.asarray(rs.get(y_src, name), dtype=float).ravel()
    return y


# =====================================================================
# 단계
# =====================================================================
@dataclass
class SteadyResult:
    built: BuiltSystem
    y: np.ndarray
    residual_norm: float
    worst_residual_id: str
    newton_iterations: int
    wall_s: float


def solve_steady(built: BuiltSystem, settings: RunSettings) -> SteadyResult:
    system = built.system
    t0 = time.perf_counter()
    y0 = system.registry.initial_state()
    zero = np.zeros_like(y0)
    pattern = jacobian_pattern(system, 0.0, y0, zero) if settings.sparse_jacobian else None
    sol = newton_solve(lambda y: system.residual(0.0, y, zero), y0, options=settings.newton,
                       residual_scale=system.residual_scales(), jacobian_sparsity=pattern)
    norm, worst = scaled_norm(system, 0.0, sol.x, zero)
    if norm > settings.steady_residual_max:
        raise RuntimeError(f"steady state not converged: scaled residual {norm:.3e} at {worst} "
                           f"(> {settings.steady_residual_max:g})")
    return SteadyResult(built, sol.x, norm, worst, sol.iterations, time.perf_counter() - t0)


@dataclass
class Verification:
    built: BuiltSystem                 # 과도 계 (원래 질량식, 보정 U 고정)
    y: np.ndarray
    residual_norm: float
    worst_residual_id: str
    U_sg: float | None


def verify_steady(steady: SteadyResult, settings: RunSettings, *, power=None) -> Verification:
    """정상해를 교체 전 과도 잔차에 ydot = 0으로 넣어 다시 확인한다 (steady_state.yaml solution.verification)."""
    reg = steady.built.system.registry
    U = float(reg.get(steady.y, SG_U_KEY)) if SG_U_KEY in reg.names else None
    bt = build_system(steady.built.inputs, mode="transient", U_sg=U, power=power)
    y = transfer_state(steady.built, steady.y, bt)
    norm, worst = scaled_norm(bt.system, 0.0, y, np.zeros_like(y))
    if norm > settings.steady_residual_max:
        raise RuntimeError(f"steady solution fails the transient residual: {norm:.3e} at {worst}")
    return Verification(bt, y, norm, worst, U)


@dataclass
class HoldResult:
    result: object                     # solvers.IntegrationResult
    ydot0_differential_max: float      # 일치 초기화 후 미분 상태 변화율 최대 (scaled, 1/s)
    drift_scaled_max: float            # 유지 적분 끝 상태 변화 (상태 척도) 최대
    balance: dict
    wall_s: float


def run_hold(ver: Verification, hold: HoldSettings, settings: RunSettings) -> HoldResult:
    """검증된 정상해에서 무외란 과도를 적분한다."""
    system = ver.built.system
    reg = system.registry
    t0 = time.perf_counter()
    mask = reg.differential_mask()
    scales = system.residual_scales()
    pattern = (jacobian_pattern(system, 0.0, ver.y, np.zeros_like(ver.y)) if settings.sparse_jacobian else None)
    # consistent_initialize의 Newton 미지수는 z = [y(대수), ẏ(미분)] 순서라 패턴 열도 그 순서로 바꾼다
    z_order = np.concatenate([np.flatnonzero(~mask), np.flatnonzero(mask)])
    y0, yp0 = consistent_initialize(system.residual, 0.0, ver.y, mask, residual_scale=scales,
                                    jacobian_sparsity=None if pattern is None else pattern[:, z_order])
    s = reg.scale_vector()
    ydot_max = float(np.max(np.abs(yp0[mask] / s[mask]))) if mask.any() else 0.0
    solver = DASSL(dt=hold.dt_initial_s, max_order=hold.max_order, rtol=hold.rtol, atol=hold.atol,
                   max_step=hold.max_step_s, residual_scale=scales)
    result = solver.solve(system.residual, t_span=(0.0, hold.t_end_s), y0=y0, ydot0=yp0, differential_mask=mask,
                          jacobian_sparsity=pattern,
                          preflight=lambda t, y, yd, cj: system.validate_for_solve(t, y, yd, cj=cj,
                                                                                   jacobian_sparsity=pattern))
    drift = float(np.max(np.abs(result.y[-1] - result.y[0]) / s))
    balance = conservation.transient_balance(ver.built, result.t, result.y)
    return HoldResult(result, ydot_max, drift, balance, time.perf_counter() - t0)


# =====================================================================
# 결과 요약과 비교
# =====================================================================
def _port_p(built, key, side):
    comp = built.parts[key]
    port = comp.port(side) if hasattr(comp, "port") else comp.ports()[f"helium_{side}"]
    return port.p


def summary(built: BuiltSystem, t, y) -> dict:
    """시나리오 outputs·comparison에 대응하는 주요 결과 (SI)."""
    system, reg = built.system, built.system.registry
    out = system.outputs(t, y)
    get = lambda k: float(reg.get(y, k))  # noqa: E731
    s = {"p_sys_Pa": get("p_sys"),
         "fuel_temperature_max_K": float(out["core.solid.T_fuel_max"]),
         "core_outlet_temperature_K": float(out["core.fluid.T_outlet_mixed"]),
         "rpv_temperature_max_K": float(out["rpv.T_max"]),
         "helium_inventory_kg": conservation.helium_mass(built, t, y)}
    bal = conservation.boundary_rates(built, t, y)
    s["rccs_heat_removal_W"] = bal["rccs_removed_W"]
    order = built.flow_order
    up = next(k for k in order if order[(order.index(k) + 1) % len(order)] == "core")
    down = order[(order.index("core") + 1) % len(order)]
    s["core_inlet_temperature_K"] = float(out[f"{up}.T_he"]) if f"{up}.T_he" in out else None
    s["pressure_drop_plenum_to_plenum_Pa"] = get(f"{up}.p_fluid_cell") - get(f"{down}.p_fluid_cell")
    s["core_pressure_drop_Pa"] = s["pressure_drop_plenum_to_plenum_Pa"]
    for ch in built.parts["fluid"].channels:
        s[f"channel_mass_flow_kg_s.{ch.label}"] = float(out[f"core.fluid.{ch.label}.m_dot_inlet"])
    dz = np.abs(np.diff(np.asarray(built.inputs.geometry["core"]["elevation_profile_m"], dtype=float)))
    T_barrel = np.asarray(out["core.solid.core_barrel.T"], dtype=float)
    s["core_barrel_temperature_avg_K"] = float(T_barrel @ dz / dz.sum()) if T_barrel.shape == dz.shape else None
    if "sg" in built.parts:
        s["sg_heat_removal_W"] = float(out["sg.secondary.Q_removed"])
        s["sg_overall_U_W_m2K"] = float(out["sg.secondary.U_overall"])
        s["sg_helium_outlet_temperature_K"] = float(np.asarray(out["sg.helium.T_he"])[-1])
    if "main_circulator" in built.parts:
        s["circulator_pressure_rise_Pa"] = float(out["main_circulator.pressure_rise"])
        s["circulator_heat_W"] = bal["circulator_heat_W"]
        s["circulator_inlet_pressure_Pa"] = get(_port_p(built, "main_circulator", "inlet"))
    if "sg_inlet_plenum" in built.parts:
        s["sg_inlet_pressure_Pa"] = get("sg_inlet_plenum.p_fluid_cell")
    return s


def compare(built: BuiltSystem, values: dict) -> dict:
    """조건의 comparison(참고값)과 acceptance(판정)를 계산값과 나란히 둔다."""
    cond = built.inputs.condition
    rows = {}
    for k, ref in (cond.get("comparison") or {}).items():
        if isinstance(ref, (int, float)) and k in values and values[k] is not None:
            rows[k] = {"reference": float(ref), "computed": values[k], "difference": values[k] - float(ref)}
    accept = {}
    for k, limit in (cond.get("acceptance") or {}).items():
        base = k.removesuffix("_max_K") + "_max_K" if k.endswith("_max_K") else k
        if base in values:
            accept[k] = {"limit": float(limit), "computed": values[base], "pass": bool(values[base] <= float(limit))}
    return {"comparison": rows, "acceptance": accept}


# =====================================================================
# 출력
# =====================================================================
def flatten_names(registry):
    """Registry 이름/shape → CSV 열 이름과 단위. 평탄화한 이름이 겹치면 오류."""
    cols, units = [], []
    for name in registry.names:
        spec = registry.spec(name)
        if spec.shape:
            for idx in np.ndindex(*spec.shape):
                cols.append(f"{name}[{','.join(map(str, idx))}]")
                units.append(spec.unit)
        else:
            cols.append(name)
            units.append(spec.unit)
    if len(set(cols)) != len(cols):
        dup = sorted({c for c in cols if cols.count(c) > 1})
        raise ValueError(f"duplicate flattened output names: {dup[:5]}")
    return cols, units


def _flatten_outputs(out: dict):
    flat = {}
    for k, v in out.items():
        a = np.asarray(v, dtype=float)
        if a.ndim == 0:
            flat[k] = float(a)
        else:
            for idx in np.ndindex(*a.shape):
                flat[f"{k}[{','.join(map(str, idx))}]"] = float(a[idx])
    return flat


def _output_unit(name: str):
    leaf = name.split("[")[0].rsplit(".", 1)[-1]
    for prefix, unit in (("T_", "K"), ("Q_", "W"), ("m_dot", "kg/s"), ("U_overall", "W/(m2 K)"),
                         ("pressure_rise", "Pa"), ("p", "Pa")):
        if leaf.startswith(prefix):
            return unit
    return None


def _code_version():
    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    commit = git("rev-parse", "HEAD")
    return {"commit": commit or None, "branch": git("rev-parse", "--abbrev-ref", "HEAD") or None,
            "dirty": bool(git("status", "--porcelain", "--untracked-files=no")) if commit else None}


def _sha256(rel):
    p = ROOT / rel
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def _jsonable(v):
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    return v


def write_outputs(report, out_dir: Path) -> dict:
    """상태·출력 시계열 CSV 두 개와 메타데이터 JSON을 쓴다."""
    out_dir.mkdir(parents=True, exist_ok=True)
    built = report.transient.built if report.transient is not None else report.steady.built
    reg, system = built.system.registry, built.system
    if report.hold is not None:
        ts, ys = report.hold.result.t, report.hold.result.y
    else:
        ts, ys = np.array([0.0]), np.array([report.transient.y if report.transient else report.steady.y])
    base = f"{built.inputs.scenario_name}_{built.inputs.condition_name}"
    cols, units = flatten_names(reg)
    paths = {"states": out_dir / f"{base}_states.csv", "outputs": out_dir / f"{base}_outputs.csv",
             "meta": out_dir / f"{base}_meta.json"}
    with paths["states"].open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_s", *cols])
        for t, y in zip(ts, ys):
            w.writerow([repr(float(t)), *(repr(float(v)) for v in y)])
    flat_rows = [_flatten_outputs(system.outputs(t, y)) for t, y in zip(ts, ys)]
    out_cols = list(flat_rows[0])
    with paths["outputs"].open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_s", *out_cols])
        for t, row in zip(ts, flat_rows):
            w.writerow([repr(float(t)), *(repr(row[c]) for c in out_cols)])
    meta = report.metadata()
    meta["files"] = {k: p.name for k, p in paths.items()}
    meta["units"] = {"time_s": "s", **dict(zip(cols, units))}
    meta["output_units_inferred"] = {c: _output_unit(c) for c in out_cols}
    paths["meta"].write_text(json.dumps(_jsonable(meta), ensure_ascii=False, indent=2), encoding="utf-8")
    return paths


# =====================================================================
# 묶음 실행
# =====================================================================
@dataclass
class RunReport:
    steady: SteadyResult
    transient: Verification | None
    hold: HoldResult | None
    settings: RunSettings
    summary: dict = field(default_factory=dict)
    checks: dict = field(default_factory=dict)
    energy_balance: dict = field(default_factory=dict)
    started_utc: str = ""
    paths: dict = field(default_factory=dict)

    def metadata(self) -> dict:
        b = self.steady.built
        relevant = b.describe()
        hold = None
        if self.hold is not None:
            r = self.hold.result
            hold = {"t_end_s": float(r.t[-1]), "accepted_steps": int(len(r.orders)), "rejected_steps": int(r.rejected_steps),
                    "error_test_failures": int(r.error_test_failures), "newton_failures": int(r.newton_failures),
                    "max_order": int(max(r.orders)) if len(r.orders) else 0,
                    "newton_iterations_total": int(np.sum(r.newton_iterations)),
                    "ydot0_differential_scaled_max": self.hold.ydot0_differential_max,
                    "drift_scaled_max": self.hold.drift_scaled_max,
                    "mass_error_rel_max": self.hold.balance["mass_error_rel_max"],
                    "energy_error_rel_max": self.hold.balance["energy_error_rel_max"],
                    "wall_s": self.hold.wall_s, "method": r.method}
        return {"kind": "HTGR primary system run (S0 steady state)", "started_utc": self.started_utc,
                "code_version": _code_version(),
                "input_sha256": {k: _sha256(v) for k, v in b.inputs.files.items()},
                **relevant,
                "unknowns": {"steady": b.system.registry.size,
                             "transient": self.transient.built.system.registry.size if self.transient else None},
                "solver": {"steady": {"newton_iterations": self.steady.newton_iterations,
                                      "scaled_residual_max": self.steady.residual_norm,
                                      "worst_residual_id": self.steady.worst_residual_id,
                                      "sparse_jacobian": self.settings.sparse_jacobian, "wall_s": self.steady.wall_s},
                           "verification": None if self.transient is None else
                           {"scaled_residual_max": self.transient.residual_norm,
                            "worst_residual_id": self.transient.worst_residual_id, "U_sg_fixed": self.transient.U_sg},
                           "hold_transient": hold},
                "summary": self.summary, "energy_balance_steady": self.energy_balance, **self.checks,
                "note": "PSID·벤치마크 값은 참고 비교이며 판정은 acceptance만. 가정 목록(assumptions, model_assumptions) 확인"}


def run_scenario(scenario="steady_state", condition=None, *, assumptions=(), config_dir=None, hold=True,
                 output_dir=None, write=True, log=print) -> RunReport:
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    scen_path = (Path(config_dir) if config_dir else ROOT / "config") / "scenarios" / f"{scenario}.yaml"
    head = yaml.safe_load(scen_path.read_text(encoding="utf-8")) if scen_path.is_file() else None
    if head is None:
        raise FileNotFoundError(f"scenario file not found: {scen_path}")
    if head.get("id") != "S0":
        raise NotImplementedError(f"{scen_path.name}: only the S0 steady-state runner is implemented "
                                  f"(scenario id {head.get('id')!r} needs system/scenario.py, STEP C)")
    inputs = load_inputs(scenario, condition, config_dir=config_dir, assumptions=assumptions)
    settings = settings_from_scenario(inputs)
    built = build_system(inputs, mode="steady")
    d = built.describe()
    log(f"[{scenario}/{inputs.condition_name}] built: {built.system.registry.size} unknowns, "
        f"{len(d['assumptions'])} [가정] input(s), {len(d['calibration'])} calibrated K from "
        f"{[v for k, v in inputs.files.items() if k.startswith(('assumptions', 'calibration'))]}")
    steady = solve_steady(built, settings)
    log(f"  steady: {steady.newton_iterations} Newton iterations, scaled residual {steady.residual_norm:.2e}, "
        f"{steady.wall_s:.1f} s")
    ver = verify_steady(steady, settings)
    log(f"  verification (transient residual, ydot=0): {ver.residual_norm:.2e} at {ver.worst_residual_id}")
    report = RunReport(steady, ver, None, settings, started_utc=started)
    if hold and settings.hold is not None:
        report.hold = run_hold(ver, settings.hold, settings)
        r = report.hold.result
        log(f"  hold transient 0–{r.t[-1]:g} s: {len(r.orders)} steps ({r.rejected_steps} rejected), "
            f"drift {report.hold.drift_scaled_max:.2e}, mass err {report.hold.balance['mass_error_rel_max']:.1e}, "
            f"energy err {report.hold.balance['energy_error_rel_max']:.1e}, {report.hold.wall_s:.1f} s")
    report.summary = summary(ver.built, 0.0, ver.y)
    report.energy_balance = conservation.steady_energy_balance(ver.built, 0.0, ver.y)
    report.checks = compare(ver.built, report.summary)
    if write:
        report.paths = write_outputs(report, Path(output_dir) if output_dir else ROOT / "output" / "results")
        log(f"  wrote: {', '.join(str(p) for p in report.paths.values())}")
    return report
