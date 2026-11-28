"""전역 보존 진단 (Role 5).

``balance_residual``은 단순 산술 유틸리티다. 아래 함수들은 builder가 만든 실제 계통(BuiltSystem)에서
전역 질량·에너지를 진단한다. 진단값이며 과도 잔차식에 들어가지 않는다 (국부 질량·에너지식이 해를 정한다).

정의 (저마하수, 운동·위치에너지 무시, h = c_p (T − T_ref))
    M_he_loop    = Σ_i ρ(p_sys, h_i) V_i                         He 체적마다 한 번 (builder.inventory)
    U_he_loop    = Σ_i M_i h_i − p_sys Σ_i V_i                    He 내부에너지 [J] (T_ref 기준)
    U_solid_total= Σ_j ∫_{T_ref}^{T_j} C_j(T) dT                  노심 고체 + RPV [J]. 흑연 c_p(T)가 온도에 따라 크게
                                                                 변하므로 C·T가 아니라 적분 (Gauss–Legendre 8점, T_ref 300 K)
전역 에너지 수지 (내부 열교환 — 고체↔He, 배럴↔riser, riser↔RPV — 은 서로 상쇄되어 나타나지 않는다)
    d(U_he + U_solid)/dt = P_core + Q_circulator − Q_SG − Q_RCCS + Σ_bc ṁ h        (개방 영역만 경계 엔탈피 유속)
    순환기 입열 Q_circulator는 축일이 He에 준 에너지 전체(모델에서 He 엔탈피 상승)이고, 압력일 V ṗ_sys는
    U_he 정의(−p_sys V)에 들어 있다. SG-0 물과 RCCS-0 냉각면은 계통 밖 열침원이다.
"""
from __future__ import annotations

import numpy as np

from models._inputs import ModelInputMissing
from system.boundaries import PortBoundary


def balance_residual(accumulation_rate, source_rate, sink_rate):
    """d(저장량)/dt - 유입원 + 유출원. 단위 일치가 선행되어야 함."""
    return accumulation_rate - source_rate + sink_rate


# ---------------------------------------------------------------------------------------------
# 저장량
# ---------------------------------------------------------------------------------------------
def _enthalpy_keys(component):
    """He 체적 장치의 셀 엔탈피 상태 키 (helium_volumes와 같은 순서)."""
    from models.core.fluid import ParallelChannelCoreFluid
    from models.heat_exchanger.helium_side import SGHeliumSide

    if isinstance(component, ParallelChannelCoreFluid):
        return [component._ck(ch, "h_he") for ch in component.channels]
    if isinstance(component, SGHeliumSide):
        return [component.pipe.key("h_he")]
    return [component.key("h_he")]


def helium_mass(built, t, y) -> float:
    """M_he_loop [kg]: 모델 유로 체적의 He. 과도 질량 보존 진단은 이 값만 쓴다."""
    return built.total_inventory(t, y)


def stagnant_helium_mass(built, t, y) -> dict:
    """정체 공간 He 추정 [kg] (비교용, 상태 없음).

    체적은 geometry.yaml ``flow_network.stagnant_volumes``, 온도는 operation.yaml
    ``helium_inventory.stagnant_volume_temperature_K`` [가정], 압력은 현재 p_sys다.
    """
    inputs = built.inputs
    volumes = (inputs.geometry.get("flow_network") or {}).get("stagnant_volumes") or {}
    T = (inputs.operation.get("helium_inventory") or {}).get("stagnant_volume_temperature_K")
    if T is None:
        raise ModelInputMissing(f"{inputs.files['operation']}: missing helium_inventory.stagnant_volume_temperature_K "
                                f"(scenario '{inputs.scenario_name}', condition '{inputs.condition_name}')")
    V = 0.0
    for name, v in volumes.items():
        if v.get("volume_m3") is None:
            raise ModelInputMissing(f"{inputs.files['geometry']}: missing flow_network.stagnant_volumes.{name}.volume_m3")
        V += float(v["volume_m3"])
    he = built.parts["he"]
    p_sys = float(built.system.context.get(y, "p_sys"))
    rho = float(he.density(p_sys, he.enthalpy(p_sys, float(T))))
    return {"volume_m3": V, "temperature_K": float(T), "density_kg_m3": rho, "mass_kg": rho * V}


def helium_energy(built, t, y) -> float:
    """U_he_loop = Σ M h − p_sys V [J]."""
    ctx = built.system.context
    p_sys = float(ctx.get(y, "p_sys"))
    total = 0.0
    for comp in built.inventory:
        masses = np.asarray(comp.helium_inventory(t, y, ctx), dtype=float)
        h = np.concatenate([np.atleast_1d(np.asarray(ctx.get(y, k), dtype=float)) for k in _enthalpy_keys(comp)])
        if h.shape != masses.shape:
            raise ValueError(f"{comp.name}: enthalpy cells {h.shape} do not match inventory cells {masses.shape}")
        total += float(masses @ h) - p_sys * float(np.sum(comp.helium_volumes()))
    return total


SOLID_T_REF_K = 300.0
_GL_X, _GL_W = np.polynomial.legendre.leggauss(8)


def _sensible_energy(heat_capacity, T):
    """Σ_j ∫_{T_ref}^{T_j} C_j(T') dT' — 노드별 온도 의존 열용량 C(T) [J/K] 배열 함수의 적분."""
    T = np.asarray(T, dtype=float)
    half = 0.5 * (T - SOLID_T_REF_K)
    total = np.zeros_like(T)
    for x, w in zip(_GL_X, _GL_W):
        total += w * np.asarray(heat_capacity(SOLID_T_REF_K + half * (1.0 + x)), dtype=float)
    return float(np.sum(half * total))


def solid_energy(built, t, y) -> float:
    """U_solid_total = Σ ∫ C(T) dT [J] (노심 고체 + RPV, 기준 T_ref = 300 K)."""
    ctx = built.system.context
    solid, rpv = built.parts["solid"], built.parts["rpv"]
    return (_sensible_energy(solid.heat_capacity, solid.temperatures(y, ctx))
            + _sensible_energy(rpv.heat_capacity, np.asarray(ctx.get(y, rpv.T_key), dtype=float)))


# ---------------------------------------------------------------------------------------------
# 경계 열률과 유속
# ---------------------------------------------------------------------------------------------
def boundary_rates(built, t, y) -> dict:
    """전역 경계 항 [W]: 노심 출력, 순환기 입열, SG 제거열, RCCS 제거열, 개방 경계 엔탈피 유속."""
    system = built.system
    out = system.outputs(t, y)
    rates = system._trial_context(t, y).heat_rates
    ctx = system.context
    power = built.parts["power"]
    P = float(power.total_power(t) if callable(power.total_power) else power.total_power)
    terms = {"core_power_W": P,
             "circulator_heat_W": sum(float(v) for k, v in out.items() if k.endswith(".Q_to_helium")),
             "sg_removed_W": float(out.get("sg.secondary.Q_removed", 0.0)),
             "rccs_removed_W": -float(np.sum(rates["rpv_rccs.to_a"]))}
    flux = 0.0
    for comp in system.components:
        if isinstance(comp, PortBoundary):
            m = float(ctx.get(y, comp.port.m_dot))
            h = float(ctx.get(y, comp.port.h_instream if m >= 0 else comp.port.h_outflow))
            flux += m * h
    terms["boundary_enthalpy_flux_W"] = flux
    terms["net_W"] = (terms["core_power_W"] + terms["circulator_heat_W"] - terms["sg_removed_W"]
                      - terms["rccs_removed_W"] + flux)
    return terms


def boundary_mass_flux(built, t, y) -> float:
    """개방 영역 경계 순유입 [kg/s] (폐루프는 0)."""
    ctx = built.system.context
    return float(sum(float(ctx.get(y, comp.port.m_dot)) for comp in built.system.components
                     if isinstance(comp, PortBoundary)))


def steady_energy_balance(built, t, y) -> dict:
    """정상해 전역 에너지 수지: 순 열률과 노심 출력 대비 비."""
    terms = boundary_rates(built, t, y)
    terms["relative_imbalance"] = abs(terms["net_W"]) / max(abs(terms["core_power_W"]), 1.0)
    return terms


# ---------------------------------------------------------------------------------------------
# 과도 이력 수지
# ---------------------------------------------------------------------------------------------
def _trapezoid(values, t):
    values, t = np.asarray(values, dtype=float), np.asarray(t, dtype=float)
    return np.concatenate([[0.0], np.cumsum(0.5 * (values[1:] + values[:-1]) * np.diff(t))])


def transient_balance(built, t, ys) -> dict:
    """적분 이력의 전역 질량·에너지 수지.

    질량   M(t) − M(0) − ∫ Σ ṁ_bc dt
    에너지 [U_he + U_solid](t) − (0) − ∫ (P + Q_circ − Q_SG − Q_RCCS + Σ ṁ h) dt   (사다리꼴 적분)
    오차는 시간 적분·사다리꼴 근사·매끄러운 upwind 근사의 합이다. 노심 출력 × 시간 구간으로 정규화한다.
    출력 이력의 꺾임(램프 끝)을 걸친 스텝에서는 사다리꼴 오차가 커졌다가 다음 스텝에서 상쇄되므로 끝값도 함께 본다.
    """
    t = np.asarray(t, dtype=float)
    M = np.array([helium_mass(built, ti, yi) for ti, yi in zip(t, ys)])
    m_bc = np.array([boundary_mass_flux(built, ti, yi) for ti, yi in zip(t, ys)])
    E = np.array([helium_energy(built, ti, yi) + solid_energy(built, ti, yi) for ti, yi in zip(t, ys)])
    net = np.array([boundary_rates(built, ti, yi)["net_W"] for ti, yi in zip(t, ys)])
    P = np.array([boundary_rates(built, ti, yi)["core_power_W"] for ti, yi in zip(t, ys)])
    mass_err = M - M[0] - _trapezoid(m_bc, t)
    energy_err = E - E[0] - _trapezoid(net, t)
    span = max(float(t[-1] - t[0]), 1e-30)
    e_ref = float(np.max(np.abs(P))) * span
    return {"M_he_kg": M, "U_total_J": E,
            "mass_error_rel_max": float(np.max(np.abs(mass_err)) / M[0]),
            "energy_error_rel_max": float(np.max(np.abs(energy_err)) / e_ref),
            "energy_error_rel_end": float(abs(energy_err[-1]) / e_ref),
            "mass_error_kg": mass_err, "energy_error_J": energy_err}
