"""1D 헬륨 유한체적(FVM) 보존 커널 (Role 2).

팀 보고서 [2] 「헬륨 유동 및 열전달」 파트 B §3(엇갈림 격자·상류차분·DASSL 잔차)와 변수표·방정식표의
저마하수 계약(docs/interface_spec.md §4)을 구현한다. Registry·Component를 모르는 순수 함수다.

격자와 부호
    셀 i = 0..N−1: 엔탈피 h_i, 체적 V_i 등 셀량은 shape (N,)
    면 j = 0..N  : 면 i+½의 질량유량 ṁ_j는 shape (N+1,). 기하학적 축 방향(정격 유동 방향)이 양수.
                    면 0은 장치 입구, 면 N은 장치 출구다.
    포트 부호     : Component 포트 유량은 "장치로 들어오면 양수"다. 입구 포트 = +ṁ_0, 출구 포트 = −ṁ_N
                    (port_inflows_from_faces).
저마하수 근사
    밀도·물성은 전역 p_sys와 셀 h_i로 평가한다 (ρ = ρ(p_sys, h)). 국부 압력 p_i는 운동량식에만 쓴다.

잔차 (보고서 B §3.2)
    질량   R_mass,i   = V_i (ρ_p,i ṗ_sys + ρ_h,i ḣ_i) − (ṁ_{i−½} − ṁ_{i+½})                      [kg/s]
    에너지 R_energy,i = ρ_i V_i ḣ_i − V_i ṗ_sys + ṁ_{i+½}(h*_{i+½} − h_i) − ṁ_{i−½}(h*_{i−½} − h_i) − Q_i   [W]
    운동량 R_mom,j    = δ_inert (l/A)_j dṁ_j/dt + (p_{i+1} − p_i) + ΔP_acc + ΔP_fric + ΔP_form + ΔP_grav − ΔP_comp  [Pa]
    R_energy + h_i R_mass 는 보존형 d(ρVh)/dt − Vṗ + Σ(면 엔탈피 플럭스) − Q 와 항등이다 (T_ref 무관).
면 엔탈피 (C∞ 상류차분)
    h*_j = w h_left + (1 − w) h_right,  w = ½[1 + tanh(ṁ_j/ṁ_ε)],  ṁ_ε = eps_mdot_rel × ṁ_nom,path
운동량 관성
    δ_inert = 0 (기본, 변수표 delta_inert = false): ṁ_j는 대수(algebraic) 상태, Index-1 DAE.
    δ_inert = 1이면 ṁ_j가 미분(differential) 상태가 된다 (face_flow_state_kind).
"""
from __future__ import annotations

import numpy as np

from .pressure_loss import direction_weight


def _f(x):
    return np.asarray(x, dtype=float)


def _out(x):
    return x[()] if isinstance(x, np.ndarray) and x.ndim == 0 else x


# ----------------------------------------------------------------- 면 수송값
def smooth_upwind(m_dot, val_left, val_right, m_eps):
    """면 수송값 w·left + (1 − w)·right. |ṁ| ≫ ṁ_ε에서 1차 상류차분과 같고 ṁ = 0에서 평균."""
    w = _f(direction_weight(m_dot, m_eps))
    return _out(w * _f(val_left) + (1.0 - w) * _f(val_right))


def face_enthalpies(m_faces, h_cells, h_inlet_boundary, h_outlet_boundary, m_eps):
    """면 엔탈피 h*_j, shape (N+1,).

    면 0의 왼쪽은 입구 경계(상류 포트 stream) 엔탈피, 면 N의 오른쪽은 출구 경계(역류 시 들어오는 stream) 엔탈피다.
    """
    m_faces, h_cells = _f(m_faces), _f(h_cells)
    _check_shapes(m_faces, h_cells)
    left = np.concatenate(([float(h_inlet_boundary)], h_cells))
    right = np.concatenate((h_cells, [float(h_outlet_boundary)]))
    return _f(smooth_upwind(m_faces, left, right, m_eps))


def smooth_positive(x, eps):
    """max(x, 0)의 C∞ 근사 ½(x + √(x² + ε²)) (항상 양수). 포트 유입량 가중 등에 쓴다."""
    x = _f(x)
    eps = _f(eps)
    return _out(0.5 * (x + np.sqrt(x * x + eps * eps)))


def port_inflows_from_faces(m_faces):
    """(입구 포트 유입, 출구 포트 유입) = (+ṁ_0, −ṁ_N). 포트는 장치로 들어오는 방향이 양수."""
    m_faces = _f(m_faces)
    return _out(m_faces[0]), _out(-m_faces[-1])


def face_flow_state_kind(delta_inert):
    """면 유량의 StateRegistry kind: 관성 없음(기본) → "algebraic", 관성 포함 → "differential"."""
    if not isinstance(delta_inert, (bool, np.bool_)):
        raise TypeError("delta_inert must be a bool")
    return "differential" if delta_inert else "algebraic"


# ----------------------------------------------------------------- 셀 잔차
def mass_residual(V, drho_dp, drho_dh, p_sys_dot, h_dot, m_in, m_out):
    """R_mass = V (ρ_p ṗ_sys + ρ_h ḣ) − (ṁ_in − ṁ_out) [kg/s]. ṁ_in/out은 축 방향 면 유량(부호 있음)."""
    return _out(_f(V) * (_f(drho_dp) * _f(p_sys_dot) + _f(drho_dh) * _f(h_dot)) - (_f(m_in) - _f(m_out)))


def energy_residual(rho, V, h, h_dot, p_sys_dot, m_in, h_in_star, m_out, h_out_star, q_wall):
    """R_energy = ρVḣ − Vṗ_sys + ṁ_out(h*_out − h) − ṁ_in(h*_in − h) − Q_wall [W]. Q_wall은 유체로 들어오면 양수."""
    rho, V, h = _f(rho), _f(V), _f(h)
    return _out(rho * V * _f(h_dot) - V * _f(p_sys_dot)
                + _f(m_out) * (_f(h_out_star) - h) - _f(m_in) * (_f(h_in_star) - h) - _f(q_wall))


def cell_residuals(*, V, rho, drho_dp, drho_dh, h, h_dot, p_sys_dot, q_wall,
                   m_faces, h_inlet_boundary, h_outlet_boundary, m_eps):
    """한 유로(N셀)의 질량·에너지 잔차를 한 번에 계산한다.

    셀 인자 shape (N,) (스칼라는 broadcast), m_faces shape (N+1,).
    반환: dict(R_mass (N,), R_energy (N,), h_face (N+1,))
    """
    h = _f(h)
    m_faces = _f(m_faces)
    _check_shapes(m_faces, h)
    h_face = face_enthalpies(m_faces, h, h_inlet_boundary, h_outlet_boundary, m_eps)
    m_in, m_out = m_faces[:-1], m_faces[1:]
    R_mass = _f(mass_residual(V, drho_dp, drho_dh, p_sys_dot, h_dot, m_in, m_out))
    R_energy = _f(energy_residual(rho, V, h, h_dot, p_sys_dot, m_in, h_face[:-1], m_out, h_face[1:], q_wall))
    shape = h.shape
    return {"R_mass": np.broadcast_to(R_mass, shape).copy(),
            "R_energy": np.broadcast_to(R_energy, shape).copy(),
            "h_face": h_face}


def energy_flux_balance(m_faces, h_face):
    """셀별 순 엔탈피 유출 Σ(ṁ h*)_out − (ṁ h*)_in [W], shape (N,). 보존 검산용."""
    flux = _f(m_faces) * _f(h_face)
    return flux[1:] - flux[:-1]


# ----------------------------------------------------------------- 면 잔차
def momentum_residual(p_left, p_right, dp_acc=0.0, dp_fric=0.0, dp_form=0.0, dp_grav=0.0, dp_comp=0.0,
                      *, include_inertia=False, l_over_A=0.0, m_dot_dot=0.0):
    """R_mom = δ (l/A) dṁ/dt + (p_right − p_left) + ΔP_acc + ΔP_fric + ΔP_form + ΔP_grav − ΔP_comp [Pa].

    ΔP_comp는 순환기가 있는 면의 승압(양수면 하류 압력을 올림). 관성은 include_inertia=True일 때만 넣는다.
    """
    if not isinstance(include_inertia, (bool, np.bool_)):
        raise TypeError("include_inertia must be a bool")
    r = (_f(p_right) - _f(p_left) + _f(dp_acc) + _f(dp_fric) + _f(dp_form) + _f(dp_grav) - _f(dp_comp))
    if include_inertia:
        r = r + _f(l_over_A) * _f(m_dot_dot)
    return _out(r)


def face_inertia_length(dz_left, A_left, dz_right, A_right):
    """(l/A)_j = Δz_i/(2A_i) + Δz_{i+1}/(2A_{i+1}) [1/m]. 경계면은 없는 쪽 Δz에 0을 준다."""
    return _out(_f(dz_left) / (2.0 * _f(A_left)) + _f(dz_right) / (2.0 * _f(A_right)))


# ----------------------------------------------------------------- 내부
def _check_shapes(m_faces, h_cells):
    if h_cells.ndim != 1 or m_faces.ndim != 1:
        raise ValueError("h_cells must be (N,) and m_faces must be (N+1,) 1-D arrays")
    if m_faces.shape[0] != h_cells.shape[0] + 1:
        raise ValueError(f"face flows must have N+1 = {h_cells.shape[0] + 1} entries, got {m_faces.shape[0]}")
