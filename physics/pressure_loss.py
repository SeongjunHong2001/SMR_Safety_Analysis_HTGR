"""1D 헬륨 유로 압력손실 항 (Role 2).

팀 보고서 [2] 「헬륨 유동 및 열전달」 파트 B §3.2(Face 운동량식)와 파트 C §2(마찰·형상손실)를 구현한다.
모든 함수는 Registry를 모르는 순수 계산이고 SI다. 압력 항은 Pa, 정상 유동 방향(면 유량 양수) 기준으로
"하류 압력을 낮추는 방향"이 양수다. Face 운동량 잔차는

    R_mom,j = (p_{i+1} − p_i) + ΔP_acc + ΔP_fric + ΔP_form + ΔP_grav − ΔP_comp   (flow.momentum_residual)

벽면 마찰 (Churchill 1977, 비원형 형상계수 Φs 확장)
    f = 8 [ (8 / (Φs Re))^12 + (A + B)^−1.5 ]^(1/12)
    A = [2.457 ln(1 / ((7/Re)^0.9 + 0.27 ε/D_h))]^16,  B = (37530 / Re)^16
    ΔP_fric = f (L/D_h) ṁ|ṁ| / (2 ρ A²) = (f Re) μ L ṁ / (2 ρ D_h² A)
  f·Re = 8 [ (8/Φs)^12 + Re^12 (A + B)^−1.5 ]^(1/12) 은 Re → 0에서 64/Φs로 유한하므로, |ṁ|가 해석적으로 약분된
  선형 점성식(층류)을 분기 없이 그대로 얻는다 (보고서 C §2.2의 "Re ≤ 100 선형식"과 같은 값, 0 나눗셈 없음).
국부 형상손실 (보고서 C §2.4, B §3.1)
    ΔP_form = K_eff ṁ |ṁ|_sm / (2 ρ A_ref²),  |ṁ|_sm = √(ṁ² + ṁ_ε²),  K_eff = w K⁺ + (1 − w) K⁻,  w = ½[1 + tanh(ṁ/ṁ_ε)]
중력 (반셀 가중, B §3.2)
    ΔP_grav,j = ½ g (ρ_i Δz_elev,i + ρ_{i+1} Δz_elev,{i+1})       Δz_elev = 셀 하나의 (출구 고도 − 입구 고도)
열가속 (B §3.2 주 1)
    ΔP_acc,j = ṁ² / A² (1/ρ_{i+1} − 1/ρ_i)          일정 단면 가열 채널
    ΔP_acc,junc = ½ (ṁ²/(ρ_R A_R²) − ṁ²/(ρ_L A_L²))  단면 급변 접합부 (베르누이 동압 변화)

입력 정책: 기하(길이, 면적, D_h)·상수가 음수·비유한이면 ValueError. 상태 물성(ρ, μ)이 NaN이면 결과도 NaN
(물성 provider의 비물리 상태 정책), 유한한 0 이하이면 ValueError.
"""
from __future__ import annotations

import math

import numpy as np

GRAVITY = 9.80665                 # m/s²
PHI_S_CIRCULAR = 1.0              # 원관
PHI_S_PARALLEL_PLATES = 2.0 / 3.0  # 무한 평행평판 틈새 (f_lam = 96/Re)
SUPPORTED_FRICTION_MODES = ("churchill",)


# ----------------------------------------------------------------- 내부 도구
def _f(x):
    return np.asarray(x, dtype=float)


def _out(x):
    return x[()] if isinstance(x, np.ndarray) and x.ndim == 0 else x


def _geometry(name, value, *, allow_zero=False):
    value = _f(value)
    bad = ~np.isfinite(value) | ((value < 0) if allow_zero else (value <= 0))
    if np.any(bad):
        raise ValueError(f"{name} must be finite and {'>= 0' if allow_zero else '> 0'}; got {value[bad].ravel()[:3]!r}")
    return value


def _state(name, value):
    value = _f(value)
    finite = np.isfinite(value)
    if np.any(finite & (value <= 0)):
        raise ValueError(f"{name} must be > 0; got {value[finite & (value <= 0)].ravel()[:3]!r}")
    return value


def _flow(m_dot):
    return _f(m_dot)


# ----------------------------------------------------------------- 무차원수
def reynolds(m_dot, D_h, A, mu):
    """Re = |ṁ| D_h / (μ A) (항상 0 이상)."""
    D_h, A, mu = _geometry("D_h", D_h), _geometry("A", A), _state("mu", mu)
    with np.errstate(invalid="ignore"):
        return _out(np.abs(_flow(m_dot)) * D_h / (mu * A))


def churchill_f_re(Re, rel_roughness=0.0, phi_s=PHI_S_CIRCULAR):
    """f·Re (Churchill 1977, Φs 확장). Re = 0에서 64/Φs. Re < 0은 ValueError."""
    Re = _f(Re)
    finite = np.isfinite(Re)
    if np.any(finite & (Re < 0)):
        raise ValueError("Reynolds number must be >= 0 (use |ṁ|)")
    eps = _geometry("relative roughness", rel_roughness, allow_zero=True)
    phi_s = _geometry("shape factor phi_s", phi_s)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        inner = (7.0 / Re) ** 0.9 + 0.27 * eps
        A = (2.457 * np.log(1.0 / inner)) ** 16
        B = (37530.0 / Re) ** 16
        turb = Re ** 12 * (A + B) ** -1.5
        turb = np.where(Re == 0, 0.0, turb)          # Re → 0 극한 (A, B → ∞)
        f_re = 8.0 * ((8.0 / phi_s) ** 12 + turb) ** (1.0 / 12.0)
    return _out(np.where(finite, f_re, np.nan))


def darcy_friction_factor(Re, rel_roughness=0.0, phi_s=PHI_S_CIRCULAR):
    """Darcy f (진단용). Re = 0에서는 inf. 압력손실 계산에는 friction_dp를 쓴다."""
    Re = _f(Re)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _out(_f(churchill_f_re(Re, rel_roughness, phi_s)) / Re)


# ----------------------------------------------------------------- 압력 항
def friction_dp(m_dot, rho, mu, length, D_h, A, *, roughness=0.0, phi_s=PHI_S_CIRCULAR,
                mode_friction="churchill"):
    """길이 L 구간의 벽면 마찰 압력손실 [Pa], ṁ과 같은 부호. 층류에서 ṁ에 선형, ṁ = 0에서 매끄럽다."""
    if mode_friction not in SUPPORTED_FRICTION_MODES:
        raise ValueError(f"unsupported mode_friction {mode_friction!r}; supported: {SUPPORTED_FRICTION_MODES}")
    m_dot = _flow(m_dot)
    rho, mu = _state("rho", rho), _state("mu", mu)
    length = _geometry("length", length, allow_zero=True)
    D_h, A = _geometry("D_h", D_h), _geometry("A", A)
    rough = _geometry("roughness", roughness, allow_zero=True)
    with np.errstate(invalid="ignore"):
        Re = np.abs(m_dot) * D_h / (mu * A)
        f_re = _f(churchill_f_re(Re, rough / D_h, phi_s))
        return _out(f_re * mu * length * m_dot / (2.0 * rho * D_h ** 2 * A))


def face_friction_dp(m_dot, left, right, *, roughness=0.0, phi_s=PHI_S_CIRCULAR):
    """Face j = i + ½의 반셀 가중 마찰 [Pa] (보고서 B §3.2).

    left/right: 인접 셀의 dict(rho, mu, dz, D_h, A). 각 셀의 절반 길이(dz/2)에 그 셀의 물성·기하를 쓴다.
    경계 Face처럼 한쪽 셀만 있으면 그쪽에 None을 준다.
    """
    total = 0.0
    for cell in (left, right):
        if cell is None:
            continue
        total = total + _f(friction_dp(m_dot, cell["rho"], cell["mu"], _f(cell["dz"]) / 2.0, cell["D_h"], cell["A"],
                                       roughness=roughness, phi_s=cell.get("phi_s", phi_s)))
    return _out(_f(total))


def smooth_abs(m_dot, m_eps):
    """|ṁ|_sm = √(ṁ² + ṁ_ε²)."""
    m_eps = _geometry("m_eps", m_eps)
    m_dot = _flow(m_dot)
    return _out(np.sqrt(m_dot * m_dot + m_eps * m_eps))


def direction_weight(m_dot, m_eps):
    """w = ½[1 + tanh(ṁ/ṁ_ε)]: 정방향 1, 역류 0, ṁ = 0에서 ½."""
    m_eps = _geometry("m_eps", m_eps)
    return _out(0.5 * (1.0 + np.tanh(_flow(m_dot) / m_eps)))


def form_loss_dp(m_dot, rho, A_ref, K_forward, K_reverse=None, *, m_eps):
    """국부 형상손실 [Pa], ṁ과 같은 부호. K_reverse가 없으면 정·역류 대칭."""
    rho, A_ref = _state("rho", rho), _geometry("A_ref", A_ref)
    K_f = _geometry("K_forward", K_forward, allow_zero=True)
    K_r = K_f if K_reverse is None else _geometry("K_reverse", K_reverse, allow_zero=True)
    m_dot = _flow(m_dot)
    w = _f(direction_weight(m_dot, m_eps))
    K = w * K_f + (1.0 - w) * K_r
    with np.errstate(invalid="ignore"):
        return _out(K * m_dot * _f(smooth_abs(m_dot, m_eps)) / (2.0 * rho * A_ref ** 2))


def sudden_expansion_K(A_small, A_large):
    """급확대 K = (1 − β)², β = A_small/A_large. 기준 면적은 작은 관(상류)."""
    beta = _geometry("A_small", A_small) / _geometry("A_large", A_large)
    if np.any(beta > 1):
        raise ValueError("sudden expansion requires A_small <= A_large")
    return _out((1.0 - beta) ** 2)


def sudden_contraction_K(A_large, A_small):
    """급축소 K = (1/μ_c − 1)², μ_c = 0.62 + 0.38 β³, β = A_small/A_large. 기준 면적은 작은 관(하류).

    NUREG/CR-5535 Vol. 1 Eq. (3.4-42)~(3.4-44) (보고서 C §2.4). β → 0 극한 0.376.
    """
    beta = _geometry("A_small", A_small) / _geometry("A_large", A_large)
    if np.any(beta > 1):
        raise ValueError("sudden contraction requires A_small <= A_large")
    mu_c = 0.62 + 0.38 * beta ** 3
    return _out((1.0 / mu_c - 1.0) ** 2)


def gravity_dp(rho, dz_elev, *, g=GRAVITY):
    """ρ g Δz_elev [Pa]. Δz_elev는 유동 방향으로 오른 고도 (하향 유동이면 음수)."""
    rho = _state("rho", rho)
    with np.errstate(invalid="ignore"):
        return _out(rho * g * _f(dz_elev))


def face_gravity_dp(rho_left, dz_elev_left, rho_right, dz_elev_right, *, g=GRAVITY):
    """Face j의 반셀 가중 중력항 ½ g (ρ_i Δz_i + ρ_{i+1} Δz_{i+1}) [Pa]. 한쪽 셀이 없으면 그쪽 Δz에 0을 준다."""
    return _out(0.5 * (_f(gravity_dp(rho_left, dz_elev_left, g=g)) + _f(gravity_dp(rho_right, dz_elev_right, g=g))))


def acceleration_dp(m_dot, A, rho_left, rho_right):
    """일정 단면 가열 채널의 열가속 ṁ²/A² (1/ρ_R − 1/ρ_L) [Pa] (정·역류 모두 하류 쪽 밀도가 'right'인지 호출자가 정함)."""
    A = _geometry("A", A)
    rho_l, rho_r = _state("rho_left", rho_left), _state("rho_right", rho_right)
    m_dot = _flow(m_dot)
    with np.errstate(invalid="ignore"):
        return _out(m_dot * m_dot / (A * A) * (1.0 / rho_r - 1.0 / rho_l))


def junction_acceleration_dp(m_dot, A_left, A_right, rho_left, rho_right):
    """단면 급변 접합부의 동압 변화 ½ (ṁ²/(ρ_R A_R²) − ṁ²/(ρ_L A_L²)) [Pa]."""
    A_l, A_r = _geometry("A_left", A_left), _geometry("A_right", A_right)
    rho_l, rho_r = _state("rho_left", rho_left), _state("rho_right", rho_right)
    m2 = _flow(m_dot) ** 2
    with np.errstate(invalid="ignore"):
        return _out(0.5 * (m2 / (rho_r * A_r * A_r) - m2 / (rho_l * A_l * A_l)))


def face_acceleration_dp(m_dot, A_left, A_right, rho_left, rho_right):
    """내부 Face의 가속 압력 변화 [Pa]: 열가속(계수 1) + 단면 변화 베르누이(계수 ½).

    ΔP = ṁ²/A_m² (1/ρ_R − 1/ρ_L) + ½ ṁ²/ρ_m (1/A_R² − 1/A_L²),  A_m = ½(A_L + A_R), ρ_m = ½(ρ_L + ρ_R)
    A_L = A_R이면 acceleration_dp(일정 단면, 보고서 B 주 1), ρ_L = ρ_R이면 junction_acceleration_dp와 같다.
    """
    A_l, A_r = _geometry("A_left", A_left), _geometry("A_right", A_right)
    rho_l, rho_r = _state("rho_left", rho_left), _state("rho_right", rho_right)
    m2 = _flow(m_dot) ** 2
    A_m, rho_m = 0.5 * (A_l + A_r), 0.5 * (rho_l + rho_r)
    with np.errstate(invalid="ignore"):
        return _out(m2 / (A_m * A_m) * (1.0 / rho_r - 1.0 / rho_l)
                    + 0.5 * m2 / rho_m * (1.0 / (A_r * A_r) - 1.0 / (A_l * A_l)))


def plenum_boundary_acceleration_dp(m_dot, A_channel, rho_channel, *, side):
    """플레넘(속도 ≈ 0)과 유로 사이 Face의 베르누이 동압 변화 [Pa] (보고서 [2] D §4.3).

    side = "inlet" (플레넘 → 유로 입구 Face): +ṁ²/(2ρA²)
    side = "outlet" (유로 출구 → 플레넘 Face): −ṁ²/(2ρA²)
    유동 방향과 무관하다 (정압 + 동압 보존). 손실은 form_loss_dp의 정·역류 K로 따로 넣는다:
    입구 (1 + K_in) q 강하, 출구 (−1 + K_out) q = 0 (K_out = 1).
    """
    if side not in ("inlet", "outlet"):
        raise ValueError("side must be 'inlet' or 'outlet'")
    A = _geometry("A_channel", A_channel)
    rho = _state("rho_channel", rho_channel)
    with np.errstate(invalid="ignore"):
        q = _flow(m_dot) ** 2 / (2.0 * rho * A * A)
    return _out(q if side == "inlet" else -q)


BOUNDARY_KINDS = ("plenum", "pipe")


def channel_face_terms(m_faces, rho, mu, A, D_h, dz, dz_elev, *, m_eps, roughness=0.0, phi_s=PHI_S_CIRCULAR,
                       inlet=("plenum", 0.5, 1.0), outlet=("plenum", 1.0, 0.5), face_K=None, g=GRAVITY):
    """N셀 1D 유로의 면별 운동량 압력 항 (N+1,) [Pa] (보고서 [2] B §3.2, D §4.3).

    셀 인자 (N,): rho, mu, A, D_h, dz(셀 길이), dz_elev(유동 방향 고도 변화). m_faces (N+1,).
    inlet/outlet = (경계 종류, K_정방향, K_역류). 경계 종류 "plenum"이면 플레넘(속도 0)과 유로 사이 베르누이 ±q를
    넣고, "pipe"면 동압 변화 없이 K 손실만 넣는다. face_K (N+1,)는 위치별 추가 K(정·역 대칭).
    반환 dict(acc, fric, form, grav): flow.momentum_residual의 dp_* 인자로 그대로 쓴다.
    """
    m = _flow(m_faces)
    rho, mu = _f(rho), _f(mu)
    A, D_h, dz, dz_elev = _f(A), _f(D_h), _f(dz), _f(dz_elev)
    n = rho.shape[0]
    if m.shape != (n + 1,):
        raise ValueError(f"m_faces must have N+1 = {n + 1} entries")
    for kind, *_ in (inlet, outlet):
        if kind not in BOUNDARY_KINDS:
            raise ValueError(f"boundary kind must be one of {BOUNDARY_KINDS}")
    # 반셀 마찰: 셀 i의 앞쪽 절반은 면 i(유량 m_i), 뒤쪽 절반은 면 i+1(유량 m_{i+1})
    half = _f(friction_dp(np.concatenate((m[:-1], m[1:])), np.tile(rho, 2), np.tile(mu, 2), np.tile(dz / 2.0, 2),
                          np.tile(D_h, 2), np.tile(A, 2), roughness=roughness, phi_s=phi_s))
    fric = np.zeros(n + 1)
    fric[:-1] += half[:n]
    fric[1:] += half[n:]
    grav_half = 0.5 * _f(gravity_dp(rho, dz_elev, g=g))
    grav = np.zeros(n + 1)
    grav[:-1] += grav_half
    grav[1:] += grav_half
    acc = np.zeros(n + 1)
    if n > 1:
        acc[1:-1] = face_acceleration_dp(m[1:-1], A[:-1], A[1:], rho[:-1], rho[1:])
    if inlet[0] == "plenum":
        acc[0] = plenum_boundary_acceleration_dp(m[0], A[0], rho[0], side="inlet")
    if outlet[0] == "plenum":
        acc[-1] = plenum_boundary_acceleration_dp(m[-1], A[-1], rho[-1], side="outlet")
    rho_face = np.concatenate(([rho[0]], 0.5 * (rho[:-1] + rho[1:]), [rho[-1]]))
    A_face = np.concatenate(([A[0]], np.minimum(A[:-1], A[1:]), [A[-1]]))
    K_f = np.zeros(n + 1)
    K_r = np.zeros(n + 1)
    K_f[0], K_r[0] = float(inlet[1]), float(inlet[2])
    K_f[-1], K_r[-1] = float(outlet[1]), float(outlet[2])
    if face_K is not None:
        extra = _geometry("face_K", face_K, allow_zero=True)
        if extra.shape != (n + 1,):
            raise ValueError(f"face_K must have N+1 = {n + 1} entries")
        K_f, K_r = K_f + extra, K_r + extra
    form = _f(form_loss_dp(m, rho_face, A_face, K_f, K_r, m_eps=m_eps))
    return {"acc": acc, "fric": fric, "form": form, "grav": grav}


def mdot_epsilon(m_nominal_path, eps_rel=1.0e-4):
    """경로별 스무딩 폭 ṁ_ε = eps_mdot_rel × ṁ_nom,path (변수표 eps_mdot_rel 1e-4)."""
    m_nom = abs(float(m_nominal_path))
    if not (math.isfinite(m_nom) and m_nom > 0 and math.isfinite(eps_rel) and eps_rel > 0):
        raise ValueError("m_nominal_path and eps_rel must be positive and finite")
    return eps_rel * m_nom
