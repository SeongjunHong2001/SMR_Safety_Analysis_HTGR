"""헬륨 유로 대류 Nusselt 상관식 (Role 3).

팀 보고서 [2] 「헬륨 유동 및 열전달」 파트 C §3의 결정 트리를 구현한다. 변수표 매개변수 이름을 그대로 쓴다.

강제대류 (난류, mode_nu = "mceligot")
    Nu_turb = 0.021 Re_safe^0.8 Pr^0.4 r_eff^-0.5            R1 Table AIV.13 (McEligot)
    r_eff   = ½ (r_T + 1 + √((r_T − 1)² + eps_T_ratio_clamp²)),  r_T = T_wall / T_bulk
              (냉각 구간 T_wall < T_bulk에서 보정이 1로 매끄럽게 수렴)
    Re_safe = √(Re² + Re_min_safe²)                          (Re → 0에서 Re^0.8 도함수 발산 방지, 역류 대칭)
층류 바닥값
    Nu_lam  = 48/11 ≈ 4.364 (원관 등열유속), 평행평판 틈새 8.235
층류–난류 천이 (C∞)
    w_tr    = ½ [1 + tanh((Re_safe − Re_tr_mid) / Re_tr_width)]
    Nu_fc   = (1 − w_tr) Nu_lam + w_tr Nu_turb
자연대류 (include_nc, 수직 유로 Churchill–Chu, 대표길이 D_h)
    Nu_nc   = [0.825 + 0.387 Ra^(1/6) / (1 + (0.492/Pr)^(9/16))^(8/27)]²
    Ra      = g β |ΔT|_sm D_h³ ρ² c_p / (μ k),  β = 1/T_bulk (이상기체),  |ΔT|_sm = √((T_w − T_b)² + dT_eps_nc²)
결합 (smooth max)
    Nu      = (Nu_fc^m + Nu_nc^m)^(1/m),  m = m_smooth_max
진단
    Ri      = Gr / Re_safe²   (0.1 미만 강제, 0.1–10 혼합, 10 초과 자연대류)

물성(μ, k, c_p, ρ)은 벌크 온도에서 평가한 값을 state로 받는다 (Role 2 provider).
주의: 관내 자연순환 유량은 1D 운동량식 중력항이 계산하므로 Nu_nc를 더하면 P-LOFC에서 Nu가 커지는
비보수 편향이 있다. 안전해석에서는 include_nc=False 민감도를 함께 계산한다 (보고서 [2] C §3.2).
"""
from dataclasses import dataclass

import numpy as np

from interfaces.contracts import HeatTransferCorrelation

from ._checks import as_float, bad_temperature, mask, require_geometry
from .convection import coefficient_from_nusselt

GRAVITY = 9.80665                   # m/s²
NU_LAM_CIRCULAR = 48.0 / 11.0       # Sellars et al. (1956) 원관 등열유속 완전발달 층류
NU_LAM_PARALLEL_PLATES = 8.235      # 평행평판 양면 대칭 등열유속
SHAPES = ("circular", "parallel_plates")
SUPPORTED_NU_MODES = ("mceligot",)

_STATE_KEYS = ("m_dot", "T_bulk", "T_wall", "mu", "k", "cp")


def _state_value(state, key):
    try:
        return as_float(state[key])
    except KeyError as exc:
        raise KeyError(f"convection state is missing {key!r} (required: {_STATE_KEYS} + 'rho' when include_nc)") from exc


def _require_positive_state(name, value):
    finite = np.isfinite(value)
    if np.any(finite & (value <= 0)):
        raise ValueError(f"state {name} must be > 0; got {value[finite & (value <= 0)].ravel()[:3]!r}")
    return value


@dataclass(frozen=True)
class McEligotConvection:
    """``HeatTransferCorrelation``: ``coefficient(state, geometry)`` → h_conv [W/(m²·K)].

    state:    m_dot [kg/s, 부호 허용], T_bulk, T_wall [K], mu [Pa·s], k [W/(m·K)], cp [J/(kg·K)],
              rho [kg/m³] (include_nc일 때)
    geometry: D_h [m], A_flow [m²] (m_dot이 지나는 유로면적), shape ("circular" 기본 | "parallel_plates")
    """

    Re_min_safe: float = 100.0
    Re_tr_mid: float = 3150.0
    Re_tr_width: float = 350.0
    eps_T_ratio_clamp: float = 0.01
    include_nc: bool = True
    m_smooth_max: float = 6.0
    dT_eps_nc: float = 0.1
    gravity: float = GRAVITY
    Nu_lam_circular: float = NU_LAM_CIRCULAR
    Nu_lam_parallel_plates: float = NU_LAM_PARALLEL_PLATES

    def __post_init__(self):
        for name in ("Re_min_safe", "Re_tr_mid", "Re_tr_width", "eps_T_ratio_clamp", "m_smooth_max",
                     "dT_eps_nc", "gravity", "Nu_lam_circular", "Nu_lam_parallel_plates"):
            value = getattr(self, name)
            if not (np.isfinite(value) and value > 0):
                raise ValueError(f"McEligotConvection.{name} must be positive and finite, got {value!r}")
        if self.m_smooth_max < 1:
            raise ValueError("m_smooth_max must be >= 1")
        if not isinstance(self.include_nc, (bool, np.bool_)):
            raise TypeError("include_nc must be a bool")

    def _laminar(self, geometry):
        shape = geometry.get("shape", "circular")
        if shape not in SHAPES:
            raise ValueError(f"geometry shape must be one of {SHAPES}, got {shape!r}")
        return self.Nu_lam_circular if shape == "circular" else self.Nu_lam_parallel_plates

    def evaluate(self, state, geometry):
        """모든 중간값을 dict로 돌려준다 (진단·시험용). 비물리 온도 원소는 NaN."""
        m_dot = _state_value(state, "m_dot")
        T_b = _state_value(state, "T_bulk")
        T_w = _state_value(state, "T_wall")
        mu = _require_positive_state("mu", _state_value(state, "mu"))
        k = _require_positive_state("k", _state_value(state, "k"))
        cp = _require_positive_state("cp", _state_value(state, "cp"))
        try:
            D_h = require_geometry("geometry D_h", geometry["D_h"])
            A = require_geometry("geometry A_flow", geometry["A_flow"])
        except KeyError as exc:
            raise KeyError(f"convection geometry is missing {exc} (required: D_h, A_flow)") from exc
        Nu_lam = self._laminar(geometry)
        bad = bad_temperature(T_b, T_w)

        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            Re = m_dot * D_h / (A * mu)
            Re_safe = np.sqrt(Re * Re + self.Re_min_safe ** 2)
            Pr = mu * cp / k
            r_T = T_w / T_b
            r_eff = 0.5 * (r_T + 1.0 + np.sqrt((r_T - 1.0) ** 2 + self.eps_T_ratio_clamp ** 2))
            Nu_turb = 0.021 * Re_safe ** 0.8 * Pr ** 0.4 * r_eff ** -0.5
            w_tr = 0.5 * (1.0 + np.tanh((Re_safe - self.Re_tr_mid) / self.Re_tr_width))
            Nu_fc = (1.0 - w_tr) * Nu_lam + w_tr * Nu_turb

            if "rho" in state:
                rho = _require_positive_state("rho", as_float(state["rho"]))
                dT = np.sqrt((T_w - T_b) ** 2 + self.dT_eps_nc ** 2)
                Gr = self.gravity * (1.0 / T_b) * dT * D_h ** 3 * rho * rho / (mu * mu)
                Ra = Gr * Pr
                Nu_nc = (0.825 + 0.387 * Ra ** (1.0 / 6.0)
                         / (1.0 + (0.492 / Pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)) ** 2
                Ri = Gr / Re_safe ** 2
            elif self.include_nc:
                raise KeyError("convection state is missing 'rho' (required when include_nc=True)")
            else:
                Gr = Ra = Nu_nc = Ri = np.full(np.broadcast(Re, T_b).shape, np.nan)

            if self.include_nc:
                m = self.m_smooth_max
                Nu = (Nu_fc ** m + Nu_nc ** m) ** (1.0 / m)
            else:
                Nu = Nu_fc

        result = {"Re": Re, "Re_safe": Re_safe, "Pr": Pr, "r_T": r_T, "w_tr": w_tr,
                  "Nu_lam": np.full(np.shape(Nu), Nu_lam), "Nu_turb": Nu_turb, "Nu_fc": Nu_fc,
                  "Gr": Gr, "Ra": Ra, "Nu_nc": Nu_nc, "Ri": Ri, "Nu": Nu}
        result = {key: mask(np.broadcast_to(value, np.shape(Nu)), bad) for key, value in result.items()}
        result["h_conv"] = coefficient_from_nusselt(result["Nu"], k, D_h)
        result["mceligot_nominal"] = (np.abs(as_float(result["Re"])) >= 1.0e4) & (as_float(result["r_T"]) >= 1.0)
        return result

    def nusselt(self, state, geometry):
        return self.evaluate(state, geometry)["Nu"]

    def coefficient(self, state, geometry):
        return self.evaluate(state, geometry)["h_conv"]


def make_convection_correlation(mode_nu="mceligot", **parameters):
    """``mode_nu`` 설정으로 대류 상관식을 고른다. 보고서 [2] C의 민감도 대안(taylor, gnielinski,
    dittus_boelter)은 아직 구현하지 않았으므로 명시적으로 오류를 낸다."""
    if mode_nu not in SUPPORTED_NU_MODES:
        raise ValueError(f"unsupported or unimplemented mode_nu {mode_nu!r}; supported: {SUPPORTED_NU_MODES}")
    return McEligotConvection(**parameters)


__all__ = ["HeatTransferCorrelation", "McEligotConvection", "make_convection_correlation", "NU_LAM_CIRCULAR",
           "NU_LAM_PARALLEL_PLATES", "GRAVITY"]
