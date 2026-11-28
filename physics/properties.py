"""헬륨 물성 provider (Role 2).

``R1HeliumProperties``는 ``interfaces.contracts.HeliumProperties``의 9개 메서드를
OECD/NEA MHTGR-350 벤치마크 물성(R1 = NEA/NSC/R(2017)4, Table AIV.13)으로 구현한다.
계수는 ``data/properties/r1_benchmark_properties.yaml``의 ``helium`` 절에서 읽으며
코드에 숫자를 따로 두지 않는다 ([제안] 2026-11-26: 헬륨 물성은 R1 세트로 전 계통 통일).

식 (SI, T [K], p [Pa], h [J/kg])
    상태방정식   rho = p / (R T)                       이상기체
    엔탈피       h = c_p (T − T_ref),  T = T_ref + h / c_p
    비열         c_p = 일정,  c_v = c_p − R
    점도         mu = a T^b
    열전도도     k = a (1 + b p) T^(c (1 − d p))      (R1의 bar 계수를 Pa로 환산해 저장)
    도함수       drho/dp|_h = 1 / (R T),  drho/dh|_p = −p / (R T² c_p)

엔탈피 기준 온도 ``T_ref``의 기본값은 변수표 ``T_he_enthalpy_ref`` 제안값 0 K다.
기준을 바꿔도 같은 (p, T)의 밀도·물성은 변하지 않는다.

입력 정의역
    p > 0, T > 0이 아니면 (또는 값이 유한하지 않으면) 해당 원소에 NaN을 돌려준다.
    Newton 선탐색은 비유한 잔차를 거부하므로, 시도 상태가 잠시 비물리 영역으로 가도
    예외로 적분이 멈추지 않는다.
    문서화된 적용 범위(``valid_T_K``, ``valid_p_Pa``)를 벗어난 양의 값은 그대로 계산한다.
    범위 확인은 ``in_valid_range``/``check_valid_range``로 따로 한다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

DEFAULT_R1_PROPERTY_FILE = (
    Path(__file__).resolve().parents[1] / "data" / "properties" / "r1_benchmark_properties.yaml"
)

# [가정] 적용 범위.
#   온도: 보고서[1] §6.1 "벤치마크 물성 적용 범위는 300–1900 K 수준".
#   압력: 시나리오 범위 0.1 MPa(S3 감압 목표) ~ 7.18 MPa(안전밸브 설정압)를 여유 있게 포함.
DEFAULT_VALID_T_K = (300.0, 1900.0)
DEFAULT_VALID_P_PA = (5.0e4, 8.0e6)

SUPPORTED_EOS_MODES = ("ideal_gas",)


def _float_array(x):
    return np.asarray(x, dtype=float)


def _out(x):
    """0-d 배열은 NumPy 스칼라로, 그 외에는 배열 그대로 돌려준다."""
    return x[()] if isinstance(x, np.ndarray) and x.ndim == 0 else x


@dataclass(frozen=True)
class R1HeliumProperties:
    """R1 Table AIV.13 헬륨 물성. 생성자 인자는 모두 SI."""

    R: float                      # 비기체상수 [J/(kg·K)]
    cp_const: float               # 정압비열 [J/(kg·K)]
    mu_a: float                   # 점도 계수 [Pa·s·K^-b]
    mu_b: float
    k_a: float                    # 열전도도 계수 [W/(m·K)·K^-c]
    k_b: float                    # [1/Pa]
    k_c: float
    k_d: float                    # [1/Pa]
    T_ref: float = 0.0            # 엔탈피 영점 온도 [K]
    valid_T_K: tuple = DEFAULT_VALID_T_K
    valid_p_Pa: tuple = DEFAULT_VALID_P_PA
    source: str = "R1 Table AIV.13"

    def __post_init__(self):
        for name in ("R", "cp_const", "mu_a", "mu_b", "k_a", "k_c"):
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0):
                raise ValueError(f"R1HeliumProperties.{name} must be positive and finite, got {value!r}")
        for name in ("k_b", "k_d", "T_ref"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"R1HeliumProperties.{name} must be finite")
        if self.T_ref < 0:
            raise ValueError("T_ref is an absolute temperature and must be >= 0 K")
        if self.cp_const <= self.R:
            raise ValueError("cp must exceed R so that cv = cp - R is positive")
        for name in ("valid_T_K", "valid_p_Pa"):
            lo, hi = getattr(self, name)
            if not (0 < lo < hi and math.isfinite(hi)):
                raise ValueError(f"{name} must be (low, high) with 0 < low < high")

    # ------------------------------------------------------------------ 생성
    @classmethod
    def from_yaml(cls, path=None, *, T_ref=0.0, valid_T_K=DEFAULT_VALID_T_K,
                  valid_p_Pa=DEFAULT_VALID_P_PA):
        """``r1_benchmark_properties.yaml``의 ``helium`` 절에서 계수를 읽는다."""
        path = Path(path) if path is not None else DEFAULT_R1_PROPERTY_FILE
        with path.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
        try:
            he = data["helium"]
            if he.get("density_model") != "ideal_gas":
                raise ValueError(f"{path}: helium.density_model must be 'ideal_gas'")
            mu, k = he["viscosity"], he["conductivity"]
            values = dict(R=he["specific_gas_constant_J_kgK"], cp_const=he["cp_J_kgK"],
                          mu_a=mu["a"], mu_b=mu["b"], k_a=k["a"], k_b=k["b"], k_c=k["c"], k_d=k["d"])
        except KeyError as exc:
            raise KeyError(f"{path}: missing helium property key {exc}") from exc
        for key, value in values.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{path}: helium coefficient {key} must be a number, got {value!r}")
        return cls(**{k: float(v) for k, v in values.items()}, T_ref=float(T_ref),
                   valid_T_K=tuple(valid_T_K), valid_p_Pa=tuple(valid_p_Pa),
                   source=str(he.get("source", "R1 Table AIV.13")))

    # ------------------------------------------------------------ 내부 도구
    def _p_T(self, p, T):
        p, T = np.broadcast_arrays(_float_array(p), _float_array(T))
        bad = ~(np.isfinite(p) & np.isfinite(T) & (p > 0) & (T > 0))
        return p, T, bad

    def _T_from_h(self, h):
        return self.T_ref + _float_array(h) / self.cp_const

    @staticmethod
    def _mask(value, bad):
        value = np.array(value, dtype=float, copy=True)
        value[bad] = np.nan
        return _out(value)

    # ------------------------------------------------- HeliumProperties 계약
    def temperature(self, p, h):
        p, T, bad = self._p_T(p, self._T_from_h(h))
        return self._mask(T, bad)

    def enthalpy(self, p, T):
        p, T, bad = self._p_T(p, T)
        return self._mask(self.cp_const * (T - self.T_ref), bad)

    def density(self, p, h):
        p, T, bad = self._p_T(p, self._T_from_h(h))
        with np.errstate(divide="ignore", invalid="ignore"):
            return self._mask(p / (self.R * T), bad)

    def viscosity(self, p, T):
        p, T, bad = self._p_T(p, T)
        with np.errstate(invalid="ignore"):
            return self._mask(self.mu_a * np.power(T, self.mu_b), bad)

    def conductivity(self, p, T):
        p, T, bad = self._p_T(p, T)
        with np.errstate(invalid="ignore"):
            k = self.k_a * (1.0 + self.k_b * p) * np.power(T, self.k_c * (1.0 - self.k_d * p))
        return self._mask(k, bad)

    def cp(self, p, T):
        p, T, bad = self._p_T(p, T)
        return self._mask(np.full(p.shape, self.cp_const), bad)

    def cv(self, p, T):
        p, T, bad = self._p_T(p, T)
        return self._mask(np.full(p.shape, self.cp_const - self.R), bad)

    def drho_dp(self, p, h):
        p, T, bad = self._p_T(p, self._T_from_h(h))
        with np.errstate(divide="ignore", invalid="ignore"):
            return self._mask(1.0 / (self.R * T), bad)

    def drho_dh(self, p, h):
        p, T, bad = self._p_T(p, self._T_from_h(h))
        with np.errstate(divide="ignore", invalid="ignore"):
            return self._mask(-p / (self.R * T * T * self.cp_const), bad)

    # ------------------------------------------------------------ 적용 범위
    def in_valid_range(self, p, T):
        """문서화된 적용 범위 안이면 True (broadcast shape의 bool)."""
        p, T, bad = self._p_T(p, T)
        (tlo, thi), (plo, phi) = self.valid_T_K, self.valid_p_Pa
        ok = ~bad & (T >= tlo) & (T <= thi) & (p >= plo) & (p <= phi)
        return _out(ok)

    def check_valid_range(self, p, T, *, where="helium state"):
        """적용 범위를 벗어난 원소가 있으면 ValueError. 초기화·진단용."""
        ok = np.asarray(self.in_valid_range(p, T))
        if not np.all(ok):
            p_b, T_b = np.broadcast_arrays(_float_array(p), _float_array(T))
            idx = np.flatnonzero(~ok.ravel())
            first = idx[0]
            raise ValueError(
                f"{where}: {idx.size} value(s) outside R1 helium range "
                f"T {self.valid_T_K} K, p {self.valid_p_Pa} Pa; "
                f"first at flat index {first}: T={T_b.ravel()[first]!r} K, p={p_b.ravel()[first]!r} Pa"
            )


def make_helium_properties(mode="ideal_gas", *, path=None, T_ref=0.0):
    """``mode_he_eos`` 설정으로 헬륨 provider를 고른다. 지원하지 않는 모드는 오류."""
    if mode not in SUPPORTED_EOS_MODES:
        raise ValueError(f"unsupported mode_he_eos {mode!r}; supported: {SUPPORTED_EOS_MODES}")
    return R1HeliumProperties.from_yaml(path, T_ref=T_ref)
