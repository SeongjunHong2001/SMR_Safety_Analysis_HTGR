"""노심 고체 재료 물성 (Role 3).

``R1SolidProperties``는 OECD/NEA MHTGR-350 벤치마크(R1) 물성을 ``data/properties/r1_benchmark_properties.yaml``에서 읽는다.

열전도도 (mode)
    "ex2b_fixed"  R1 Table AIV.1 고정값 (Ex. 2a·2b): 컴팩트 20, 균질 연료블록 37, H-451 66, 2020 35, Alloy 800H 17.8 W/(m·K)
비열 (R1 Table AIV.3)
    흑연(H-451, 2020, 연료블록·컴팩트): cp = 4184 (c0 + c1 T + c2/T + c3/T² + c4/T³ + c5/T⁴), T는 300–2000 K로 제한
    Alloy 800H 515, SA533B(pressure_vessel) 580 J/(kg·K)
밀도
    H-451 1850 (PCEA AG 대체값), 2020 1780, Alloy 800H 7800, SA533B 7500 kg/m³. RPV 열전도도 40 W/(m·K) (AIV.1).
    [가정] 연료 컴팩트·균질 연료블록의 밀도·비열은 R1 파일에 옮기지 않아 H-451 값을 쓴다 (정상해와 무관, 과도 열용량만 영향).
방사율: 고체 공통 0.85 (R1 Table AIV.1).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

DEFAULT_FILE = Path(__file__).resolve().parents[2] / "data" / "properties" / "r1_benchmark_properties.yaml"
MATERIALS = ("fuel_compact", "fuel_block", "graphite_H451", "graphite_2020", "alloy_800H", "pressure_vessel")
GRAPHITE_LIKE = ("fuel_compact", "fuel_block", "graphite_H451", "graphite_2020")


class R1SolidProperties:
    """R1 고체 물성. ``conductivity(material, T)``, ``cp(material, T)``, ``density(material)``, ``emissivity(material)``."""

    def __init__(self, path=None, *, mode="ex2b_fixed"):
        if mode != "ex2b_fixed":
            raise ValueError("only mode='ex2b_fixed' (R1 Table AIV.1) is implemented; temperature/fluence-dependent "
                             "conductivity (Ex. 2c) needs the participant-only fluence map")
        self.mode = mode
        path = Path(path) if path is not None else DEFAULT_FILE
        with path.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
        fixed = data["fixed_properties_ex2ab"]
        self._k = {m: float(fixed["conductivity_W_mK"][m]) for m in MATERIALS}
        self._fixed_k = {key: float(v) for key, v in fixed["conductivity_W_mK"].items()}   # air, coolant 포함
        self._emissivity = float(fixed["emissivity"])
        h451 = data["graphite_H451"]
        self._cp_mult = float(h451["cp"]["multiplier"])
        self._cp_c = np.asarray(h451["cp"]["c"], dtype=float)
        self._rho = {"graphite_H451": float(h451["density_kg_m3"]),
                     "graphite_2020": float(data["graphite_2020"]["density_kg_m3"]),
                     "alloy_800H": float(data["alloy_800H"]["density_kg_m3"])}
        self._rho["fuel_compact"] = self._rho["graphite_H451"]     # [가정]
        self._rho["fuel_block"] = self._rho["graphite_H451"]       # [가정]
        self._cp_800H = float(data["alloy_800H"]["cp_J_kgK"])
        rpv = data["pressure_vessel_SA533B"]                       # RPV (R1 Table AIV.12)
        self._rho["pressure_vessel"] = float(rpv["density_kg_m3"])
        self._cp_rpv = float(rpv["cp_J_kgK"])

    @staticmethod
    def _check(material):
        if material not in MATERIALS:
            raise ValueError(f"unknown solid material {material!r}; known: {MATERIALS}")

    def conductivity(self, material, T=None):
        self._check(material)
        value = self._k[material]
        return value if T is None else np.full(np.shape(T), value)[()]

    def cp(self, material, T):
        self._check(material)
        T = np.asarray(T, dtype=float)
        if material == "alloy_800H":
            return np.full(T.shape, self._cp_800H)[()]
        if material == "pressure_vessel":
            return np.full(T.shape, self._cp_rpv)[()]
        Tc = np.clip(T, 300.0, 2000.0)
        c = self._cp_c
        return (self._cp_mult * (c[0] + c[1] * Tc + c[2] / Tc + c[3] / Tc ** 2 + c[4] / Tc ** 3
                                 + c[5] / Tc ** 4))[()]

    def fixed_conductivity(self, key):
        """R1 Table AIV.1 고정 열전도도 표의 임의 항목 (예: 'air' 0.068, 'coolant' 0.41) [W/(m·K)]."""
        if key not in self._fixed_k:
            raise ValueError(f"no fixed conductivity {key!r} in R1 Table AIV.1; known: {sorted(self._fixed_k)}")
        return self._fixed_k[key]

    def density(self, material):
        self._check(material)
        return self._rho[material]

    def emissivity(self, material):
        self._check(material)
        return self._emissivity
