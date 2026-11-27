"""Role 4 SG 헬륨 측.

``HeliumSide`` 계약과 ``SGHeliumSide`` 구현. SG He 측(shell, 하향 유동)을 1D 유한체적 유로로 둔다
(팀 보고서 [2] 파트 B 커널, Role 2 ``Piping``을 그대로 쓴다).

geometry.yaml ``heat_exchanger.helium_side``:
    elevation_profile_m (10,)  셀 경계 고도 (FSH 1셀 + EES 8셀, FSH–EES 간격은 FSH 셀에 포함, 팀 결정)
    shell_flow_area_total_m2   shroud 환형 전면 단면 (관 단면 미차감) [계산] 9.23 m²
    mean_flow_area_m2          다발 공극률 × 전면 단면 = 자유체적/길이 [계산] 6.66 m² — 1D 셀 체적·평균 속도에 쓴다
    min_free_flow_area_m2      관군 최소 자유유로면적 [계산] 3.85 m² (Zukauskas 교차류 상관식용, SG-2)
    hydraulic_diameter_m       4 × 자유체적 / 관 외면적 [계산] 0.0577 m. ``hydraulic_diameter_m`` 인자가 있으면 그것을 쓴다
    면적 선택 순서: flow_area_m2 인자 → mean_flow_area_m2 → shell_flow_area_total_m2.
    (최소 면적을 1D 면적으로 쓰면 He 재고량이 자유체적보다 작아지므로 쓰지 않는다.)
포트 helium_inlet (SG 입구, 노드 300 쪽), helium_outlet (330 쪽). 열률은 상위 HeatExchanger가 heat_rate로 넘긴다.
"""
from __future__ import annotations

from typing import Protocol

import numpy as np

from models._inputs import GEOMETRY, ModelInputMissing
from models.piping import Piping


class HeliumSide(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D primary-side residuals with current-trial W heat input."""
        ...
    def ports(self): ...
    def residual_ids(self): ...
    def residual_scales(self): ...
    def outputs(self, t, y, context): ...


class SGHeliumSide:
    """SG He 측 1D 유로 (``HeliumSide`` 구현, ``HeliumInventory`` 제공)."""

    def __init__(self, config, properties, *, initial, name="sg.helium", hydraulic_diameter_m=None,
                 flow_area_m2=None, inlet=("pipe", 0.0, 0.0), outlet=("pipe", 0.0, 0.0), face_K=None):
        hx = (config or {}).get("heat_exchanger") if hasattr(config, "get") else None
        if hx is None or "helium_side" not in hx:
            raise ModelInputMissing(f"{GEOMETRY}: missing heat_exchanger.helium_side")
        hs = hx["helium_side"]
        where = "heat_exchanger.helium_side"
        z = hs.get("elevation_profile_m")
        if z is None:
            raise ModelInputMissing(f"{GEOMETRY}: missing required input {where}.elevation_profile_m")
        D = hydraulic_diameter_m if hydraulic_diameter_m is not None else hs.get("hydraulic_diameter_m")
        if D is None:
            raise ModelInputMissing(f"{GEOMETRY}: missing required input {where}.hydraulic_diameter_m "
                                    "(M11 SG He 측 수력직경) — geometry 값이나 hydraulic_diameter_m 인자로 명시")
        A = flow_area_m2 if flow_area_m2 is not None else (hs.get("mean_flow_area_m2")
                                                           or hs.get("shell_flow_area_total_m2"))
        if A is None:
            raise ModelInputMissing(f"{GEOMETRY}: missing required input {where}.shell_flow_area_total_m2")
        n = len(z) - 1
        if n != int(hx["axial_cells"]):
            raise ValueError(f"{GEOMETRY}: {where}.elevation_profile_m must have axial_cells + 1 entries")
        self.area_heat = np.asarray(hx["heat_transfer_area_m2"], dtype=float)
        if self.area_heat.shape != (n,):
            raise ValueError(f"{GEOMETRY}: heat_exchanger.heat_transfer_area_m2 must have {n} entries")
        self.pipe = Piping({"z_faces_m": list(z), "flow_area_m2": float(A), "hydraulic_diameter_m": float(D)},
                           properties, name=name, initial=initial, inlet=inlet, outlet=outlet, face_K=face_K)
        self.name = name
        self.n = n

    # 상위 HeatExchanger 계약
    def register(self, registry):
        self.pipe.register(registry)

    def residual(self, t, y, ydot, context, *, heat_rate):
        return self.pipe.residual_with_heat(t, y, ydot, context, heat_rate)

    def ports(self):
        return self.pipe.ports()

    def residual_ids(self):
        return self.pipe.residual_ids()

    def residual_scales(self):
        return self.pipe.residual_scales()

    def outputs(self, t, y, context):
        out = self.pipe.outputs(t, y, context)
        out.pop(self.pipe.key("Q_wall"), None)            # 벽 열은 결합 모델이 기록
        return out

    # 재고량·결합용
    def helium_volumes(self):
        return self.pipe.helium_volumes()

    def helium_inventory(self, t, y, context):
        return self.pipe.helium_inventory(t, y, context)

    def coupling_state(self, t, y, context):
        return self.pipe.coupling_state(t, y, context)

    def face_pressure_terms(self, t, y, context):
        return self.pipe.face_pressure_terms(t, y, context)

    @property
    def h_key(self):
        return self.pipe.key("h_he")
