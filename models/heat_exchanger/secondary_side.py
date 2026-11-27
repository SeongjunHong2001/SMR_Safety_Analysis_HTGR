"""Role 4 SG 2차측.

``SecondarySide`` 계약과 SG-0용 ``FixedTemperatureWaterSink``.

SG-0 ([제안] 2026-11-26, operation.yaml ``boundaries.sg_heat_sink``)
    2차측 = 일정량·일정 온도의 액체 물 (상변화·2차측 유동 없음). He 셀 i의 제거열 Q_i = U A_i (T_He,i − T_물).
    물 포트가 없고(ports 빈 dict), 상태는 U 보정 모드에서만 하나다:
        U 보정 (정상 초기화)  ``{name}.U_overall`` (algebraic, W/(m²·K))와 잔차 T_He,출구 − T_목표 = 0
        U 고정 (과도)          상태·잔차 없음, U = 보정값
    받은 열률(heat_rate = 제거열)은 출력 ``{name}.Q_removed``로 기록한다 (상위 HeatExchanger가 넘긴 값을 그대로 합산).
"""
from __future__ import annotations

from typing import Protocol

import numpy as np

from models._inputs import ModelInputMissing


class SecondarySide(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D secondary-side residuals with current-trial W heat input."""
        ...
    def ports(self): ...
    def residual_ids(self): ...
    def residual_scales(self): ...
    def outputs(self, t, y, context): ...


class FixedTemperatureWaterSink:
    """SG-0 열침원 (``SecondarySide`` 구현)."""

    def __init__(self, sink_config, *, helium_side, U=None, calibrate_outlet_K=None, U_initial=None,
                 name="sg.secondary"):
        if not sink_config:
            raise ModelInputMissing("config/operation.yaml: missing boundaries.sg_heat_sink")
        if sink_config.get("mode") != "fixed_temperature_liquid_water" or sink_config.get("phase_change") is not False:
            raise ValueError("SG-0 sink requires mode fixed_temperature_liquid_water and phase_change: false")
        T_w = sink_config.get("water_temperature_K")
        if T_w is None:
            raise ModelInputMissing("config/operation.yaml: missing boundaries.sg_heat_sink.water_temperature_K")
        if (U is None) == (calibrate_outlet_K is None):
            raise ValueError("give exactly one of U (fixed, transient) or calibrate_outlet_K (steady calibration)")
        self.name = name
        self.T_water = float(T_w)
        self.helium_side = helium_side
        self.U_fixed = None if U is None else float(U)
        self.target = None if calibrate_outlet_K is None else float(calibrate_outlet_K)
        est = sink_config.get("overall_U_estimate_W_m2K")
        self.U0 = float(U_initial if U_initial is not None else (est if est is not None else 400.0))
        self.water_mass = sink_config.get("water_mass_kg")

    @property
    def calibrating(self):
        return self.target is not None

    @property
    def U_key(self):
        return f"{self.name}.U_overall"

    def U(self, y, context):
        return float(context.get(y, self.U_key)) if self.calibrating else self.U_fixed

    def register(self, registry):
        if self.calibrating:
            registry.register(self.U_key, kind="algebraic", unit="W/(m2*K)", initial=self.U0, scale=100.0)

    def residual(self, t, y, ydot, context, *, heat_rate):
        if not self.calibrating:
            return np.empty(0)
        fs = self.helium_side.coupling_state(t, y, context)
        return np.array([fs["T_bulk"][-1] - self.target])

    def ports(self):
        return {}

    def residual_ids(self):
        return (f"{self.name}.U_calibration",) if self.calibrating else ()

    def residual_scales(self):
        return np.array([1.0]) if self.calibrating else np.empty(0)

    def outputs(self, t, y, context):
        U = self.U(y, context)
        fs = self.helium_side.coupling_state(t, y, context)
        Q = U * self.helium_side.area_heat * (fs["T_bulk"] - self.T_water)
        return {f"{self.name}.U_overall": U, f"{self.name}.Q_removed": float(Q.sum()),
                f"{self.name}.T_water": self.T_water}
