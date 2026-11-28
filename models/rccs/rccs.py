"""RCCS (Role 4).

``RCCS``        RCCS-1 공기 순환 모델 자리. 공기 경로 입력(패널·굴뚝, chimney_height_m 상충)이 없어 미구현.
``RCCS0Sink``   RCCS-0 고정 온도 열침원 (독립 상태 없음). RPV 외면 → 원자로 공동 → RCCS 경계(r 4.633 m)로
                복사와 정체 공기 전도를 한 경로로 계산하는 ThermalLink 전달 함수다 (side a = RPV, side b = 열침원).
                Q_i = ε_eff σ A_i (T_i⁴ − T_s⁴) + 2π k_air Δz_i (T_i − T_s) / ln(r_s / r_o)
                ε_eff = 1/(1/ε_RPV + (r_o/r_s)(1/ε_RCCS − 1)) (동심 원통), A_i = 2π r_o Δz_i
                R1 Table I.10: 경계 30 °C, ε 0.80, 하부 단열. RPV ε 0.85, 공기 0.068 W/(m·K) (R1 Table AIV.1).
"""
from __future__ import annotations

import math

import numpy as np

from interfaces.component import Component
from models._inputs import GEOMETRY, ModelInputMissing, require
from models.thermal_link import ThermalLink
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN, effective_emissivity_gray_enclosure


class RCCS(Component):
    def __init__(self, config=None, air_properties=None):
        self.config = {} if config is None else config
        self.air_properties = air_properties

    def register(self, registry):
        raise NotImplementedError("RCCS requires models/rccs/rccs.py Role 4 implementation and air-path inputs")

    def residual(self, t, y, ydot, context):
        raise NotImplementedError("RCCS air-loop residual requires Role 4 implementation")


class RCCS0Sink:
    """RCCS-0 고정 온도 경계. ``config``는 geometry.yaml ``rccs`` dict, ``vessel``은 ``Vessel``."""

    def __init__(self, config, vessel, properties, *, T_sink, name="rpv_rccs", include_air_conduction=True):
        if not config or "boundary" not in config:
            raise ModelInputMissing(f"{GEOMETRY}: missing rccs.boundary")
        b = config["boundary"]
        self.r_s = require(b, "radius_m", "rccs.boundary")
        self.eps_s = require(b, "emissivity", "rccs.boundary")
        if T_sink is None:
            raise ModelInputMissing("config/operation.yaml: missing boundaries.rccs.boundary_temperature_K")
        self.T_s = float(T_sink)
        self.vessel = vessel
        if self.r_s <= vessel.r_o:
            raise ValueError("RCCS boundary radius must exceed the RPV outer radius")
        self.eps_eff = float(effective_emissivity_gray_enclosure(properties.emissivity("pressure_vessel"),
                                                                 self.eps_s, vessel.r_o / self.r_s))
        self.area = 2.0 * math.pi * vessel.r_o * vessel.dz
        self.G_air = (2.0 * math.pi * properties.fixed_conductivity("air") * vessel.dz / math.log(self.r_s / vessel.r_o)
                      if include_air_conduction else np.zeros(vessel.n))
        self.name = name

    def heat_removed(self, T):
        T = np.asarray(T, dtype=float)
        return (self.eps_eff * STEFAN_BOLTZMANN * self.area * (T ** 4 - self.T_s ** 4)
                + self.G_air * (T - self.T_s))

    def transfer(self, t, y, context):
        Q = self.heat_removed(context.get(y, self.vessel.T_key))
        return -Q, np.array([Q.sum()])

    def link(self):
        return ThermalLink(self.name, self.transfer)
