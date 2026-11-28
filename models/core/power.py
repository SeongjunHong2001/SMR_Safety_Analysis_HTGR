"""Role 3 노심 출력 provider.

``CorePower``(상태 없는 provider)와 ``StatefulCorePower``(점동특성 상태 소유) 계약, 그리고 지정 총출력 구현
``SpecifiedCorePower``를 둔다.

``SpecifiedCorePower`` (mode_core_power = "specified_total")
    P_i(t) = P_total(t) · f_channel[유로] · f_axial[셀]  →  고체 격자의 컴팩트(연료) 노드에만 W로 넣는다.
    분포는 operation.yaml ``reactor.power_distribution`` (현재 균일 분포 [가정], 팀 결정 2026-11-26).
    분율 합은 각각 1이어야 하고, 연료가 없는 셀·유로(우회 등)에 분율이 있으면 오류.
    핵분열 + 붕괴열 모드(mode_core_power = "fission_plus_decay")는 붕괴열 곡선(D7)과 점동특성 계수가
    없어 아직 지원하지 않는다 (사고 시나리오 단계).
"""
from __future__ import annotations

from typing import Protocol

import numpy as np

from models._inputs import ModelInputMissing

MODES = ("specified_total",)


class CorePower(Protocol):
    def evaluate(self, t, y, context):
        """Return a finite 1-D W array on the solid grid for this trial."""
        ...


class StatefulCorePower(CorePower, Protocol):
    """Optional point-kinetics state owner, with its own residual ledger."""

    def register(self, registry): ...
    def residual(self, t, y, ydot, context): ...
    def residual_ids(self): ...
    def residual_scales(self): ...


class SpecifiedCorePower:
    """지정 총출력을 고체 격자 컴팩트 노드에 배분하는 상태 없는 provider."""

    def __init__(self, solid, distribution, total_power, *, mode="specified_total"):
        if mode not in MODES:
            raise ModelInputMissing(f"mode_core_power {mode!r} is not implemented; 'fission_plus_decay' needs the "
                                    "decay-heat curve (D7) and point-kinetics coefficients")
        self.mode = mode
        self.solid = solid
        self.total_power = total_power
        if distribution is None:
            raise ModelInputMissing("config/operation.yaml: missing reactor.power_distribution")
        axial = np.asarray(distribution["axial_fraction"], dtype=float)
        channel = dict(distribution["channel_fraction"])
        if axial.shape != (solid.n,) or np.any(axial < 0) or not np.isclose(axial.sum(), 1.0, atol=1e-9):
            raise ValueError("power_distribution.axial_fraction must be nonnegative, sum to 1 and match the axial grid")
        if not np.isclose(sum(channel.values()), 1.0, atol=1e-6):
            raise ValueError("power_distribution.channel_fraction must sum to 1")
        heated = {col.label for col in solid.columns}
        for label, frac in channel.items():
            if frac > 0 and label not in heated:
                raise ValueError(f"power fraction {frac} given to channel {label!r} that has no fuel in the solid model")
        self.shape = np.zeros(solid.grid_size)
        for col in solid.columns:
            if np.any(axial[~col.fuel] > 0):
                raise ValueError("axial power fraction on a non-fuel cell is not supported")
            self.shape[solid.slices[(col.label, "fuel")]] = channel.get(col.label, 0.0) * axial[col.fuel]
        if not np.isclose(self.shape.sum(), 1.0, atol=1e-6):
            raise ValueError("power fractions of the modeled heated channels must sum to 1")

    def total(self, t):
        value = float(self.total_power(t) if callable(self.total_power) else self.total_power)
        if not np.isfinite(value) or value < 0:
            raise ValueError("total core power must be finite and nonnegative")
        return value

    def evaluate(self, t, y, context):
        return self.total(t) * self.shape
