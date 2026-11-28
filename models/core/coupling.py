"""Role 3 노심 고체–유체 열 결합.

``CoreCoupling`` 계약과 ``UnitCellCoreCoupling`` 구현. 같은 trial에서 경계마다 열률을 한 번만 계산해
(Q_to_fluid, Q_from_solid)를 돌려준다 (합이 같음, Core가 검사).

가열 기둥 유로 (average, hot)
    블록 노드 → 냉각공 벽 G_bw (solid.wall_conductance)와 벽 → 헬륨 대류 G_h = h A_w를 직렬로 쓴다.
    A_w = 젖은 둘레 4A/D_h × Δz (유체 유로 기하). h는 ``McEligotConvection.coefficient``이고 벽온도가 필요하므로
    두 번 평가한다: ① T_wall = T_block으로 h₀ → T_w = (G_bw T_b + G_h₀ T_he)/(G_bw + G_h₀) ② h(T_w).
    Q = (T_block − T_he) G_bw G_h / (G_bw + G_h).  (이전 평가값 캐시 없음, 같은 trial에서 결정적)
우회 유로 (bypass)
    R1 Type I 우회 영역 분율(geometry core.channels.bypass.type_I_fraction)로 우회 젖은 면적을 고체에 나눈다 [가정]:
      intra_core → 평균 연료 블록, inner_reflector + inner_control_rod → 내측 반사체,
      outer_control_rod + outer_reflector_row1·2 → 외측 반사체, permanent_reflector → 영구 반사체.
    Q = h A_w,share (T_solid − T_bypass), 벽온도 = 고체 노드 온도.
열률 부호: 고체 → 유체 양수 (Q_to_fluid 양수, Q_from_solid 양수).
"""
from __future__ import annotations

from typing import Protocol

import numpy as np


class CoreCoupling(Protocol):
    def evaluate(self, t, y, context):
        """Return (Q_to_fluid, Q_from_solid) in W with equal sums."""
        ...


BYPASS_REGION_TO_SOLID = {
    "intra_core": ("average", "block"),
    "inner_reflector": ("inner_reflector",),
    "inner_control_rod": ("inner_reflector",),
    "outer_control_rod": ("outer_reflector",),
    "outer_reflector_row1": ("outer_reflector",),
    "outer_reflector_row2": ("outer_reflector",),
    "permanent_reflector": ("permanent_reflector",),
}


class UnitCellCoreCoupling:
    """``ParallelChannelCoreFluid``와 ``UnitCellCoreSolid``의 대류 결합."""

    def __init__(self, fluid, solid, correlation, *, bypass_regions=None, bypass_channel="bypass"):
        self.fluid, self.solid, self.corr = fluid, solid, correlation
        labels = [label for label, _ in fluid.cell_order]
        self.fluid_slices = {}
        for label in dict.fromkeys(labels):
            idx = [i for i, lab in enumerate(labels) if lab == label]
            self.fluid_slices[label] = slice(idx[0], idx[-1] + 1)
        for col in solid.columns:
            if col.label not in self.fluid_slices:
                raise ValueError(f"solid column {col.label!r} has no matching fluid channel")
        self.bypass_channel = bypass_channel if bypass_channel in self.fluid_slices else None
        self.bypass_shares = []
        if self.bypass_channel is not None and bypass_regions:
            total = sum(float(v) for v in bypass_regions.values())
            if total <= 0:
                raise ValueError("bypass region fractions must sum to a positive value")
            merged = {}
            for region, frac in bypass_regions.items():
                if region not in BYPASS_REGION_TO_SOLID:
                    raise ValueError(f"unknown bypass region {region!r}")
                target = BYPASS_REGION_TO_SOLID[region]
                if target not in solid.slices:
                    continue                      # 링 없는 모델이면 그 몫은 단열로 남김
                merged[target] = merged.get(target, 0.0) + float(frac) / total
            self.bypass_shares = list(merged.items())

    def _h(self, fs, sl, T_wall):
        state = dict(m_dot=fs["m_dot"][sl], T_bulk=fs["T_bulk"][sl], T_wall=T_wall, mu=fs["mu"][sl],
                     k=fs["k"][sl], cp=fs["cp"][sl], rho=fs["rho"][sl])
        shapes = set(fs["shape"][sl])
        if len(shapes) != 1:
            raise ValueError("one fluid channel must have a single flow shape")
        geometry = dict(D_h=fs["D_h"][sl], A_flow=fs["A_flow"][sl], shape=shapes.pop())
        return np.asarray(self.corr.coefficient(state, geometry), dtype=float)

    def evaluate(self, t, y, context):
        fs = self.fluid.coupling_state(t, y, context)
        T_s = self.solid.temperatures(y, context)
        q_fluid = np.zeros(len(self.fluid.cell_order))
        q_solid = np.zeros(self.solid.grid_size)
        for col in self.solid.columns:
            sl = self.fluid_slices[col.label]
            bs = self.solid.slices[(col.label, "block")]
            T_b, T_he = T_s[bs], fs["T_bulk"][sl]
            A_w = 4.0 * fs["A_flow"][sl] / fs["D_h"][sl] * fs["dz"][sl]
            G_bw = self.solid.wall_conductance(col.label)
            G_h0 = self._h(fs, sl, T_b) * A_w
            T_w = (G_bw * T_b + G_h0 * T_he) / (G_bw + G_h0)
            G_h = self._h(fs, sl, T_w) * A_w
            Q = (T_b - T_he) * G_bw * G_h / (G_bw + G_h)
            q_fluid[sl] += Q
            q_solid[bs] += Q
        if self.bypass_channel is not None:
            sl = self.fluid_slices[self.bypass_channel]
            T_bp = fs["T_bulk"][sl]
            P_dz = 4.0 * fs["A_flow"][sl] / fs["D_h"][sl] * fs["dz"][sl]
            for target, share in self.bypass_shares:
                ts = self.solid.slices[target]
                T_src = T_s[ts]
                Q = self._h(fs, sl, T_src) * share * P_dz * (T_src - T_bp)
                q_fluid[sl] += Q
                q_solid[ts] += Q
        return q_fluid, q_solid
