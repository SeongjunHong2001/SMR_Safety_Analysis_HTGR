"""Role 3 노심 고체 모델.

``CoreSolid``는 상위 ``Core``의 계약(Protocol)이고 ``UnitCellCoreSolid``가 실제 구현이다.

구조 (팀 보고서 [3] 항목 1·2, R1 §I.5·AIV, geometry.yaml heat_structures HS-01~HS-06·HS-08)
    가열 기둥 유로마다(기본 average 65기둥, hot 1기둥) 축방향 17셀:
        연료 셀(10개)  노드 2개 — 컴팩트 T_fuel, 블록 흑연 T_block (단위격자: 냉각공 1 : 연료공 N_f/N_c)
        그 밖의 셀     블록 노드 1개 (CRE, 상·하부 반사체, 전이·유량분배·포스트 블록)
    반경 링 4개 × 17셀: 내측 반사체(HS-03), 외측 교체 반사체(HS-04), 영구 반사체(HS-05), 노심 배럴(HS-06).
    평균 유로의 블록 노드가 연료 영역(r 0.824–1.743 m) 링 역할을 하며 내·외측 반사체와 반경 전도한다.
    고온 유로는 1기둥 단독 (반경 연결 없음, HS-02 outer: null).
열 경로와 저항 (냉각공 1개당 단위 길이)
    컴팩트 평균 → 표면     R_c   = (N_c/N_f) / (8π k_compact)            (균일 발열 원통의 평균–표면)
    컴팩트 → 연료공 벽      R_gap = (N_c/N_f) / (2π r_compact h_gap),  h_gap = k_He/δ + 4σ ε_eff T̄³
    연료공 벽 → 냉각공 벽   R_web = 1 / (n_pair k_web S_pair),  S_pair = 2π / acosh((p² − r_f² − r_c²)/(2 r_f r_c))
                            (두 평행 원통 형상계수, 삼각 격자에서 냉각공 하나에 이웃 연료공 n_pair = 6개 [가정])
    블록 노드는 web 가운데에 둔다: 컴팩트–블록 G_fb = N_c Δz / (R_c + R_gap + R_web/2),
    블록–냉각공 벽 G_bw = N_c Δz / (R_web/2) (결합 모델이 대류와 직렬로 쓴다).
에너지식 (노드마다, W)
    ρ c_p(T) V dT/dt − heat_rate − Σ 전도 유입 − (배럴) 외부 링크 열률 = 0
    heat_rate는 Core가 넘기는 power − Q_from_solid (고체 격자 순서 ``grid_order``).
물성: ``physics.heat_transfer.solid_properties.R1SolidProperties`` (Ex. 2b 고정 열전도도), 간극 헬륨은 Role 2 물성.
"""
from __future__ import annotations

import math
from typing import Protocol

import numpy as np

from interfaces.component import Component
from interfaces.ports import ThermalSurfacePort
from models._inputs import GEOMETRY, ModelInputMissing
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN, effective_emissivity_gray_enclosure


class CoreSolid(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Solid-cell residual with net current-trial heat received in W."""
        ...
    def ports(self): ...
    def residual_ids(self): ...
    def residual_scales(self): ...
    def outputs(self, t, y, context): ...


# 축방향 셀 영역 → 블록 재료, 공극률 키 (geometry core.cell_region, core.block_porosity) [가정: HS-08 재료 배정]
REGION = {
    "core_restraint_element": ("alloy_800H", "top_reflector_standard"),
    "upper_reflector_block": ("graphite_H451", "top_reflector_standard"),
    "upper_reflector_half_block": ("graphite_H451", "top_reflector_standard"),
    "fuel": ("fuel_block", None),
    "lower_reflector": ("graphite_H451", "bottom_reflector_standard"),
    "lower_transition_reflector": ("graphite_H451", "transition_block"),
    "flow_distribution_block": ("graphite_2020", "flow_distribution_block"),
    "support_post_block": ("graphite_2020", "post_block"),
}
RINGS = (("inner_reflector", "graphite_H451"), ("outer_reflector", "graphite_H451"),
         ("permanent_reflector", "graphite_2020"), ("core_barrel", "alloy_800H"))


def _req(mapping, key, where):
    value = mapping.get(key) if hasattr(mapping, "get") else None
    if value is None:
        raise ModelInputMissing(f"{GEOMETRY}: missing required input {where}.{key}")
    return value


class _Column:
    """가열 기둥 유로 하나의 고체 기하."""

    def __init__(self, label, n_std, n_rsc, core, z):
        fb = core["fuel_block"]
        where = "core.fuel_block"
        self.label = label
        self.dz = -np.diff(z)
        self.n = self.dz.size
        regions = list(_req(core, "cell_region", "core"))
        if len(regions) != self.n:
            raise ValueError(f"{GEOMETRY}: core.cell_region must have {self.n} entries")
        unknown = [r for r in regions if r not in REGION]
        if unknown:
            raise ValueError(f"{GEOMETRY}: unknown core.cell_region entries {unknown}")
        self.fuel = np.array([r == "fuel" for r in regions])
        self.n_fuel = int(self.fuel.sum())
        self.material = [REGION[r][0] for r in regions]
        hex_area = float(_req(core, "hex_block_area_m2", "core"))
        cols = n_std + n_rsc
        r_fh, r_comp = float(_req(fb, "fuel_hole_radius_m", where)), float(_req(fb, "compact_radius_m", where))
        r_cl, r_cs = float(_req(fb, "coolant_hole_large_radius_m", where)), float(_req(fb, "coolant_hole_small_radius_m", where))
        nf = _req(fb, "fuel_holes", where)
        ncl, ncs = _req(fb, "coolant_holes_large", where), _req(fb, "coolant_holes_small", where)
        self.N_f = n_std * nf["standard"] + n_rsc * nf["rsc"]                     # 층당 연료공
        self.N_c = n_std * (ncl["standard"] + ncs["standard"]) + n_rsc * (ncl["rsc"] + ncs["rsc"])
        coolant_area = math.pi * (r_cl ** 2 * (n_std * ncl["standard"] + n_rsc * ncl["rsc"])
                                  + r_cs ** 2 * (n_std * ncs["standard"] + n_rsc * ncs["rsc"]))
        self.r_fh, self.r_comp, self.r_cl = r_fh, r_comp, r_cl
        self.gap = float(_req(fb, "compact_graphite_gap_m", where))
        self.pitch = float(_req(fb, "hole_pitch_m", where))
        porosity = core.get("block_porosity", {}) or {}
        area_block = np.empty(self.n)
        for i, r in enumerate(regions):
            key = REGION[r][1]
            if key is None:          # 연료 셀: 흑연 = 육각 − 냉각공 − 연료공 (컴팩트는 따로)
                area_block[i] = cols * hex_area - coolant_area - self.N_f * math.pi * r_fh ** 2
            else:
                phi = porosity.get(key)
                if phi is None:
                    phi = porosity.get("top_reflector_standard")   # [가정] CRE 공극률 미정 → 상부반사체 값
                area_block[i] = cols * hex_area * (1.0 - float(phi))
        self.V_block = area_block * self.dz
        self.V_fuel = self.N_f * math.pi * r_comp ** 2 * self.dz[self.fuel]
        self.A_axial = cols * hex_area - coolant_area                            # 축방향 전도 고체 단면


class UnitCellCoreSolid(Component):
    """단위격자 2노드 연료 기둥 + 반경 링 노심 고체 (``CoreSolid`` 구현)."""

    def __init__(self, config=None, properties=None, *, helium=None, initial=None, channels=("average", "hot"),
                 rings=True, n_pair=6, outer_link=None, gap_pressure_key="p_sys", name="core.solid"):
        if config is None or properties is None or helium is None:
            raise ModelInputMissing("UnitCellCoreSolid requires config (geometry.yaml), solid properties and "
                                    "helium properties (gap conductance)")
        if initial is None or "T" not in initial:
            raise ValueError("UnitCellCoreSolid requires initial=dict(T=K or (17,) profile)")
        self.name = name
        self.props, self.helium = properties, helium
        self.gap_pressure_key = gap_pressure_key
        self.outer_link = outer_link
        core = _req(config, "core", "")
        z = np.asarray(_req(core, "elevation_profile_m", "core"), dtype=float)
        layout = _req(core, "layout", "core")
        ch_cfg = _req(core, "channels", "core")
        n_std_total, n_rsc_total = int(layout["fuel_columns_standard"]), int(layout["fuel_columns_rsc"])
        # 고온 유로 등 'average'가 아닌 가열 유로는 표준 기둥으로 보고, average가 나머지를 가진다 [가정]
        others = sum(int(_req(ch_cfg[c], "fuel_columns", f"core.channels.{c}")) for c in channels if c != "average")
        self.columns = []
        for label in channels:
            if label == "average":
                self.columns.append(_Column(label, n_std_total - others, n_rsc_total, core, z))
            else:
                self.columns.append(_Column(label, int(_req(ch_cfg[label], "fuel_columns", f"core.channels.{label}")),
                                            0, core, z))
        if "average" not in channels and rings:
            raise ValueError("radial rings connect to the 'average' channel; include it or set rings=False")
        self.n = self.columns[0].n
        self.dz = self.columns[0].dz
        self.rings = rings
        if rings:
            rad = _req(core, "radial_m", "core")
            r = {k: float(_req(rad, k, "core.radial_m")) for k in
                 ("inner_reflector_outer", "active_core_outer", "replaceable_reflector_outer",
                  "permanent_reflector_outer", "core_barrel_inner", "core_barrel_outer")}
            self.ring_bounds = {"inner_reflector": (0.0, r["inner_reflector_outer"]),
                                "outer_reflector": (r["active_core_outer"], r["replaceable_reflector_outer"]),
                                "permanent_reflector": (r["replaceable_reflector_outer"], r["permanent_reflector_outer"]),
                                "core_barrel": (r["core_barrel_inner"], r["core_barrel_outer"])}
            self.r_fuel = (r["inner_reflector_outer"], r["active_core_outer"])
        # 단위격자 형상계수 (냉각공 하나당 이웃 연료공 n_pair개)
        c0 = self.columns[0]
        arg = (c0.pitch ** 2 - c0.r_fh ** 2 - c0.r_cl ** 2) / (2.0 * c0.r_fh * c0.r_cl)
        if arg <= 1.0:
            raise ValueError("hole pitch too small for the two-cylinder shape factor")
        self.S_pair = 2.0 * math.pi / math.acosh(arg)
        self.n_pair = float(n_pair)
        self.eps_gap = float(effective_emissivity_gray_enclosure(
            properties.emissivity("fuel_compact"), properties.emissivity("graphite_H451"), c0.r_comp / c0.r_fh))
        self._build_grid()
        T0 = np.broadcast_to(np.asarray(initial["T"], dtype=float), (self.n,))
        self.T0 = T0.copy()

    # ------------------------------------------------------------ 격자
    def _build_grid(self):
        self.slices, order, start = {}, [], 0
        for col in self.columns:
            for kind, size in (("fuel", col.n_fuel), ("block", col.n)):
                self.slices[(col.label, kind)] = slice(start, start + size)
                order += [(col.label, kind, i) for i in range(size)]
                start += size
        if self.rings:
            for ring, _ in RINGS:
                self.slices[(ring,)] = slice(start, start + self.n)
                order += [(ring, "ring", i) for i in range(self.n)]
                start += self.n
        self.grid_size = start
        self.grid_order = tuple(order)

    def _state_key(self, group):
        if len(group) == 1:
            return self.key(f"{group[0]}.T")
        return self.key(f"{group[0]}.T_{group[1]}")

    def column(self, label):
        for col in self.columns:
            if col.label == label:
                return col
        raise KeyError(label)

    # ------------------------------------------------------------ 등록
    def register(self, registry):
        for group, sl in self.slices.items():
            size = sl.stop - sl.start
            if len(group) == 2 and group[1] == "fuel":
                init = self.T0[self.column(group[0]).fuel] + 100.0
            elif len(group) == 2:
                init = self.T0 + 50.0
            else:
                init = self.T0
            registry.register(self._state_key(group), kind="differential", unit="K", shape=(size,), initial=init,
                              scale=1000.0)

    def temperatures(self, y, context):
        """고체 격자 순서의 온도 벡터 (grid_size,)."""
        T = np.empty(self.grid_size)
        for group, sl in self.slices.items():
            T[sl] = np.asarray(context.get(y, self._state_key(group)), dtype=float)
        return T

    def _ydot(self, ydot, context):
        Td = np.empty(self.grid_size)
        for group, sl in self.slices.items():
            Td[sl] = np.asarray(context.get(ydot, self._state_key(group)), dtype=float)
        return Td

    # ------------------------------------------------------------ 물성·컨덕턴스
    def _gap_pressure(self, y, context):
        try:
            return float(context.get(y, self.gap_pressure_key))
        except KeyError as exc:
            raise ModelInputMissing(f"{self.name}: gap pressure key {self.gap_pressure_key!r} is not registered "
                                    "(core gas pressure for the compact–hole and reflector–barrel helium gaps)") from exc

    def _r_web(self, col):
        return 1.0 / (self.n_pair * self.props.conductivity("graphite_H451") * self.S_pair)

    def fuel_block_conductance(self, col, T_fuel, T_block, p_gap):
        """연료 셀별 컴팩트–블록 컨덕턴스 G_fb [W/K] (n_fuel,)."""
        k_c = self.props.conductivity("fuel_compact")
        T_mean = 0.5 * (T_fuel + T_block)
        h_gap = (np.asarray(self.helium.conductivity(p_gap, T_mean), dtype=float) / col.gap
                 + 4.0 * STEFAN_BOLTZMANN * self.eps_gap * T_mean ** 3)
        ratio = col.N_c / col.N_f
        R = ratio / (8.0 * math.pi * k_c) + ratio / (2.0 * math.pi * col.r_comp * h_gap) + 0.5 * self._r_web(col)
        return col.N_c * col.dz[col.fuel] / R

    def wall_conductance(self, label):
        """블록 노드–냉각공 벽 컨덕턴스 G_bw [W/K] (17,). 결합 모델이 대류 저항과 직렬로 쓴다."""
        col = self.column(label)
        return col.N_c * col.dz / (0.5 * self._r_web(col))

    def _block_k(self, col):
        return np.array([self.props.conductivity(m) for m in col.material])

    def _edges(self, T, p_gap):
        """(i, j, G) 컨덕턴스 목록. 비선형 간극(헬륨 전도 + 복사)은 현재 온도로 평가한다."""
        I, J, G = [], [], []
        dz = self.dz
        for col in self.columns:
            fs, bs = self.slices[(col.label, "fuel")], self.slices[(col.label, "block")]
            fuel_cells = np.flatnonzero(col.fuel)
            g_fb = self.fuel_block_conductance(col, T[fs], T[bs][fuel_cells], p_gap)
            I += list(range(fs.start, fs.stop))
            J += list(bs.start + fuel_cells)
            G += list(g_fb)
            k = self._block_k(col)
            g_ax = 1.0 / (dz[:-1] / (2 * k[:-1] * col.A_axial) + dz[1:] / (2 * k[1:] * col.A_axial))
            I += list(range(bs.start, bs.stop - 1))
            J += list(range(bs.start + 1, bs.stop))
            G += list(g_ax)
        if self.rings:
            avg = self.column("average")
            k_avg = self._block_k(avg)
            k_ring = {name: self.props.conductivity(mat) for name, mat in RINGS}
            # 노드 반경: 링 중간 (내측 반사체는 반지름의 절반) [가정]
            mid = {name: 0.5 * (lo + hi) for name, (lo, hi) in self.ring_bounds.items()}
            r_avg = 0.5 * sum(self.r_fuel)
            two_pi_dz = 2.0 * math.pi * dz

            def shell(r_from, r_to, k):
                return np.log(max(r_from, r_to) / min(r_from, r_to)) / (k * two_pi_dz)

            avg_b = self.slices[("average", "block")]
            # 내측 반사체 ↔ 연료 영역, 연료 영역 ↔ 외측 반사체
            for ring, r_if in (("inner_reflector", self.r_fuel[0]), ("outer_reflector", self.r_fuel[1])):
                R = shell(mid[ring], r_if, k_ring[ring]) + shell(r_avg, r_if, k_avg)
                sl = self.slices[(ring,)]
                I += list(range(sl.start, sl.stop)); J += list(range(avg_b.start, avg_b.stop)); G += list(1.0 / R)
            # 외측 ↔ 영구 반사체
            r_if = self.ring_bounds["outer_reflector"][1]
            R = shell(mid["outer_reflector"], r_if, k_ring["outer_reflector"]) + shell(
                mid["permanent_reflector"], r_if, k_ring["permanent_reflector"])
            a, b = self.slices[("outer_reflector",)], self.slices[("permanent_reflector",)]
            I += list(range(a.start, a.stop)); J += list(range(b.start, b.stop)); G += list(1.0 / R)
            # 영구 반사체 ↔ 배럴: 헬륨 간극 전도 + 복사 (현재 온도)
            r_o, r_b = self.ring_bounds["permanent_reflector"][1], self.ring_bounds["core_barrel"][0]
            a, b = self.slices[("permanent_reflector",)], self.slices[("core_barrel",)]
            T_mean = 0.5 * (T[a] + T[b])
            eps = float(effective_emissivity_gray_enclosure(self.props.emissivity("graphite_2020"),
                                                            self.props.emissivity("alloy_800H"), r_o / r_b))
            h_gap = (np.asarray(self.helium.conductivity(p_gap, T_mean), dtype=float) / (r_b - r_o)
                     + 4.0 * STEFAN_BOLTZMANN * eps * T_mean ** 3)
            R = (shell(mid["permanent_reflector"], r_o, k_ring["permanent_reflector"])
                 + 1.0 / (h_gap * r_o * two_pi_dz) + shell(mid["core_barrel"], r_b, k_ring["core_barrel"]))
            I += list(range(a.start, a.stop)); J += list(range(b.start, b.stop)); G += list(1.0 / R)
            # 링 축방향 전도
            for ring, mat in RINGS:
                lo, hi = self.ring_bounds[ring]
                A = math.pi * (hi ** 2 - lo ** 2)
                k = k_ring[ring]
                g = 1.0 / (dz[:-1] / (2 * k * A) + dz[1:] / (2 * k * A))
                sl = self.slices[(ring,)]
                I += list(range(sl.start, sl.stop - 1)); J += list(range(sl.start + 1, sl.stop)); G += list(g)
        return np.asarray(I), np.asarray(J), np.asarray(G, dtype=float)

    def conduction_net(self, T, p_gap):
        """노드별 내부 전도 순유입 [W] (grid_size,). 합은 0 (보존)."""
        I, J, G = self._edges(T, p_gap)
        q = G * (T[I] - T[J])
        net = np.zeros(self.grid_size)
        np.add.at(net, I, -q)
        np.add.at(net, J, q)
        return net

    def heat_capacity(self, T):
        """노드별 ρ c_p(T) V [J/K] (grid_size,)."""
        C = np.empty(self.grid_size)
        for col in self.columns:
            fs, bs = self.slices[(col.label, "fuel")], self.slices[(col.label, "block")]
            C[fs] = self.props.density("fuel_compact") * self.props.cp("fuel_compact", T[fs]) * col.V_fuel
            rho = np.array([self.props.density(m) for m in col.material])
            cp = np.array([self.props.cp(m, T[bs][i]) for i, m in enumerate(col.material)])
            C[bs] = rho * cp * col.V_block
        if self.rings:
            for ring, mat in RINGS:
                lo, hi = self.ring_bounds[ring]
                sl = self.slices[(ring,)]
                C[sl] = self.props.density(mat) * self.props.cp(mat, T[sl]) * math.pi * (hi ** 2 - lo ** 2) * self.dz
        return C

    def external_rate(self, context):
        """배럴 외면 링크 열률 (grid_size,) [W, 고체로 들어오면 양수]."""
        q = np.zeros(self.grid_size)
        if self.outer_link is None:
            return q
        rates = getattr(context, "heat_rates", None)
        key = f"{self.outer_link}.to_a"
        if rates is None or key not in rates:
            raise ModelInputMissing(f"{self.name}: outer_link heat rate {key!r} is missing from the context")
        rate = np.asarray(rates[key], dtype=float)
        if rate.shape != (self.n,):
            raise ValueError(f"{self.name}: outer link rate must have shape ({self.n},)")
        q[self.slices[("core_barrel",)]] = rate
        return q

    # ------------------------------------------------------------ 잔차
    def residual(self, t, y, ydot, context, *, heat_rate):
        heat_rate = np.asarray(heat_rate, dtype=float)
        if heat_rate.shape != (self.grid_size,):
            raise ValueError(f"{self.name}: heat_rate must have shape ({self.grid_size},) in grid_order")
        T = self.temperatures(y, context)
        p_gap = self._gap_pressure(y, context)
        return (self.heat_capacity(T) * self._ydot(ydot, context) - heat_rate
                - self.conduction_net(T, p_gap) - self.external_rate(context))

    def residual_ids(self):
        return tuple(self.key(f"{g[0]}.{g[1] if len(g) > 1 else 'ring'}_energy[{i}]")
                     for g, sl in self.slices.items() for i in range(sl.stop - sl.start))

    def residual_scales(self):
        return np.full(self.grid_size, 1.0e5)

    def ports(self):
        if not self.rings:
            return {}
        key = self._state_key(("core_barrel",))
        return {f"outer_surface[{i}]": ThermalSurfacePort(key, index=i) for i in range(self.n)}

    # ------------------------------------------------------------ 출력
    def compact_center_temperatures(self, label, y, context):
        """연료 셀별 컴팩트 중심온도 추정 T_fuel + q'/(8π k_c), q' = G_fb (T_fuel − T_block)/(N_f Δz) (정상 근사)."""
        col = self.column(label)
        T = self.temperatures(y, context)
        fs, bs = self.slices[(label, "fuel")], self.slices[(label, "block")]
        Tf, Tb = T[fs], T[bs][col.fuel]
        g = self.fuel_block_conductance(col, Tf, Tb, self._gap_pressure(y, context))
        q_line = g * (Tf - Tb) / (col.N_f * col.dz[col.fuel])
        return Tf + q_line / (8.0 * math.pi * self.props.conductivity("fuel_compact"))

    def outputs(self, t, y, context):
        T = self.temperatures(y, context)
        out = {}
        peak = -np.inf
        for col in self.columns:
            out[self.key(f"{col.label}.T_fuel")] = T[self.slices[(col.label, "fuel")]].copy()
            out[self.key(f"{col.label}.T_block")] = T[self.slices[(col.label, "block")]].copy()
            centers = self.compact_center_temperatures(col.label, y, context)
            out[self.key(f"{col.label}.T_compact_center")] = centers
            peak = max(peak, float(np.max(centers)))
        out[self.key("T_fuel_max")] = peak
        if self.rings:
            for ring, _ in RINGS:
                out[self.key(f"{ring}.T")] = T[self.slices[(ring,)]].copy()
        out[self.key("C_total")] = float(np.sum(self.heat_capacity(T)))   # 전체 열용량 [J/K]
        return out
