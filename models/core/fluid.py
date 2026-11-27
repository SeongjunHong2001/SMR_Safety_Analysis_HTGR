"""Role 2 노심 유체 모델.

``CoreFluid``는 상위 ``Core``가 요구하는 계약(Protocol)이고, ``ParallelChannelCoreFluid``가 그 실제 구현이다.

모델 (팀 보고서 [2] 파트 B §3, 파트 C §2·§4, 파트 D §4.3)
    상부 플레넘(입구 포트)과 하부 플레넘(출구 포트) 사이의 병렬 1D 유로 M개(기본: average, hot, bypass) × N셀.
    셀 엔탈피 h_he (미분), 셀 국부 압력 p_fluid_cell (대수), 면 유량 m_dot_fluid_face (기본 대수, N+1개).
    밀도·물성은 전역 p_sys(지정 하부 플레넘이 등록)와 셀 엔탈피로 평가한다 (저마하수).
    잔차: 유로마다 질량 N + 에너지 N + 운동량 N+1, 포트마다 유량 합·유출 엔탈피 혼합 2개,
          Type I 고정 우회 분율을 쓰면 분율식 1개 (대응 미지수: 우회 보정 압력손실 dp_bypass_cal).
    포트 4개 상태(p, m_dot, h_outflow, h_instream) 중 나머지 2개 식은 접속부(StreamJunction)가 준다.

면 운동량 (flow.momentum_residual)
    내부 면: (p_{i+1} − p_i) + 열가속·단면변화 + 반셀 마찰 + 형상손실(face_K) + 반셀 중력 = 0
    입구 면 0: (p_0 − p_in) + ṁ²/(2ρA²) + K_in 손실 + 셀 0 반셀 마찰·중력 = 0   (정방향 1.5 q 강하, K_in = 0.5)
    출구 면 N: (p_out − p_{N−1}) − ṁ²/(2ρA²) + K_out 손실 + 셀 N−1 반셀 마찰·중력 = 0 (정방향 순변화 0, K_out = 1)
    역류 K는 입구 1.0, 출구 0.5 (보고서 [2] C 표, 날카로운 입구·출구).

열률
    ``residual(..., heat_rate=...)``의 heat_rate는 유체 셀로 들어오는 W 배열이고, 유로 순서대로 이어 붙인 (M·N,)이다
    (``cell_order``). 고체 온도·출력은 등록하지 않는다.
"""
from __future__ import annotations

import math
from typing import Protocol

import numpy as np

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from physics import flow
from physics import pressure_loss as pl


class CoreFluid(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Fluid-cell residual with current-trial heat received in W."""
        ...
    def ports(self): ...
    def residual_ids(self): ...
    def residual_scales(self): ...
    def outputs(self, t, y, context): ...


BYPASS_MODES = ("type_I_fixed_fraction", "calibrated_parallel")
_GEO = "config/geometry.yaml"


def _required(mapping, key, where):
    value = mapping.get(key) if hasattr(mapping, "get") else None
    if value is None:
        raise ValueError(f"{_GEO}: missing required input {where}.{key}")
    return value


def _smooth_relu(x, eps):
    """max(x, 0)의 매끄러운 근사 ½(x + √(x² + ε²)) (항상 양수)."""
    return 0.5 * (x + np.sqrt(x * x + eps * eps))


class _Channel:
    """한 병렬 유로의 셀·면 기하 (SI)."""

    def __init__(self, label, cfg, z_faces, overrides):
        cfg = {**dict(cfg), **dict(overrides or {})}
        where = f"core.channels.{label}"
        self.label = label
        n = len(z_faces) - 1
        dz = -np.diff(z_faces)
        if np.any(dz <= 0):
            raise ValueError(f"{_GEO}: core.elevation_profile_m must decrease in flow order")
        self.dz = dz                                   # 셀 길이 [m]
        self.dz_elev = np.diff(z_faces)                # 유동 방향 고도 변화 (하향 음수)
        area = cfg.get("cell_flow_area_m2")
        if area is None:
            area = np.full(n, float(_required(cfg, "flow_area_m2", where)))
        self.A = np.asarray(area, dtype=float)
        dh = cfg.get("cell_hydraulic_diameter_m")
        if dh is None:
            dh = np.full(n, float(_required(cfg, "hydraulic_diameter_m", where)))
        self.D_h = np.asarray(dh, dtype=float)
        if self.A.shape != (n,) or self.D_h.shape != (n,):
            raise ValueError(f"{_GEO}: {where} cell arrays must have {n} entries")
        if np.any(self.A <= 0) or np.any(self.D_h <= 0):
            raise ValueError(f"{_GEO}: {where} areas and hydraulic diameters must be positive")
        shape = cfg.get("shape", "circular")
        if shape not in ("circular", "parallel_plates"):
            raise ValueError(f"{where}.shape must be circular or parallel_plates")
        self.shape = shape
        self.phi_s = pl.PHI_S_CIRCULAR if shape == "circular" else pl.PHI_S_PARALLEL_PLATES
        self.roughness = float(cfg.get("roughness_m", 0.0) or 0.0)
        self.V = self.A * self.dz
        self.n = n


class ParallelChannelCoreFluid(Component):
    """병렬 유로 노심 유체 (``CoreFluid`` 구현, ``HeliumInventory`` 제공).

    config: geometry.yaml 전체 dict (``core.elevation_profile_m``, ``core.channels.*`` 사용)
    properties: ``HeliumProperties`` (예: ``physics.properties.make_helium_properties()``)
    initial: dict(m_dot [kg/s, 노심 총유량], T_inlet [K], T_outlet [K], p_inlet [Pa]) — 초기 추정값
    """

    def __init__(self, config=None, properties=None, *, initial=None, channels=None,
                 mode_bypass="type_I_fixed_fraction", bypass_channel="bypass", bypass_fraction=None,
                 bypass_K_cal=None, channel_overrides=None, inlet_K=(0.5, 1.0), outlet_K=(1.0, 0.5),
                 face_K=None, eps_mdot_rel=1.0e-4, delta_inert=False, p_sys_key="p_sys",
                 name="core.fluid"):
        if config is None or properties is None:
            raise ValueError("ParallelChannelCoreFluid requires config (geometry.yaml) and properties")
        if initial is None:
            raise ValueError("ParallelChannelCoreFluid requires initial=dict(m_dot, T_inlet, T_outlet, p_inlet)")
        self.name = name
        self.properties = properties
        self.p_sys_key = p_sys_key
        if not isinstance(delta_inert, (bool, np.bool_)):
            raise TypeError("delta_inert must be a bool")
        self.delta_inert = bool(delta_inert)
        core = config.get("core") if hasattr(config, "get") else None
        if core is None:
            raise ValueError(f"{_GEO}: missing core section")
        z = np.asarray(_required(core, "elevation_profile_m", "core"), dtype=float)
        ch_cfg = _required(core, "channels", "core")
        labels = tuple(channels) if channels is not None else tuple(ch_cfg)
        overrides = channel_overrides or {}
        self.channels = tuple(_Channel(lab, _required(ch_cfg, lab, "core.channels"), z, overrides.get(lab))
                              for lab in labels)
        self.n = self.channels[0].n
        self.inlet_K = tuple(float(v) for v in inlet_K)
        self.outlet_K = tuple(float(v) for v in outlet_K)
        self.face_K = {} if face_K is None else {k: np.asarray(v, dtype=float) for k, v in face_K.items()}
        for lab, arr in self.face_K.items():
            if arr.shape != (self.n + 1,) or np.any(arr < 0):
                raise ValueError(f"face_K[{lab}] must be nonnegative with {self.n + 1} entries")

        # 우회 유로 처리
        if mode_bypass not in BYPASS_MODES:
            raise ValueError(f"mode_bypass must be one of {BYPASS_MODES}")
        self.mode_bypass = mode_bypass
        self.bypass_index = labels.index(bypass_channel) if bypass_channel in labels else None
        if self.bypass_index is not None and mode_bypass == "type_I_fixed_fraction":
            if bypass_fraction is None:
                frac = _required(ch_cfg[bypass_channel], "type_I_fraction", f"core.channels.{bypass_channel}")
                bypass_fraction = sum(float(v) for v in frac.values())
            if not 0 < bypass_fraction < 1:
                raise ValueError("bypass_fraction must be in (0, 1)")
        elif self.bypass_index is not None and mode_bypass == "calibrated_parallel":
            if bypass_K_cal is None:
                raise ValueError("mode_bypass='calibrated_parallel' requires bypass_K_cal (S0에서 보정한 우회 손실계수)")
            arr = self.face_K.get(bypass_channel, np.zeros(self.n + 1)).copy()
            arr[-1] += float(bypass_K_cal)
            self.face_K[bypass_channel] = arr
        self.bypass_fraction = bypass_fraction if mode_bypass == "type_I_fixed_fraction" else None

        # 초기 추정값과 스케일
        for key in ("m_dot", "T_inlet", "T_outlet", "p_inlet"):
            if key not in initial or not math.isfinite(float(initial[key])):
                raise ValueError(f"initial[{key!r}] is required")
        self.initial = {k: float(v) for k, v in initial.items()}
        areas = np.array([ch.A[0] for ch in self.channels])
        share = areas / areas.sum()
        if self.bypass_fraction is not None:
            others = [i for i in range(len(self.channels)) if i != self.bypass_index]
            share[others] = share[others] / share[others].sum() * (1.0 - self.bypass_fraction)
            share[self.bypass_index] = self.bypass_fraction
        self.m_nominal = np.abs(self.initial["m_dot"]) * share
        self.m_eps = np.array([pl.mdot_epsilon(m, eps_mdot_rel) for m in self.m_nominal])
        cp = float(properties.cp(self.initial["p_inlet"], self.initial["T_inlet"]))
        self._h_scale = cp * max(abs(self.initial["T_outlet"] - self.initial["T_inlet"]), 10.0)

    # ------------------------------------------------------------------ 키
    def _ck(self, ch, var):
        return self.key(f"{ch.label}.{var}")

    def _pk(self, port, field):
        return self.key(f"{port}.{field}")

    @property
    def cell_order(self):
        """heat_rate·출력 배열의 셀 순서: ((유로, 셀 인덱스), ...)."""
        return tuple((ch.label, i) for ch in self.channels for i in range(ch.n))

    # -------------------------------------------------------------- 등록
    def register(self, registry):
        he = self.properties
        p0 = self.initial["p_inlet"]
        kind_m = flow.face_flow_state_kind(self.delta_inert)
        T_face = np.linspace(self.initial["T_inlet"], self.initial["T_outlet"], self.n + 1)
        h0 = np.asarray(he.enthalpy(p0, 0.5 * (T_face[1:] + T_face[:-1])), dtype=float)
        for ch, m0 in zip(self.channels, self.m_nominal):
            registry.register(self._ck(ch, "h_he"), kind="differential", unit="J/kg", shape=(ch.n,),
                              initial=h0, scale=float(np.max(np.abs(h0))))
            registry.register(self._ck(ch, "p_fluid_cell"), kind="algebraic", unit="Pa", shape=(ch.n,),
                              initial=p0, scale=1.0e5)
            registry.register(self._ck(ch, "m_dot_fluid_face"), kind=kind_m, unit="kg/s", shape=(ch.n + 1,),
                              initial=m0 * np.sign(self.initial["m_dot"] or 1.0), scale=max(m0, 1e-6))
        if self.bypass_fraction is not None:
            registry.register(self._ck(self.channels[self.bypass_index], "dp_bypass_cal"), kind="algebraic",
                              unit="Pa", initial=0.0, scale=1.0e3)
        h_in = float(he.enthalpy(p0, self.initial["T_inlet"]))
        h_out = float(he.enthalpy(p0, self.initial["T_outlet"]))
        m_tot = self.initial["m_dot"]
        for port, m_port, h_port in (("helium_inlet", m_tot, h_in), ("helium_outlet", -m_tot, h_out)):
            registry.register(self._pk(port, "p"), kind="algebraic", unit="Pa", initial=p0, scale=1.0e5)
            registry.register(self._pk(port, "m_dot"), kind="algebraic", unit="kg/s", initial=m_port,
                              scale=max(abs(m_tot), 1e-6))
            registry.register(self._pk(port, "h_outflow"), kind="algebraic", unit="J/kg", initial=h_port,
                              scale=abs(h_port) or 1.0)
            registry.register(self._pk(port, "h_instream"), kind="algebraic", unit="J/kg", initial=h_port,
                              scale=abs(h_port) or 1.0)

    def ports(self):
        return {port: FluidStreamPort(self._pk(port, "p"), self._pk(port, "m_dot"),
                                      self._pk(port, "h_outflow"), self._pk(port, "h_instream"),
                                      medium="helium")
                for port in ("helium_inlet", "helium_outlet")}

    # ------------------------------------------------------------- 상태 읽기
    def _channel_state(self, ch, y, ydot, context, p_sys):
        he = self.properties
        h = np.asarray(context.get(y, self._ck(ch, "h_he")), dtype=float)
        T = np.asarray(he.temperature(p_sys, h), dtype=float)
        return dict(
            h=h, h_dot=np.asarray(context.get(ydot, self._ck(ch, "h_he")), dtype=float),
            p=np.asarray(context.get(y, self._ck(ch, "p_fluid_cell")), dtype=float),
            m=np.asarray(context.get(y, self._ck(ch, "m_dot_fluid_face")), dtype=float),
            m_dot_dot=np.asarray(context.get(ydot, self._ck(ch, "m_dot_fluid_face")), dtype=float),
            T=T, rho=np.asarray(he.density(p_sys, h), dtype=float),
            drho_dp=np.asarray(he.drho_dp(p_sys, h), dtype=float),
            drho_dh=np.asarray(he.drho_dh(p_sys, h), dtype=float),
            mu=np.asarray(he.viscosity(p_sys, T), dtype=float))

    def _face_terms(self, ch, s, k_index):
        return pl.channel_face_terms(s["m"], s["rho"], s["mu"], ch.A, ch.D_h, ch.dz, ch.dz_elev,
                                     m_eps=self.m_eps[k_index], roughness=ch.roughness, phi_s=ch.phi_s,
                                     inlet=("plenum", *self.inlet_K), outlet=("plenum", *self.outlet_K),
                                     face_K=self.face_K.get(ch.label))

    def face_pressure_terms(self, t, y, context):
        """유로별 면 운동량 압력 항 [Pa] {유로: dict(acc, fric, form, grav, m, rho, A)} (잔차와 같은 계산).

        비가역 손실 = fric + form. Type I 고정 우회 분율이면 우회 유로에 ``dp_cal``(출구 면 보정 압력손실)도 준다.
        """
        p_sys = float(context.get(y, self.p_sys_key))
        out = {}
        for k, ch in enumerate(self.channels):
            s = self._channel_state(ch, y, np.zeros_like(y), context, p_sys)
            terms = {**self._face_terms(ch, s, k), "m": s["m"], "rho": s["rho"],
                     "A": np.broadcast_to(ch.A, (ch.n,)).copy()}
            if self.bypass_fraction is not None and k == self.bypass_index:
                terms["dp_cal"] = float(context.get(y, self._ck(ch, "dp_bypass_cal")))
            out[ch.label] = terms
        return out

    def _momentum(self, ch, s, k_index, p_in, p_out):
        n = ch.n
        p_left = np.concatenate(([p_in], s["p"]))
        p_right = np.concatenate((s["p"], [p_out]))
        terms = self._face_terms(ch, s, k_index)
        acc, fric, form, grav = terms["acc"], terms["fric"], terms["form"], terms["grav"]
        dp_cal = 0.0
        if self.bypass_fraction is not None and k_index == self.bypass_index:
            dp_cal = np.zeros(n + 1)
            dp_cal[-1] = s["dp_cal"]
        l_over_A = flow.face_inertia_length(np.concatenate(([0.0], ch.dz)), np.concatenate(([1.0], ch.A)),
                                            np.concatenate((ch.dz, [0.0])), np.concatenate((ch.A, [1.0])))
        r = flow.momentum_residual(p_left, p_right, dp_acc=acc, dp_fric=fric, dp_form=form + dp_cal,
                                   dp_grav=grav, include_inertia=self.delta_inert, l_over_A=l_over_A,
                                   m_dot_dot=s["m_dot_dot"])
        return np.asarray(r, dtype=float)

    # -------------------------------------------------------------- 잔차
    def residual(self, t, y, ydot, context, *, heat_rate):
        heat_rate = np.asarray(heat_rate, dtype=float)
        n_cells = sum(ch.n for ch in self.channels)
        if heat_rate.shape != (n_cells,):
            raise ValueError(f"{self.name}: heat_rate must have shape ({n_cells},) in cell_order, got {heat_rate.shape}")
        p_sys = float(context.get(y, self.p_sys_key))
        p_sys_dot = float(context.get(ydot, self.p_sys_key))
        g = lambda key: float(context.get(y, key))  # noqa: E731
        p_in, p_out = g(self._pk("helium_inlet", "p")), g(self._pk("helium_outlet", "p"))
        h_in_stream, h_out_stream = g(self._pk("helium_inlet", "h_instream")), g(self._pk("helium_outlet", "h_instream"))

        parts, offset = [], 0
        first_faces, last_faces, first_h, last_h = [], [], [], []
        for k, ch in enumerate(self.channels):
            s = self._channel_state(ch, y, ydot, context, p_sys)
            if self.bypass_fraction is not None and k == self.bypass_index:
                s["dp_cal"] = g(self._ck(ch, "dp_bypass_cal"))
            q = heat_rate[offset:offset + ch.n]
            offset += ch.n
            cell = flow.cell_residuals(V=ch.V, rho=s["rho"], drho_dp=s["drho_dp"], drho_dh=s["drho_dh"], h=s["h"],
                                       h_dot=s["h_dot"], p_sys_dot=p_sys_dot, q_wall=q, m_faces=s["m"],
                                       h_inlet_boundary=h_in_stream, h_outlet_boundary=h_out_stream,
                                       m_eps=self.m_eps[k])
            parts += [cell["R_mass"], cell["R_energy"], self._momentum(ch, s, k, p_in, p_out)]
            first_faces.append(s["m"][0])
            last_faces.append(s["m"][-1])
            first_h.append(s["h"][0])
            last_h.append(s["h"][-1])

        first_faces, last_faces = np.array(first_faces), np.array(last_faces)
        # 포트 유출 엔탈피: 포트 밖으로 나가는 유로 유량으로 가중한 혼합 (항상 정의되도록 매끄러운 가중)
        w_in = _smooth_relu(-first_faces, self.m_eps)
        w_out = _smooth_relu(last_faces, self.m_eps)
        h_mix_in = float(np.dot(w_in, first_h) / w_in.sum())
        h_mix_out = float(np.dot(w_out, last_h) / w_out.sum())
        ports = np.array([
            g(self._pk("helium_inlet", "m_dot")) - first_faces.sum(),
            g(self._pk("helium_inlet", "h_outflow")) - h_mix_in,
            g(self._pk("helium_outlet", "m_dot")) + last_faces.sum(),
            g(self._pk("helium_outlet", "h_outflow")) - h_mix_out,
        ])
        parts.append(ports)
        if self.bypass_fraction is not None:
            parts.append(np.array([first_faces[self.bypass_index] - self.bypass_fraction * first_faces.sum()]))
        return np.concatenate(parts)

    def residual_ids(self):
        ids = []
        for ch in self.channels:
            ids += [self._ck(ch, f"mass[{i}]") for i in range(ch.n)]
            ids += [self._ck(ch, f"energy[{i}]") for i in range(ch.n)]
            ids += [self._ck(ch, f"momentum[{j}]") for j in range(ch.n + 1)]
        ids += [self._pk("helium_inlet", "flow_sum"), self._pk("helium_inlet", "outflow_enthalpy"),
                self._pk("helium_outlet", "flow_sum"), self._pk("helium_outlet", "outflow_enthalpy")]
        if self.bypass_fraction is not None:
            ids.append(self._ck(self.channels[self.bypass_index], "type_I_fraction"))
        return tuple(ids)

    def residual_scales(self):
        scales = []
        for ch, m0 in zip(self.channels, self.m_nominal):
            m0 = max(m0, 1e-6)
            scales += [np.full(ch.n, m0), np.full(ch.n, m0 * self._h_scale / ch.n), np.full(ch.n + 1, 1.0e3)]
        m_tot = max(abs(self.initial["m_dot"]), 1e-6)
        scales.append(np.array([m_tot, self._h_scale, m_tot, self._h_scale]))
        if self.bypass_fraction is not None:
            scales.append(np.array([m_tot]))
        return np.concatenate(scales)

    # ------------------------------------------------------- 재고량·출력
    def helium_volumes(self):
        return np.concatenate([ch.V for ch in self.channels])

    def helium_inventory(self, t, y, context):
        p_sys = float(context.get(y, self.p_sys_key))
        rho = [np.asarray(self.properties.density(p_sys, context.get(y, self._ck(ch, "h_he"))), dtype=float)
               for ch in self.channels]
        return np.concatenate(rho) * self.helium_volumes()

    def coupling_state(self, t, y, context):
        """Role 3 결합용 현재 trial 상태 (셀 순서 = cell_order, 모두 (M·N,)).

        T_bulk, rho, mu, k, cp, m_dot(셀 양쪽 면 평균), p_sys, D_h, A_flow, dz, shape(유로별 문자열 튜플)
        """
        he = self.properties
        p_sys = float(context.get(y, self.p_sys_key))
        cols = {key: [] for key in ("T_bulk", "rho", "mu", "k", "cp", "m_dot", "D_h", "A_flow", "dz")}
        shapes = []
        for ch in self.channels:
            h = np.asarray(context.get(y, self._ck(ch, "h_he")), dtype=float)
            m = np.asarray(context.get(y, self._ck(ch, "m_dot_fluid_face")), dtype=float)
            T = np.asarray(he.temperature(p_sys, h), dtype=float)
            cols["T_bulk"].append(T)
            cols["rho"].append(np.asarray(he.density(p_sys, h), dtype=float))
            cols["mu"].append(np.asarray(he.viscosity(p_sys, T), dtype=float))
            cols["k"].append(np.asarray(he.conductivity(p_sys, T), dtype=float))
            cols["cp"].append(np.asarray(he.cp(p_sys, T), dtype=float))
            cols["m_dot"].append(0.5 * (m[:-1] + m[1:]))
            cols["D_h"].append(ch.D_h)
            cols["A_flow"].append(ch.A)
            cols["dz"].append(ch.dz)
            shapes += [ch.shape] * ch.n
        out = {key: np.concatenate(val) for key, val in cols.items()}
        out["p_sys"] = p_sys
        out["shape"] = tuple(shapes)
        return out

    def outputs(self, t, y, context):
        he = self.properties
        p_sys = float(context.get(y, self.p_sys_key))
        g = lambda key: float(context.get(y, key))  # noqa: E731
        out = {}
        flows = []
        for ch in self.channels:
            h = np.asarray(context.get(y, self._ck(ch, "h_he")), dtype=float)
            m = np.asarray(context.get(y, self._ck(ch, "m_dot_fluid_face")), dtype=float)
            out[self._ck(ch, "T_he")] = np.asarray(he.temperature(p_sys, h), dtype=float)
            out[self._ck(ch, "m_dot_inlet")] = float(m[0])
            flows.append(m[0])
        out[self.key("T_outlet_mixed")] = float(he.temperature(p_sys, g(self._pk("helium_outlet", "h_outflow"))))
        out[self.key("dp_plenum_to_plenum")] = g(self._pk("helium_inlet", "p")) - g(self._pk("helium_outlet", "p"))
        if self.bypass_index is not None:
            out[self.key("bypass_fraction")] = float(flows[self.bypass_index] / np.sum(flows))
        if self.bypass_fraction is not None:
            out[self.key("dp_bypass_cal")] = g(self._ck(self.channels[self.bypass_index], "dp_bypass_cal"))
        return out


def initial_from_operation(operation):
    """operation.yaml의 호환 키에서 노심 유체 초기 추정값을 만든다."""
    try:
        primary = operation["primary"]
        reactor = operation["reactor"]
        m = float(primary["initial_mass_flow_kg_s"])
        T_in = float(primary["initial_temperature_K"])
        p = float(primary["initial_pressure_Pa"])
        q = float(reactor["initial_thermal_power_W"])
    except (KeyError, TypeError) as exc:
        raise ValueError(f"config/operation.yaml: missing primary/reactor initial key {exc}") from exc
    cp = 5195.0  # 초기 추정용 (R1 c_p). 실제 계산은 provider 물성을 쓴다.
    return dict(m_dot=m, T_inlet=T_in, T_outlet=T_in + q / (m * cp), p_inlet=p)
