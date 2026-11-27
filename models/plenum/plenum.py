"""Plenum 저장량과 열에너지 (Role 2, mode_plenum = finite_volume).

0D 유한체적 플레넘. 셀 하나(h_he 미분, p_fluid_cell 대수)와 실제 접속 수만큼의 포트 ``helium_ports[i]``를 가진다.

잔차 (팀 보고서 [2] 파트 B §3.2, docs/interface_spec.md §5)
    질량   V(ρ_p ṗ_sys + ρ_h ḣ) − Σ_i ṁ_i = 0                        (ṁ_i: 포트 i 유입 양수)
    에너지 ρVḣ − Vṗ_sys + Σ_i s⁺(ṁ_i)(h − h_instream,i) − Q = 0      (s⁺: 매끄러운 유입량 max(ṁ, 0))
    포트 i (2개): p_port,i − [p + ρ g (z_ref − z_i)] = 0,  h_outflow,i − h = 0
    포트의 나머지 2개 식(압력 같음, 유입 stream)은 StreamJunction이 준다.
지정 하부 플레넘 (owns_p_sys=True)
    전역 ``p_sys``를 differential·owner ``primary_loop``로 한 번만 등록하고, 국부 기준식 R_p_ref = p − p_sys = 0을 더한다.
    정상 초기화용 인스턴스(steady_pressure_closure 지정)는 자기 질량식을 p_sys − p_sys,0 = 0으로 바꾼다
    (보고서 D2 결정: 기준압 폐합. 과도 인스턴스는 원래 질량식을 쓴다).
"""
from __future__ import annotations

import numpy as np

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from models._inputs import GEOMETRY, ModelInputMissing, require
from physics import flow
from physics import pressure_loss as pl


class Plenum(Component):
    """0D 유한체적 플레넘. ``config``는 geometry.yaml ``flow_network.nodes.<이름>`` dict."""

    def __init__(self, config=None, properties=None, *, name="plenum", initial=None, port_elevations=None,
                 port_labels=None, owns_p_sys=False, steady_pressure_closure=None, heat_rate=0.0,
                 eps_mdot_rel=1.0e-4, p_sys_key="p_sys"):
        self.config = dict(config or {})
        self.properties = properties
        self.name = name
        self.initial = None if initial is None else {k: float(v) for k, v in initial.items()}
        self.owns_p_sys = bool(owns_p_sys)
        self.steady_pressure_closure = (None if steady_pressure_closure is None
                                        else float(steady_pressure_closure))
        if self.steady_pressure_closure is not None and not self.owns_p_sys:
            raise ValueError("steady_pressure_closure is only for the p_sys-owning (lower) plenum")
        self.heat_rate = heat_rate
        self.eps_mdot_rel = float(eps_mdot_rel)
        self.p_sys_key = p_sys_key
        self._ready = False
        if self.config:
            self._setup(port_elevations, port_labels)

    def _setup(self, port_elevations, port_labels):
        where = f"flow_network.nodes.{self.name}"
        c = self.config
        self.V = require(c, "volume_m3", where)
        if port_elevations is None:
            port_elevations = [require(c, "z_inlet_m", where, positive=False),
                               require(c, "z_outlet_m", where, positive=False)]
        self.z_ports = np.asarray(port_elevations, dtype=float)
        if self.z_ports.ndim != 1 or self.z_ports.size < 2:
            raise ValueError(f"Plenum '{self.name}' needs at least two port elevations")
        self.n_ports = self.z_ports.size
        if c.get("z_bottom_m") is not None and c.get("z_top_m") is not None:
            self.z_ref = 0.5 * (float(c["z_bottom_m"]) + float(c["z_top_m"]))
        else:
            self.z_ref = float(np.mean(self.z_ports))
        self.port_labels = tuple(port_labels) if port_labels is not None else tuple(
            f"port{i}" for i in range(self.n_ports))
        if len(self.port_labels) != self.n_ports:
            raise ValueError("port_labels must match the number of ports")
        if self.initial is None or any(k not in self.initial for k in ("m_dot", "T", "p")):
            raise ValueError(f"Plenum '{self.name}' requires initial=dict(m_dot, T, p)")
        self.m_nominal = max(abs(self.initial["m_dot"]), 1e-6)
        self.m_eps = pl.mdot_epsilon(self.m_nominal, self.eps_mdot_rel)
        self._ready = True

    def _require_ready(self):
        if not self._ready:
            raise ModelInputMissing(f"Plenum '{self.name}': {GEOMETRY} flow_network.nodes.<name> config "
                                    "(volume_m3, port elevations) and properties are required")
        if self.properties is None:
            raise ModelInputMissing(f"Plenum '{self.name}': helium properties provider is required")

    def port_name(self, i):
        return f"helium_ports[{i}]"

    def port(self, label):
        """포트 라벨(예: 'to_core')로 FluidStreamPort를 찾는다."""
        return self.ports()[self.port_name(self.port_labels.index(label))]

    def _pk(self, i, field):
        return self.key(f"{self.port_name(i)}.{field}")

    # -------------------------------------------------------------- 등록
    def register(self, registry):
        self._require_ready()
        he, ini = self.properties, self.initial
        h0 = float(he.enthalpy(ini["p"], ini["T"]))
        if self.owns_p_sys:
            registry.register(self.p_sys_key, kind="differential", unit="Pa", initial=ini["p"], scale=1.0e5,
                              owner="primary_loop")
        registry.register(self.key("h_he"), kind="differential", unit="J/kg", initial=h0, scale=abs(h0))
        registry.register(self.key("p_fluid_cell"), kind="algebraic", unit="Pa", initial=ini["p"], scale=1.0e5)
        for i in range(self.n_ports):
            registry.register(self._pk(i, "p"), kind="algebraic", unit="Pa", initial=ini["p"], scale=1.0e5)
            registry.register(self._pk(i, "m_dot"), kind="algebraic", unit="kg/s", initial=0.0, scale=self.m_nominal)
            registry.register(self._pk(i, "h_outflow"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))
            registry.register(self._pk(i, "h_instream"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))

    def ports(self):
        return {self.port_name(i): FluidStreamPort(self._pk(i, "p"), self._pk(i, "m_dot"), self._pk(i, "h_outflow"),
                                                   self._pk(i, "h_instream"), medium="helium")
                for i in range(self.n_ports)}

    # -------------------------------------------------------------- 잔차
    def residual(self, t, y, ydot, context):
        self._require_ready()
        he = self.properties
        g = lambda key: float(context.get(y, key))  # noqa: E731
        p_sys, p_sys_dot = g(self.p_sys_key), float(context.get(ydot, self.p_sys_key))
        h, h_dot = g(self.key("h_he")), float(context.get(ydot, self.key("h_he")))
        p = g(self.key("p_fluid_cell"))
        rho = float(he.density(p_sys, h))
        m = np.array([g(self._pk(i, "m_dot")) for i in range(self.n_ports)])
        h_in = np.array([g(self._pk(i, "h_instream")) for i in range(self.n_ports)])
        q = float(self.heat_rate(t, y, context) if callable(self.heat_rate) else self.heat_rate)
        if self.steady_pressure_closure is not None:
            mass = p_sys - self.steady_pressure_closure
        else:
            mass = float(flow.mass_residual(self.V, he.drho_dp(p_sys, h), he.drho_dh(p_sys, h), p_sys_dot, h_dot,
                                            m.sum(), 0.0))
        inflow = np.asarray(flow.smooth_positive(m, self.m_eps), dtype=float)
        energy = rho * self.V * h_dot - self.V * p_sys_dot + float(np.dot(inflow, h - h_in)) - q
        parts = [mass, energy]
        if self.owns_p_sys:
            parts.append(p - p_sys)
        for i in range(self.n_ports):
            parts.append(g(self._pk(i, "p")) - (p + rho * pl.GRAVITY * (self.z_ref - self.z_ports[i])))
            parts.append(g(self._pk(i, "h_outflow")) - h)
        return np.asarray(parts, dtype=float)

    def residual_ids(self):
        if not self._ready:
            return None
        ids = [self.key("steady_pressure_closure" if self.steady_pressure_closure is not None else "mass"),
               self.key("energy")]
        if self.owns_p_sys:
            ids.append(self.key("pressure_reference"))
        for i in range(self.n_ports):
            ids += [self._pk(i, "pressure"), self._pk(i, "outflow_enthalpy")]
        return tuple(ids)

    def residual_scales(self):
        if not self._ready:
            return None
        cp = float(self.properties.cp(self.initial["p"], self.initial["T"]))
        scales = [1.0e5 if self.steady_pressure_closure is not None else self.m_nominal,
                  self.m_nominal * cp * 10.0]
        if self.owns_p_sys:
            scales.append(1.0e3)
        scales += [1.0e3, cp * 100.0] * self.n_ports
        return np.asarray(scales, dtype=float)

    # ------------------------------------------------------- 재고량·출력
    def helium_volumes(self):
        self._require_ready()
        return np.array([self.V])

    def helium_inventory(self, t, y, context):
        p_sys = float(context.get(y, self.p_sys_key))
        return np.array([float(self.properties.density(p_sys, context.get(y, self.key("h_he")))) * self.V])

    def outputs(self, t, y, context):
        p_sys = float(context.get(y, self.p_sys_key))
        h = float(context.get(y, self.key("h_he")))
        return {self.key("T_he"): float(self.properties.temperature(p_sys, h)),
                self.key("p"): float(context.get(y, self.key("p_fluid_cell")))}
