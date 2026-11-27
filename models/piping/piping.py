"""배관 1D 질량·운동량·에너지 보존 (Role 2).

팀 보고서 [2] 파트 B §3의 1D FVM 잔차(physics.flow)와 파트 C·D의 압력손실(physics.pressure_loss)을 쓰는
두 포트 배관 Component다. 덕트, riser, 환형부처럼 단면이 일정한 유로를 N셀로 나눈다.

상태 (``{name}`` = Component 이름)
    {name}.h_he (N,) differential, {name}.p_fluid_cell (N,) algebraic, {name}.m_dot_fluid_face (N+1,) algebraic
    {name}.helium_inlet.{p, m_dot, h_outflow, h_instream}, {name}.helium_outlet.{...} (스칼라 algebraic)
잔차
    질량 N, 에너지 N, 운동량 N+1, 포트 4 (입구 유량 = 면 0, 입구 유출 h = 셀 0, 출구 유량 = −면 N, 출구 유출 h = 셀 N−1).
    포트의 나머지 식(압력 같음, 유입 stream)은 StreamJunction이 준다.
경계 면
    inlet/outlet = (종류, K_정방향, K_역류). 종류 "pipe"는 동압 변화 없이 K만, "plenum"은 플레넘 베르누이 ±q + K.
열
    기본 단열. heat_source(t, y, context, pipe) → (N,) W (유체로 들어오면 양수)를 주면 에너지식에 넣는다.
"""
from __future__ import annotations

import numpy as np

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from models._inputs import GEOMETRY, ModelInputMissing, require
from physics import flow
from physics import pressure_loss as pl


class Piping(Component):
    """두 포트 1D 배관. ``config``는 geometry.yaml ``flow_network.nodes.<이름>`` dict."""

    def __init__(self, config=None, properties=None, *, name="piping", initial=None, inlet=("pipe", 0.0, 0.0),
                 outlet=("pipe", 0.0, 0.0), heat_source=None, face_K=None, eps_mdot_rel=1.0e-4,
                 delta_inert=False, p_sys_key="p_sys"):
        self.config = dict(config or {})
        self.properties = properties
        self.name = name
        self.initial = None if initial is None else {k: float(v) for k, v in initial.items()}
        self.inlet, self.outlet = tuple(inlet), tuple(outlet)
        self.heat_source = heat_source
        self.face_K = face_K
        self.eps_mdot_rel = float(eps_mdot_rel)
        if not isinstance(delta_inert, (bool, np.bool_)):
            raise TypeError("delta_inert must be a bool")
        self.delta_inert = bool(delta_inert)
        self.p_sys_key = p_sys_key
        self._ready = False
        if self.config:
            self._setup()

    # -------------------------------------------------------------- 설정
    def _setup(self):
        where = f"flow_network.nodes.{self.name}"
        c = self.config
        n = c.get("cells")
        if n is None or isinstance(n, bool) or not isinstance(n, int) or n < 1:
            raise ModelInputMissing(f"{GEOMETRY}: {where}.cells must be a positive integer")
        A = require(c, "flow_area_m2", where)
        D = require(c, "hydraulic_diameter_m", where)
        L = require(c, "length_m", where)
        z_in = require(c, "z_inlet_m", where, positive=False)
        z_out = require(c, "z_outlet_m", where, positive=False)
        self.n = n
        self.A = np.full(n, A)
        self.D_h = np.full(n, D)
        self.dz = np.full(n, L / n)
        self.dz_elev = np.full(n, (z_out - z_in) / n)
        if abs(z_out - z_in) > L * (1 + 1e-9):
            raise ValueError(f"{GEOMETRY}: {where} elevation change exceeds length")
        self.V = self.A * self.dz
        self.roughness = float(c.get("roughness_m", 0.0) or 0.0)
        shape = c.get("shape", "circular")
        self.phi_s = pl.PHI_S_CIRCULAR if shape == "circular" else pl.PHI_S_PARALLEL_PLATES
        if self.initial is None:
            raise ValueError(f"Piping '{self.name}' requires initial=dict(m_dot, T, p)")
        for key in ("m_dot", "T", "p"):
            if key not in self.initial:
                raise ValueError(f"Piping '{self.name}': initial[{key!r}] is required")
        self.m_nominal = max(abs(self.initial["m_dot"]), 1e-6)
        self.m_eps = pl.mdot_epsilon(self.m_nominal, self.eps_mdot_rel)
        self._ready = True

    def _require_ready(self):
        if not self._ready:
            raise ModelInputMissing(f"Piping '{self.name}': {GEOMETRY} flow_network.nodes.<name> config "
                                    "and properties are required")
        if self.properties is None:
            raise ModelInputMissing(f"Piping '{self.name}': helium properties provider is required")

    def _pk(self, port, field):
        return self.key(f"{port}.{field}")

    # -------------------------------------------------------------- 등록
    def register(self, registry):
        self._require_ready()
        he, ini = self.properties, self.initial
        h0 = float(he.enthalpy(ini["p"], ini["T"]))
        registry.register(self.key("h_he"), kind="differential", unit="J/kg", shape=(self.n,), initial=h0,
                          scale=abs(h0))
        registry.register(self.key("p_fluid_cell"), kind="algebraic", unit="Pa", shape=(self.n,), initial=ini["p"],
                          scale=1.0e5)
        registry.register(self.key("m_dot_fluid_face"), kind=flow.face_flow_state_kind(self.delta_inert),
                          unit="kg/s", shape=(self.n + 1,), initial=ini["m_dot"], scale=self.m_nominal)
        for port, m in (("helium_inlet", ini["m_dot"]), ("helium_outlet", -ini["m_dot"])):
            registry.register(self._pk(port, "p"), kind="algebraic", unit="Pa", initial=ini["p"], scale=1.0e5)
            registry.register(self._pk(port, "m_dot"), kind="algebraic", unit="kg/s", initial=m, scale=self.m_nominal)
            registry.register(self._pk(port, "h_outflow"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))
            registry.register(self._pk(port, "h_instream"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))

    def ports(self):
        return {port: FluidStreamPort(self._pk(port, "p"), self._pk(port, "m_dot"), self._pk(port, "h_outflow"),
                                      self._pk(port, "h_instream"), medium="helium")
                for port in ("helium_inlet", "helium_outlet")}

    # -------------------------------------------------------------- 잔차
    def _state(self, y, ydot, context):
        he = self.properties
        p_sys = float(context.get(y, self.p_sys_key))
        h = np.asarray(context.get(y, self.key("h_he")), dtype=float)
        T = np.asarray(he.temperature(p_sys, h), dtype=float)
        return dict(p_sys=p_sys, p_sys_dot=float(context.get(ydot, self.p_sys_key)), h=h, T=T,
                    h_dot=np.asarray(context.get(ydot, self.key("h_he")), dtype=float),
                    p=np.asarray(context.get(y, self.key("p_fluid_cell")), dtype=float),
                    m=np.asarray(context.get(y, self.key("m_dot_fluid_face")), dtype=float),
                    m_dot_dot=np.asarray(context.get(ydot, self.key("m_dot_fluid_face")), dtype=float),
                    rho=np.asarray(he.density(p_sys, h), dtype=float),
                    drho_dp=np.asarray(he.drho_dp(p_sys, h), dtype=float),
                    drho_dh=np.asarray(he.drho_dh(p_sys, h), dtype=float),
                    mu=np.asarray(he.viscosity(p_sys, T), dtype=float))

    def wall_heat(self, t, y, context):
        if self.heat_source is None:
            return np.zeros(self.n)
        q = np.asarray(self.heat_source(t, y, context, self), dtype=float)
        if q.shape != (self.n,):
            raise ValueError(f"Piping '{self.name}': heat_source must return shape ({self.n},)")
        return q

    def residual(self, t, y, ydot, context):
        self._require_ready()
        g = lambda key: float(context.get(y, key))  # noqa: E731
        s = self._state(y, ydot, context)
        q = self.wall_heat(t, y, context)
        cell = flow.cell_residuals(V=self.V, rho=s["rho"], drho_dp=s["drho_dp"], drho_dh=s["drho_dh"], h=s["h"],
                                   h_dot=s["h_dot"], p_sys_dot=s["p_sys_dot"], q_wall=q, m_faces=s["m"],
                                   h_inlet_boundary=g(self._pk("helium_inlet", "h_instream")),
                                   h_outlet_boundary=g(self._pk("helium_outlet", "h_instream")), m_eps=self.m_eps)
        terms = pl.channel_face_terms(s["m"], s["rho"], s["mu"], self.A, self.D_h, self.dz, self.dz_elev,
                                      m_eps=self.m_eps, roughness=self.roughness, phi_s=self.phi_s,
                                      inlet=self.inlet, outlet=self.outlet, face_K=self.face_K)
        p_left = np.concatenate(([g(self._pk("helium_inlet", "p"))], s["p"]))
        p_right = np.concatenate((s["p"], [g(self._pk("helium_outlet", "p"))]))
        l_over_A = flow.face_inertia_length(np.concatenate(([0.0], self.dz)), np.concatenate(([1.0], self.A)),
                                            np.concatenate((self.dz, [0.0])), np.concatenate((self.A, [1.0])))
        mom = flow.momentum_residual(p_left, p_right, dp_acc=terms["acc"], dp_fric=terms["fric"],
                                     dp_form=terms["form"], dp_grav=terms["grav"],
                                     include_inertia=self.delta_inert, l_over_A=l_over_A, m_dot_dot=s["m_dot_dot"])
        ports = np.array([g(self._pk("helium_inlet", "m_dot")) - s["m"][0],
                          g(self._pk("helium_inlet", "h_outflow")) - s["h"][0],
                          g(self._pk("helium_outlet", "m_dot")) + s["m"][-1],
                          g(self._pk("helium_outlet", "h_outflow")) - s["h"][-1]])
        return np.concatenate((cell["R_mass"], cell["R_energy"], np.asarray(mom, dtype=float), ports))

    def residual_ids(self):
        if not self._ready:
            return None
        return tuple([self.key(f"mass[{i}]") for i in range(self.n)]
                     + [self.key(f"energy[{i}]") for i in range(self.n)]
                     + [self.key(f"momentum[{j}]") for j in range(self.n + 1)]
                     + [self._pk("helium_inlet", "flow"), self._pk("helium_inlet", "outflow_enthalpy"),
                        self._pk("helium_outlet", "flow"), self._pk("helium_outlet", "outflow_enthalpy")])

    def residual_scales(self):
        if not self._ready:
            return None
        cp = float(self.properties.cp(self.initial["p"], self.initial["T"]))
        m = self.m_nominal
        return np.concatenate((np.full(self.n, m), np.full(self.n, m * cp * 10.0), np.full(self.n + 1, 1.0e3),
                               [m, cp * 100.0, m, cp * 100.0]))

    # ------------------------------------------------------- 재고량·출력
    def helium_volumes(self):
        self._require_ready()
        return self.V.copy()

    def helium_inventory(self, t, y, context):
        p_sys = float(context.get(y, self.p_sys_key))
        return np.asarray(self.properties.density(p_sys, context.get(y, self.key("h_he"))), dtype=float) * self.V

    def outputs(self, t, y, context):
        p_sys = float(context.get(y, self.p_sys_key))
        h = np.asarray(context.get(y, self.key("h_he")), dtype=float)
        m = np.asarray(context.get(y, self.key("m_dot_fluid_face")), dtype=float)
        return {self.key("T_he"): np.asarray(self.properties.temperature(p_sys, h), dtype=float),
                self.key("m_dot_inlet"): float(m[0]),
                self.key("dp_inlet_to_outlet"): float(context.get(y, self._pk("helium_inlet", "p"))
                                                      - context.get(y, self._pk("helium_outlet", "p"))),
                self.key("Q_wall"): float(np.sum(self.wall_heat(t, y, context)))}
