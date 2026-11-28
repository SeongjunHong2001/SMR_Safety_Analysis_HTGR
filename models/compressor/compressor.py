"""Compressor(주순환기) 압력상승과 헬륨 가열 (Role 2).

체적 없는 두 포트 장치. 포트는 ``helium_inlet``(흡입), ``helium_outlet``(토출)이다.

모드 (``mode``)
    "prescribed_flow"  S0 정상 초기화용 (팀 결정 2026-11-26: 성능곡선 미확보 → 유량 지정, 승압은 결과).
                       미지수 dp_comp (algebraic)와 유량식 ṁ_in − ṁ_set = 0. 과도 사고 해석에 쓰지 않는다.
    "performance"      승압 = performance.evaluate(speed=ω(t), mass_flow, inlet_state, geometry)["pressure_rise"].
                       ``performance``(CompressorPerformance)와 ``speed_schedule``(t → rad/s)을 반드시 준다.
                       유량은 루프 운동량으로 정해진다 (신호로 유량을 강제하지 않음).
    "coastdown_inertia" 회전자 관성 ODE (mode_coastdown = inertia_ode). 관성 I_circ_rotor가 미확보라 아직 지원하지 않는다.
잔차
    질량 ṁ_in + ṁ_out = 0, 운동량 p_out − p_in − ΔP_comp = 0,
    토출 유출 엔탈피 h_out,outlet − (h_instream,inlet + w Δh) = 0, 흡입 유출 엔탈피 h_out,inlet − (h_instream,outlet + (1 − w) Δh) = 0
    (w: 유동 방향 가중), prescribed_flow면 유량식 1개 추가.
헬륨 가열 [가정]
    Q = Q_rated (|ṁ|_sm / ṁ_rated)³ (상사법칙, 같은 밀도), Δh = Q / |ṁ|_sm = Q_rated |ṁ|_sm² / ṁ_rated³.
    정격 Q_rated = operation.yaml compressor.rated_heat_to_helium_W (3.1 MW, R19·R22).
"""
from __future__ import annotations

import numpy as np

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from models._inputs import ModelInputMissing
from physics import pressure_loss as pl

MODES = ("prescribed_flow", "performance", "coastdown_inertia")


class Compressor(Component):
    """주순환기. ``config``는 operation.yaml ``compressor`` dict (정격값)."""

    def __init__(self, config=None, properties=None, *, name="main_circulator", mode=None, m_dot_set=None,
                 m_rated=None, performance=None, speed_schedule=None, geometry=None, initial=None,
                 eps_mdot_rel=1.0e-4, p_sys_key="p_sys"):
        self.config = dict(config or {})
        self.properties = properties
        self.name = name
        self.mode = mode
        self.m_dot_set = m_dot_set
        self.m_rated = m_rated
        self.performance = performance
        self.speed_schedule = speed_schedule
        self.geometry = dict(geometry or {})
        self.initial = None if initial is None else {k: float(v) for k, v in initial.items()}
        self.eps_mdot_rel = float(eps_mdot_rel)
        self.p_sys_key = p_sys_key

    def _check(self):
        if self.mode is None or not self.config or self.properties is None:
            raise ModelInputMissing(f"Compressor '{self.name}': mode, config/operation.yaml compressor section "
                                    "and properties are required")
        if self.mode not in MODES:
            raise ValueError(f"Compressor mode must be one of {MODES}")
        if self.mode == "coastdown_inertia":
            raise ModelInputMissing("Compressor mode 'coastdown_inertia' needs the rotor inertia I_circ_rotor "
                                    "(config/geometry.yaml compressor.rotational_inertia_kg_m2 is null, M09)")
        if self.mode == "prescribed_flow" and self.m_dot_set is None:
            raise ValueError("Compressor mode 'prescribed_flow' requires m_dot_set (kg/s or callable t → kg/s)")
        if self.mode == "performance" and (self.performance is None or self.speed_schedule is None):
            raise ModelInputMissing("Compressor mode 'performance' requires a performance map provider "
                                    "(M09, no validated curve yet) and speed_schedule")
        if self.m_rated is None or not self.m_rated > 0:
            raise ValueError("Compressor requires m_rated (kg/s at which rated_heat_to_helium_W applies)")
        if self.initial is None or any(k not in self.initial for k in ("m_dot", "T", "p")):
            raise ValueError(f"Compressor '{self.name}' requires initial=dict(m_dot, T, p)")
        q = self.config.get("rated_heat_to_helium_W")
        if q is None:
            raise ModelInputMissing("config/operation.yaml: missing compressor.rated_heat_to_helium_W")
        self.Q_rated = float(q)
        self.m_eps = pl.mdot_epsilon(self.m_rated, self.eps_mdot_rel)

    def _pk(self, port, field):
        return self.key(f"{port}.{field}")

    # -------------------------------------------------------------- 등록
    def register(self, registry):
        self._check()
        he, ini = self.properties, self.initial
        h0 = float(he.enthalpy(ini["p"], ini["T"]))
        for port, m in (("helium_inlet", ini["m_dot"]), ("helium_outlet", -ini["m_dot"])):
            registry.register(self._pk(port, "p"), kind="algebraic", unit="Pa", initial=ini["p"], scale=1.0e5)
            registry.register(self._pk(port, "m_dot"), kind="algebraic", unit="kg/s", initial=m,
                              scale=float(self.m_rated))
            registry.register(self._pk(port, "h_outflow"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))
            registry.register(self._pk(port, "h_instream"), kind="algebraic", unit="J/kg", initial=h0, scale=abs(h0))
        if self.mode == "prescribed_flow":
            dp0 = float(self.config.get("rated_pressure_rise_Pa") or 1.0e4)
            registry.register(self.key("dp_comp"), kind="algebraic", unit="Pa", initial=dp0, scale=1.0e4)

    def ports(self):
        return {port: FluidStreamPort(self._pk(port, "p"), self._pk(port, "m_dot"), self._pk(port, "h_outflow"),
                                      self._pk(port, "h_instream"), medium="helium")
                for port in ("helium_inlet", "helium_outlet")}

    # -------------------------------------------------------------- 잔차
    def heat_to_helium(self, m_dot):
        """Q = Q_rated (|ṁ|_sm/ṁ_rated)³ [W] 와 Δh = Q/|ṁ|_sm [J/kg]."""
        m_abs = float(pl.smooth_abs(m_dot, self.m_eps))
        dh = self.Q_rated * m_abs * m_abs / self.m_rated ** 3
        return dh * m_abs, dh

    def pressure_rise(self, t, y, context):
        if self.mode == "prescribed_flow":
            return float(context.get(y, self.key("dp_comp")))
        he = self.properties
        p_sys = float(context.get(y, self.p_sys_key))
        h_in = float(context.get(y, self._pk("helium_inlet", "h_instream")))
        T_in = float(he.temperature(p_sys, h_in))
        inlet_state = dict(p=p_sys, h=h_in, T=T_in, rho=float(he.density(p_sys, h_in)))
        result = self.performance.evaluate(speed=float(self.speed_schedule(t)),
                                           mass_flow=float(context.get(y, self._pk("helium_inlet", "m_dot"))),
                                           inlet_state=inlet_state, geometry=self.geometry)
        return float(result["pressure_rise"])

    def residual(self, t, y, ydot, context):
        self._check()
        g = lambda key: float(context.get(y, key))  # noqa: E731
        m_in, m_out = g(self._pk("helium_inlet", "m_dot")), g(self._pk("helium_outlet", "m_dot"))
        _, dh = self.heat_to_helium(m_in)
        w = float(pl.direction_weight(m_in, self.m_eps))
        r = [m_in + m_out,
             g(self._pk("helium_outlet", "p")) - g(self._pk("helium_inlet", "p")) - self.pressure_rise(t, y, context),
             g(self._pk("helium_outlet", "h_outflow")) - (g(self._pk("helium_inlet", "h_instream")) + w * dh),
             g(self._pk("helium_inlet", "h_outflow")) - (g(self._pk("helium_outlet", "h_instream")) + (1.0 - w) * dh)]
        if self.mode == "prescribed_flow":
            target = self.m_dot_set(t) if callable(self.m_dot_set) else self.m_dot_set
            r.append(m_in - float(target))
        return np.asarray(r, dtype=float)

    def residual_ids(self):
        if self.mode not in MODES:
            return None
        ids = [self.key("mass"), self.key("momentum"), self._pk("helium_outlet", "outflow_enthalpy"),
               self._pk("helium_inlet", "outflow_enthalpy")]
        if self.mode == "prescribed_flow":
            ids.append(self.key("prescribed_flow"))
        return tuple(ids)

    def residual_scales(self):
        if self.mode not in MODES or self.m_rated is None:
            return None
        m = float(self.m_rated)
        cp = (float(self.properties.cp(self.initial["p"], self.initial["T"]))
              if self.properties is not None and self.initial is not None else 5195.0)
        s = [m, 1.0e3, cp * 10.0, cp * 10.0]
        if self.mode == "prescribed_flow":
            s.append(m)
        return np.asarray(s, dtype=float)

    def outputs(self, t, y, context):
        m_in = float(context.get(y, self._pk("helium_inlet", "m_dot")))
        q, _ = self.heat_to_helium(m_in)
        return {self.key("pressure_rise"): float(context.get(y, self._pk("helium_outlet", "p"))
                                                 - context.get(y, self._pk("helium_inlet", "p"))),
                self.key("m_dot"): m_in, self.key("Q_to_helium"): q}
