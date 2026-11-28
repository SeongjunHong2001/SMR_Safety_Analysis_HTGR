"""RPV 측벽 열용량과 열 경로 (Role 4).

``Vessel``  RPV 측벽을 축방향 N셀(반경 방향 1노드)로 둔다. 상태 ``{name}.T`` (N,) differential.
            내·외면 온도는 같은 노드라 포트 ``inner_surface[i]``, ``outer_surface[i]``가 같은 키를 가리킨다
            (별칭, 저장항 중복 없음). 열은 ThermalLink 열률(context.heat_rates)로만 받는다.
            에너지식: ρ c_p V dT/dt − Σ 링크 열률 − 축방향 전도 = 0.
            상·하부 헤드(HS-09, HS-10)는 이 모델에 없다 (단열) [가정].
``SurfaceFluidConvection``  고체 표면 온도 배열과 1D 배관 유체(Piping) 사이 대류 링크.
            축방향 격자가 달라도 고도 겹침 비율로 대응시키고(보존), 대류에 고체 측 직렬 저항(반두께 전도,
            헬륨 간극 전도 + 복사)을 더할 수 있다. 배럴 → riser 헬륨, riser 헬륨 → RPV 내면에 쓴다.
"""
from __future__ import annotations

import math

import numpy as np

from interfaces.component import Component
from interfaces.ports import ThermalSurfacePort
from models._inputs import GEOMETRY, ModelInputMissing, require
from models.thermal_link import ThermalLink
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN, effective_emissivity_gray_enclosure


def axial_overlap_matrix(z_a, z_b):
    """W[i, j] = (셀 a_i와 셀 b_j가 겹치는 고도 길이) / (셀 a_i 길이). 면 배열 순서(오름·내림)는 상관없다."""
    z_a, z_b = np.asarray(z_a, dtype=float), np.asarray(z_b, dtype=float)
    lo_a, hi_a = np.minimum(z_a[:-1], z_a[1:]), np.maximum(z_a[:-1], z_a[1:])
    lo_b, hi_b = np.minimum(z_b[:-1], z_b[1:]), np.maximum(z_b[:-1], z_b[1:])
    overlap = np.clip(np.minimum(hi_a[:, None], hi_b[None, :]) - np.maximum(lo_a[:, None], lo_b[None, :]), 0.0, None)
    return overlap / (hi_a - lo_a)[:, None]


class Vessel(Component):
    """RPV 측벽. ``config``는 geometry.yaml ``rpv`` dict."""

    def __init__(self, config=None, properties=None, *, name="rpv", z_faces=None, initial=None, links=(),
                 axial_conduction=True):
        self.config = dict(config or {})
        self.properties = properties
        self.name = name
        self.links = tuple(links)              # (링크 이름, "to_a" | "to_b")
        self.axial_conduction = bool(axial_conduction)
        self._ready = False
        if self.config:
            where = "rpv"
            self.r_i = 0.5 * require(self.config, "inner_diameter_m", where)
            self.r_o = 0.5 * require(self.config, "outer_diameter_m", where)
            if self.r_o <= self.r_i:
                raise ValueError(f"{GEOMETRY}: rpv outer diameter must exceed inner diameter")
            if z_faces is None:
                raise ValueError("Vessel requires z_faces (축방향 셀 경계 고도) for the modeled side wall")
            self.z_faces = np.asarray(z_faces, dtype=float)
            self.dz = np.abs(np.diff(self.z_faces))
            if self.dz.size < 1 or np.any(self.dz <= 0):
                raise ValueError("Vessel z_faces must be strictly monotone with at least two entries")
            self.n = self.dz.size
            if initial is None or "T" not in initial:
                raise ValueError("Vessel requires initial=dict(T=K)")
            self.T0 = np.broadcast_to(np.asarray(initial["T"], dtype=float), (self.n,)).copy()
            for link, side in self.links:
                if side not in ("to_a", "to_b"):
                    raise ValueError("Vessel link side must be 'to_a' or 'to_b'")
            self._ready = True

    def _require(self):
        if not self._ready or self.properties is None:
            raise ModelInputMissing(f"Vessel '{self.name}': {GEOMETRY} rpv config, z_faces, initial and solid "
                                    "properties are required")

    @property
    def T_key(self):
        return self.key("T")

    def register(self, registry):
        self._require()
        registry.register(self.T_key, kind="differential", unit="K", shape=(self.n,), initial=self.T0, scale=1000.0)

    def ports(self):
        if not self._ready:
            return {}
        ports = {}
        for i in range(self.n):
            ports[f"inner_surface[{i}]"] = ThermalSurfacePort(self.T_key, index=i)
            ports[f"outer_surface[{i}]"] = ThermalSurfacePort(self.T_key, index=i)
        return ports

    def heat_capacity(self, T):
        p = self.properties
        V = math.pi * (self.r_o ** 2 - self.r_i ** 2) * self.dz
        return p.density("pressure_vessel") * np.asarray(p.cp("pressure_vessel", T), dtype=float) * V

    def conduction_net(self, T):
        net = np.zeros(self.n)
        if not self.axial_conduction or self.n < 2:
            return net
        A = math.pi * (self.r_o ** 2 - self.r_i ** 2)
        k = self.properties.conductivity("pressure_vessel")
        G = 1.0 / (self.dz[:-1] / (2 * k * A) + self.dz[1:] / (2 * k * A))
        q = G * (T[:-1] - T[1:])
        net[:-1] -= q
        net[1:] += q
        return net

    def link_heat(self, context):
        total = np.zeros(self.n)
        rates = getattr(context, "heat_rates", None)
        for link, side in self.links:
            key = f"{link}.{side}"
            if rates is None or key not in rates:
                raise ModelInputMissing(f"Vessel '{self.name}': heat rate {key!r} is missing from the context")
            rate = np.asarray(rates[key], dtype=float)
            if rate.shape != (self.n,):
                raise ValueError(f"Vessel '{self.name}': link {key} must give shape ({self.n},)")
            total += rate
        return total

    def residual(self, t, y, ydot, context):
        self._require()
        T = np.asarray(context.get(y, self.T_key), dtype=float)
        T_dot = np.asarray(context.get(ydot, self.T_key), dtype=float)
        return self.heat_capacity(T) * T_dot - self.link_heat(context) - self.conduction_net(T)

    def residual_ids(self):
        return None if not self._ready else tuple(self.key(f"energy[{i}]") for i in range(self.n))

    def residual_scales(self):
        return None if not self._ready else np.full(self.n, 1.0e4)

    def outputs(self, t, y, context):
        T = np.asarray(context.get(y, self.T_key), dtype=float)
        out = {self.key("T"): T.copy(), self.key("T_max"): float(T.max())}
        rates = getattr(context, "heat_rates", None) or {}
        for link, side in self.links:
            key = f"{link}.{side}"
            if key in rates:
                out[self.key(f"Q_from_{link}")] = float(np.sum(rates[key]))
        return out


class SurfaceFluidConvection:
    """고체 표면(배열 온도) ↔ Piping 유체 대류 ThermalLink 전달 함수 (side a = 고체, side b = 유체).

    solid_T(y, context) → (n_s,) K, z_solid (n_s+1,) 고도, area (n_s,) m² (대류 면적, 고체 셀 기준)
    pipe: Piping (coupling_state, z_faces), correlation: HeatTransferCorrelation
    series(T_s, T_f, p) → (n_s,) W/K 고체 쪽 직렬 컨덕턴스 (None이면 대류만)
    Q_i = (T_s,i − T_f,i) G_h G_s/(G_h + G_s), T_f,i = Σ_j W_ij T_j, 유체 셀 j 수열 = Σ_i W_ij Q_i (보존)
    h는 벽온도 두 번 평가 (T_wall = T_s → 직렬 분배 T_w → h).
    """

    def __init__(self, name, *, solid_T, z_solid, area, pipe, correlation, series=None):
        self.name = name
        self.solid_T = solid_T
        self.area = np.asarray(area, dtype=float)
        self.pipe = pipe
        self.corr = correlation
        self.series = series
        self.W = axial_overlap_matrix(z_solid, pipe.z_faces)
        coverage = self.W.sum(axis=1)
        if not np.allclose(coverage, 1.0, atol=1e-9):
            raise ValueError(f"{name}: every solid cell must lie inside the pipe's elevation range "
                             f"(coverage {coverage.min():.3f}–{coverage.max():.3f})")
        if self.area.shape != (self.W.shape[0],):
            raise ValueError(f"{name}: area must have one entry per solid cell")

    def link(self):
        return ThermalLink(self.name, self.transfer)

    def _h(self, fs, T_wall):
        W = self.W
        state = {key: W @ fs[key] for key in ("m_dot", "T_bulk", "mu", "k", "cp", "rho")}
        state["T_wall"] = T_wall
        geometry = dict(D_h=W @ fs["D_h"], A_flow=W @ fs["A_flow"], shape=fs["shape"][0])
        return np.asarray(self.corr.coefficient(state, geometry), dtype=float), state["T_bulk"]

    def transfer(self, t, y, context):
        fs = self.pipe.coupling_state(t, y, context)
        T_s = np.asarray(self.solid_T(y, context), dtype=float)
        h0, T_f = self._h(fs, T_s)
        G_s = None if self.series is None else np.asarray(self.series(T_s, T_f, fs["p_sys"]), dtype=float)
        if G_s is None:
            T_w = T_s
        else:
            G0 = h0 * self.area
            T_w = (G_s * T_s + G0 * T_f) / (G_s + G0)
        h, _ = self._h(fs, T_w)
        G_h = h * self.area
        G = G_h if G_s is None else G_h * G_s / (G_h + G_s)
        Q = G * (T_s - T_f)                       # 고체 → 유체 양수
        return -Q, self.W.T @ Q


def helium_gap_conductance(helium, r_in, r_out, dz, eps_in, eps_out):
    """동심 원통 헬륨 간극(전도 + 선형화 복사) 컨덕턴스 함수 (T_a, T_b, p) → W/K (셀별)."""
    eps = float(effective_emissivity_gray_enclosure(eps_in, eps_out, r_in / r_out))
    dz = np.asarray(dz, dtype=float)

    def conductance(T_a, T_b, p):
        T_m = 0.5 * (np.asarray(T_a, dtype=float) + np.asarray(T_b, dtype=float))
        k = np.asarray(helium.conductivity(p, T_m), dtype=float)
        return (2.0 * math.pi * k * dz / math.log(r_out / r_in)
                + 4.0 * STEFAN_BOLTZMANN * eps * T_m ** 3 * 2.0 * math.pi * r_in * dz)
    return conductance
