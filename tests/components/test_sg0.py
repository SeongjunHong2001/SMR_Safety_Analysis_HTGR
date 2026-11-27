"""Role 4 SG-0 (일정 온도 물 열침원) 단독 시험: SGHeliumSide·NoStorageWall·FixedTemperatureWaterSink·SG0Coupling."""
import copy
from pathlib import Path

import numpy as np
import pytest
import yaml

from interfaces.component import Component
from models._inputs import ModelInputMissing
from models.heat_exchanger import HeatExchanger
from models.heat_exchanger.coupling import SG0Coupling
from models.heat_exchanger.helium_side import SGHeliumSide
from models.heat_exchanger.secondary_side import FixedTemperatureWaterSink
from models.heat_exchanger.wall import NoStorageWall
from physics.properties import make_helium_properties
from solvers import solve_steady_state
from system.assembly import SystemAssembly

ROOT = Path(__file__).resolve().parents[2]
GEO = yaml.safe_load((ROOT / "config/geometry.yaml").open(encoding="utf-8"))
OP = yaml.safe_load((ROOT / "config/operation.yaml").open(encoding="utf-8"))
SINK = OP["boundaries"]["sg_heat_sink"]
HE = make_helium_properties()
P0, M = 6.3167e6, 157.3


class Boundaries(Component):
    """시험 전용: p_sys 이완, SG He 입구 유량·온도, 출구 압력."""

    name = "bc"

    def __init__(self, ports, T_in):
        self.sg_ports, self.h_in = ports, float(HE.enthalpy(P0, T_in))

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=P0, scale=1e5, owner="primary_loop")

    def residual(self, t, y, ydot, context):
        g = lambda k: float(context.get(y, k))  # noqa: E731
        i, o = self.sg_ports["helium_inlet"], self.sg_ports["helium_outlet"]
        return np.array([float(context.get(ydot, "p_sys")) + g("p_sys") - P0, g(i.m_dot) - M, g(i.h_instream) - self.h_in,
                         g(o.p) - (P0 - 2.5e4), g(o.h_instream) - float(HE.enthalpy(P0, 528.15))])

    def residual_ids(self):
        return tuple(f"bc.{i}" for i in range(5))

    def residual_scales(self):
        return np.array([1e5, 100.0, 1e6, 1e3, 1e6])


def make_sg(**sink_kw):
    he_side = SGHeliumSide(GEO, HE, initial=dict(m_dot=M, T=750.0, p=P0))
    sink = FixedTemperatureWaterSink(SINK, helium_side=he_side, **sink_kw)
    return HeatExchanger(helium_side=he_side, wall=NoStorageWall(n_cells=he_side.n), secondary_side=sink,
                         coupling=SG0Coupling(he_side, sink), name="sg"), he_side, sink


def test_inputs_are_explicit():
    geo_null = copy.deepcopy(GEO)
    geo_null["heat_exchanger"]["helium_side"]["hydraulic_diameter_m"] = None
    with pytest.raises(ModelInputMissing, match="M11"):
        SGHeliumSide(geo_null, HE, initial=dict(m_dot=M, T=750.0, p=P0))     # D_h 없으면 명시 요구
    he_side = SGHeliumSide(GEO, HE, initial=dict(m_dot=M, T=750.0, p=P0))
    assert he_side.n == 9 and he_side.area_heat.sum() == pytest.approx(4020.0)
    # 1D 면적 = 다발 공극률 × 전면 9.23 m² (자유체적/길이), D_h = 4 V_자유/A_관 (R19 Table 5.3-8 [계산])
    assert he_side.pipe.A == pytest.approx(6.6622) and he_side.pipe.D_h == pytest.approx(0.05766)
    assert he_side.helium_volumes().sum() == pytest.approx(6.6622 * (-1.47 + 12.02), rel=1e-9)
    assert np.allclose(he_side.pipe.z_faces, GEO["heat_exchanger"]["helium_side"]["elevation_profile_m"])
    with pytest.raises(ValueError, match="exactly one"):
        FixedTemperatureWaterSink(SINK, helium_side=he_side)
    with pytest.raises(ValueError, match="phase_change"):
        FixedTemperatureWaterSink(dict(SINK, phase_change=True), helium_side=he_side, U=400.0)
    with pytest.raises(ValueError, match="no heat storage"):
        NoStorageWall(n_cells=2).residual(0, None, None, None, heat_rate=np.array([1.0e3, 0.0]))


def test_fixed_U_matches_discrete_analytic_solution():
    # 셀 값 열침원 + 1차 상류: T_i = (ṁ c_p T_{i−1} + U A_i T_w) / (ṁ c_p + U A_i)
    sg, he_side, sink = make_sg(U=421.5)
    system = SystemAssembly([sg, Boundaries(sg.ports(), 959.15)])
    y, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                              residual_scale=system.residual_scales())
    out = system.outputs(0.0, y)
    cp, U, Tw = HE.cp_const, 421.5, SINK["water_temperature_K"]
    T, expected = 959.15, []
    for A in he_side.area_heat:
        T = (M * cp * T + U * A * Tw) / (M * cp + U * A)
        expected.append(T)
    assert np.allclose(out["sg.helium.T_he"], expected, rtol=0, atol=1e-6)
    assert out["sg.secondary.Q_removed"] == pytest.approx(M * cp * (959.15 - expected[-1]), rel=1e-9)
    assert set(sg.ports()) == {"helium_inlet", "helium_outlet"}                 # SG-0: 물 포트 없음


def test_U_calibration_hits_target_outlet():
    sg, he_side, sink = make_sg(calibrate_outlet_K=528.15)
    system = SystemAssembly([sg, Boundaries(sg.ports(), 959.15)])
    y, _ = solve_steady_state(system.residual, 0.0, system.registry.initial_state(),
                              residual_scale=system.residual_scales())
    out = system.outputs(0.0, y)
    assert out["sg.helium.T_he"][-1] == pytest.approx(528.15, abs=1e-6)
    assert out["sg.secondary.Q_removed"] == pytest.approx(M * HE.cp_const * (959.15 - 528.15), rel=1e-9)  # ≈ 352.2 MW
    assert 350.0 < out["sg.secondary.U_overall"] < 600.0     # LMTD 추정 421.5와 같은 크기 (셀 값 공식이라 더 큼)


def test_coupling_is_conservative_and_signs():
    sg, he_side, sink = make_sg(U=400.0)
    system = SystemAssembly([sg, Boundaries(sg.ports(), 959.15)])
    y = system.registry.initial_state()
    q_he, q_wh, q_ws, q_sec = SG0Coupling(he_side, sink).evaluate(0.0, y, system.context)
    assert np.allclose(q_he + q_wh, 0) and np.allclose(q_ws + q_sec, 0) and np.all(q_he < 0)
