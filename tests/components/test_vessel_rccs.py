"""Role 4 RPV(Vessel)·RCCS-0·표면–유체 대류 링크 단독 시험."""
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from interfaces.contracts import SimulationContext
from interfaces.state_registry import StateRegistry
from models._inputs import ModelInputMissing
from models.rccs.rccs import RCCS0Sink
from models.vessel import Vessel
from models.vessel.vessel import axial_overlap_matrix, helium_gap_conductance
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN
from physics.heat_transfer.solid_properties import R1SolidProperties
from physics.properties import make_helium_properties

ROOT = Path(__file__).resolve().parents[2]
GEO = yaml.safe_load((ROOT / "config/geometry.yaml").open(encoding="utf-8"))
SP = R1SolidProperties()
HE = make_helium_properties()
Z = np.linspace(0.85, 13.59, 6)


def make_vessel(links=(("riser_rpv", "to_a"), ("rpv_rccs", "to_a"))):
    return Vessel(GEO["rpv"], SP, z_faces=Z, initial=dict(T=500.0), links=links)


def test_vessel_requires_inputs_and_exposes_aliased_surfaces():
    with pytest.raises(ModelInputMissing):
        Vessel().register(StateRegistry())
    with pytest.raises(ValueError, match="z_faces"):
        Vessel(GEO["rpv"], SP, initial=dict(T=500.0))
    v = make_vessel()
    ports = v.ports()
    assert ports["inner_surface[2]"].T == ports["outer_surface[2]"].T == "rpv.T"   # 한 노드 (별칭)
    assert ports["inner_surface[2]"].index == (2,)


def test_vessel_energy_residual_uses_link_rates_and_conserves_axial_conduction():
    v = make_vessel()
    reg = StateRegistry()
    v.register(reg)
    reg.freeze()
    y = reg.initial_state()
    y[:] = [450.0, 470.0, 500.0, 480.0, 460.0]
    ydot = np.full(5, 0.01)
    rates = {"riser_rpv.to_a": np.full(5, 2.0e4), "rpv_rccs.to_a": np.full(5, -1.5e4)}
    ctx = SimulationContext(reg, heat_rates=rates)
    r = v.residual(0.0, y, ydot, ctx)
    T = np.asarray(reg.get(y, "rpv.T"))
    expected = v.heat_capacity(T) * 0.01 - 5.0e3 - v.conduction_net(T)
    assert np.allclose(r, expected)
    assert abs(v.conduction_net(T).sum()) < 1e-9
    assert v.heat_capacity(np.array([500.0]))[0] == pytest.approx(
        7500 * 580 * math.pi * (3.405 ** 2 - 3.275 ** 2) * (13.59 - 0.85) / 5)
    with pytest.raises(ModelInputMissing, match="rpv_rccs.to_a"):
        v.residual(0.0, y, ydot, SimulationContext(reg, heat_rates={"riser_rpv.to_a": np.zeros(5)}))


def test_rccs0_sink_matches_hand_calculation():
    v = make_vessel()
    sink = RCCS0Sink(GEO["rccs"], v, SP, T_sink=303.15)
    # ε_eff: RPV 0.85, RCCS 0.80, r 3.405 / 4.633 → 0.735 (변수표 epsilon_rccs_effective 0.74)
    assert sink.eps_eff == pytest.approx(1 / (1 / 0.85 + (3.405 / 4.633) * (1 / 0.80 - 1)), rel=1e-12)
    assert sink.eps_eff == pytest.approx(0.735, abs=1e-3)
    T = np.array([500.0, 520.0, 540.0, 520.0, 500.0])
    dz = (13.59 - 0.85) / 5
    rad = sink.eps_eff * STEFAN_BOLTZMANN * 2 * math.pi * 3.405 * dz * (T ** 4 - 303.15 ** 4)
    air = 2 * math.pi * 0.068 * dz / math.log(4.633 / 3.405) * (T - 303.15)
    assert np.allclose(sink.heat_removed(T), rad + air, rtol=1e-12)
    assert np.all(sink.heat_removed(np.full(5, 303.15)) == 0.0)
    with pytest.raises(ModelInputMissing, match="boundary_temperature_K"):
        RCCS0Sink(GEO["rccs"], v, SP, T_sink=None)


def test_axial_overlap_matrix_is_conservative_between_grids():
    z_core = np.asarray(GEO["core"]["elevation_profile_m"], dtype=float)     # 내림차순 17셀
    W = axial_overlap_matrix(z_core, Z)                                     # 오름차순 5셀
    assert W.shape == (17, 5)
    assert np.allclose(W.sum(axis=1), 1.0)                                  # 노심 셀은 riser 고도 범위 안
    Q = np.linspace(1.0, 2.0, 17)
    assert (W.T @ Q).sum() == pytest.approx(Q.sum())                       # 유체 쪽 배분 합 보존
    W_part = axial_overlap_matrix([0.0, 2.0], [1.0, 3.0])
    assert W_part[0, 0] == pytest.approx(0.5)


def test_helium_gap_conductance_conduction_plus_radiation():
    g = helium_gap_conductance(HE, 3.201, 3.275, np.array([2.0]), 0.85, 0.85)
    T = 500.0
    k = float(HE.conductivity(6.39e6, T))
    eps = 1 / (1 / 0.85 + (3.201 / 3.275) * (1 / 0.85 - 1))
    expected = 2 * math.pi * k * 2.0 / math.log(3.275 / 3.201) + 4 * STEFAN_BOLTZMANN * eps * T ** 3 * 2 * math.pi * 3.201 * 2.0
    assert g(T, T, 6.39e6)[0] == pytest.approx(expected, rel=1e-12)
