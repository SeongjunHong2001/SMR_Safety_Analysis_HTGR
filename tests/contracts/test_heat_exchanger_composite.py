"""Assembly contract for an injected SG, using synthetic physical values only."""
import numpy as np
import pytest

from interfaces.state_registry import StateRegistry
from models.heat_exchanger import HeatExchanger
from solvers import DASSL
from solvers.jacobian import finite_difference_jacobian
from system.assembly import SystemAssembly
from tests.fixtures.synthetic_heat_exchanger import (make_composite, make_flat,
                                                     make_parts, SyntheticSGCoupling)


def test_heat_exchanger_needs_injected_parts_and_registers_each_once():
    with pytest.raises(NotImplementedError, match="helium_side.py"):
        HeatExchanger().register(StateRegistry())

    helium, wall, secondary = make_parts()
    model = HeatExchanger(helium_side=helium, wall=wall,
                          secondary_side=secondary, coupling=SyntheticSGCoupling())
    with pytest.raises(ValueError, match="Duplicate state name"):
        SystemAssembly([model, helium])


def test_heat_exchanger_composite_matches_flat_dae_and_integrates():
    composite, flat = make_composite(), make_flat()
    assert composite.registry.names == flat.registry.names
    assert composite.registry.size == 19
    assert composite.registry.differential_mask().sum() == 3
    ports = composite.components[0].ports()
    assert tuple(ports) == ("helium_inlet", "helium_outlet",
                            "secondary_inlet", "secondary_outlet")
    assert [ports[name].medium for name in ports] == ["helium", "helium", "water", "water"]

    y = composite.registry.initial_state()
    ydot = np.zeros_like(y)
    r_comp, ids_comp, scales_comp = composite.residual_details(0.0, y, ydot)
    r_flat, ids_flat, scales_flat = flat.residual_details(0.0, y, ydot)
    assert r_comp.shape == (composite.registry.size,)
    assert ids_comp == ids_flat
    assert np.array_equal(scales_comp, scales_flat)
    assert np.allclose(r_comp, r_flat)
    assert ids_comp[0] == "sg.helium.energy"
    assert ids_comp[9] == "sg.wall.energy"
    assert ids_comp[10] == "sg.secondary.energy"
    assert scales_comp[[0, 9, 10]].tolist() == [100.0, 100.0, 100.0]
    assert composite.outputs(0.0, y) == {
        "sg.helium.temperature": 350.0,
        "sg.wall.temperature": 320.0,
        "sg.secondary.temperature": 300.0,
    }
    j_comp = finite_difference_jacobian(lambda z: composite.residual(0.0, z, ydot), y)
    j_flat = finite_difference_jacobian(lambda z: flat.residual(0.0, z, ydot), y)
    assert np.allclose(j_comp, j_flat, atol=1e-8)
    assert composite.diagnose_jacobian(0.0, y, ydot, cj=1.0).rank == y.size

    def solve(system):
        return DASSL(dt=0.1, max_order=3, rtol=1e-6, atol=1e-8,
                     residual_scale=system.residual_scales()).solve(
            system.residual, t_span=(0.0, 1.0),
            y0=system.registry.initial_state(),
            differential_mask=system.registry.differential_mask(),
        )

    a, b = solve(composite), solve(flat)
    assert a.t[-1] == b.t[-1] == 1.0
    assert np.allclose(a.t, b.t, atol=1e-7)
    assert np.allclose(a.y, b.y, atol=1e-7)
    assert np.allclose(a.ydot, b.ydot, atol=1e-7)
    energy = 2.0 * a.y[:, 0] + 3.0 * a.y[:, 9] + 4.0 * a.y[:, 10]
    assert np.max(np.abs(energy - energy[0])) < 1e-5
    assert np.max(np.abs(a.y[:, 0] - a.y[0, 0])) > 0.0


def test_heat_exchanger_rejects_nonconservative_coupling():
    helium, wall, secondary = make_parts()

    class BadCoupling:
        def evaluate(self, t, y, context):
            return ([-1.0], [2.0], [-1.0], [1.0])

    system = SystemAssembly([HeatExchanger(helium_side=helium, wall=wall,
                                          secondary_side=secondary, coupling=BadCoupling())])
    with pytest.raises(ValueError, match="helium/wall.*not conservative"):
        system.residual(0.0, system.registry.initial_state(),
                        np.zeros(system.registry.size))
