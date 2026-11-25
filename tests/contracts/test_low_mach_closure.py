"""Check the shared pressure/inventory ledger on a synthetic closed loop."""
import numpy as np
import pytest

from interfaces.contracts import HeliumInventory
from solvers import DASSL, solve_steady_state
from solvers.initialization import consistent_initialize
from system.inventory import total_helium_inventory
from tests.fixtures.low_mach_loop import make_low_mach_loop


def test_low_mach_pressure_and_inventory_contract_is_square_and_independent():
    transient = make_low_mach_loop()
    steady_inventory = make_low_mach_loop(steady=True)
    steady_pressure = make_low_mach_loop(steady=True, closure="pressure")
    y = transient.registry.initial_state()
    ydot = np.zeros_like(y)
    assert transient.registry.names == (
        "p_sys", "loop.h_he", "loop.p_fluid_cell", "loop.m_dot_fluid_face")
    assert transient.registry.spec("p_sys").owner == "primary_loop"
    assert transient.registry.differential_mask().tolist() == [True] * 3 + [False] * 4
    assert isinstance(transient.components[0], HeliumInventory)
    assert total_helium_inventory(transient.components, 0, y, transient.context) == 22.2
    with pytest.raises(ValueError, match="Duplicate helium inventory"):
        total_helium_inventory(transient.components * 2, 0, y, transient.context)

    changed_local_pressure = y.copy()
    changed_local_pressure[3:5] += 100.0
    assert np.array_equal(transient.components[0].density(y, transient.context),
                          transient.components[0].density(changed_local_pressure,
                                                          transient.context))
    changed_algebraic_derivative = ydot.copy()
    changed_algebraic_derivative[3:] = 50.0
    assert np.array_equal(transient.residual(0, y, ydot),
                          transient.residual(0, y, changed_algebraic_derivative))

    for system in (transient, steady_inventory, steady_pressure):
        assert system.residual(0, y, ydot).shape == (7,)
        assert len(system.residual_ids()) == len(system.residual_scales()) == 7
    assert transient.diagnose_jacobian(0, y, ydot, cj=1.0).rank == 7
    assert transient.diagnose_jacobian(0, y, ydot, cj=0.0).rank == 6
    assert steady_inventory.diagnose_jacobian(0, y, ydot, cj=0.0).rank == 7
    assert steady_pressure.diagnose_jacobian(0, y, ydot, cj=0.0).rank == 7
    with pytest.raises(ValueError, match="not square"):
        duplicate = make_low_mach_loop(duplicate_mass=True)
        duplicate.residual(0, duplicate.registry.initial_state(), np.zeros(7))


def test_low_mach_steady_closure_consistent_initialization_and_dassl():
    steady = make_low_mach_loop(steady=True)
    transient = make_low_mach_loop()
    guess = steady.registry.initial_state()
    guess += np.array([0.2, 0.1, -0.1, 0.3, -0.2, 0.1, -0.1])
    y0, ydot_steady = solve_steady_state(
        steady.residual, 0.0, guess, residual_scale=steady.residual_scales())
    assert np.allclose(ydot_steady, 0.0)
    assert np.linalg.norm(steady.residual(0.0, y0, ydot_steady), ord=np.inf) < 1e-7
    assert np.linalg.norm(transient.residual(0.0, y0, ydot_steady), ord=np.inf) < 1e-7
    assert abs(total_helium_inventory(steady.components, 0, y0, steady.context) - 22.2) < 1e-8

    initialized_y, initialized_ydot = consistent_initialize(
        transient.residual, 0.0, y0, transient.registry.differential_mask(),
        residual_scale=transient.residual_scales())
    assert np.allclose(initialized_y, y0, atol=1e-7)
    assert np.linalg.norm(transient.residual(0, initialized_y, initialized_ydot),
                          ord=np.inf) < 1e-7
    result = DASSL(dt=0.1, rtol=1e-6, atol=1e-8,
                   residual_scale=transient.residual_scales()).solve(
        transient.residual, t_span=(0.0, 2.0), y0=initialized_y,
        ydot0=initialized_ydot,
        differential_mask=transient.registry.differential_mask())
    assert result.t[-1] == 2.0
    inventories = np.array([transient.outputs(t, state)["loop.M_he_loop"]
                            for t, state in zip(result.t, result.y)])
    assert np.max(np.abs(inventories - inventories[0])) < 1e-6
    assert abs(result.y[-1, 1] - result.y[0, 1]) > 1e-3

    # u = h - p/rho: with the local mass equations, the enthalpy-form energy
    # residual must give d(sum rho*V*u)/dt = net external heat.
    for t, state, derivative in zip(result.t, result.y, result.ydot):
        p_sys, h = state[0], state[1:3]
        pdot, hdot = derivative[0], derivative[1:3]
        rho = p_sys + 0.1 * h
        rho_dot = pdot + 0.1 * hdot
        energy_rate = np.sum(rho_dot * h + rho * hdot - pdot)
        net_heat = 0.5 * np.sin(t) - 0.5 * (h[1] - 10.0)
        assert abs(energy_rate - net_heat) < 2e-5
