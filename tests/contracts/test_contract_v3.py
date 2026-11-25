from functools import partial

import numpy as np
import pytest

from interfaces.component import Component
from interfaces.contracts import (HeliumProperties, HeatTransferCorrelation,
                                  SimulationContext)
from interfaces.ports import FluidStreamPort
from interfaces.state_registry import StateRegistry
from models.core import Core
from models.heat_exchanger import HeatExchanger
from models.compressor import Compressor
from models.piping import Piping
from models.plenum import Plenum
from models.vessel import Vessel
from models.boundary import Boundary
from models.rccs import RCCS
from models.thermal_link import ThermalLink
from physics import heat_transfer
from solvers import DASSL
from solvers.jacobian import finite_difference_jacobian
from system.assembly import SystemAssembly
from system.connections import (StreamJunction,
                                ideal_fluid_stream_junction_residual)
from tests.fixtures.providers import (SyntheticHeliumProperties,
                                      SyntheticHeatCorrelation,
                                      SyntheticSGBoundary)
from tests.fixtures.closed_loop import make_loop
from tests.fixtures.synthetic_models import (SyntheticCoupling, SyntheticFluid,
                                             SyntheticHeatConstraint, SyntheticSolid,
                                             make_composite, make_derived, make_flat)


def test_legacy_public_imports_and_unimplemented_production_fail():
    for cls in (Core, HeatExchanger, Compressor, Piping, Plenum,
                Vessel, Boundary, RCCS):
        with pytest.raises(NotImplementedError):
            cls().register(StateRegistry())
    for function, args in ((heat_transfer.conduction, (400, 300, 1, 1, 1)),
                           (heat_transfer.convection, (400, 300, 1, 1)),
                           (heat_transfer.radiation, (400, 300, 1, 1, 1)),
                           (heat_transfer.heat_transfer_coefficient, (1, 1, 1))):
        with pytest.raises(NotImplementedError):
            function(*args)


def test_context_is_read_only_and_fixture_provider_broadcasts():
    registry = StateRegistry()
    context = SimulationContext(registry, config={"core": {"cells": [1, 2]}},
                                scenario={"events": []},
                                providers={"helium": SyntheticHeliumProperties()})
    with pytest.raises(TypeError):
        context.config["core"]["cells"] = (4,)
    assert context.config["core"]["cells"] == (1, 2)
    assert context.with_heat_rates({"link": [1, 2]}).heat_rates["link"] == (1, 2)
    assert not context.heat_rates
    provider = context.providers["helium"]
    assert isinstance(provider, HeliumProperties)
    assert np.asarray(provider.density(np.array([10.0, 20.0]), 5.0)).shape == (2,)
    assert np.allclose(provider.enthalpy(1.0, [2.0, 4.0]), [10.0, 20.0])
    assert isinstance(SyntheticHeatCorrelation(), HeatTransferCorrelation)
    assert SyntheticHeatCorrelation().coefficient({"T_wall": [1, 2]}, {"A": [1, 1]}).shape == (2,)
    assert SyntheticSGBoundary().heat_sink(0, 0.5) == 5.0
    bound = SystemAssembly([make_flat().components[0]], context=context)
    assert bound.context.registry is bound.registry


def test_registry_alias_resolves_after_registration_without_new_state():
    reg = StateRegistry()
    reg.alias("core.p", "p_sys")
    reg.register("p_sys", kind="differential", initial=5e6, unit="Pa", owner="system")
    reg.freeze()
    assert reg.size == 1
    assert reg.names == ("p_sys",)
    assert reg.get(reg.initial_state(), "core.p") == 5e6
    assert reg.spec("core.p").owner == "system"
    with pytest.raises(ValueError, match="Duplicate"):
        other = StateRegistry()
        other.register("p_sys", kind="algebraic", initial=1, unit="Pa")
        other.alias("p_sys", "p_sys")


def _stream_registry(ports):
    reg = StateRegistry()
    for port in ports:
        for key, unit in ((port.p, "Pa"), (port.m_dot, "kg/s"),
                          (port.h_outflow, "J/kg"), (port.h_instream, "J/kg")):
            reg.register(key, kind="algebraic", initial=0.0, unit=unit)
    reg.freeze()
    return reg


def test_junction_metadata_scalar_keys_media_and_alias_independence():
    ports = tuple(FluidStreamPort(f"{n}.p", f"{n}.m_dot", f"{n}.h_out", f"{n}.h_in")
                  for n in "abc")
    reg = _stream_registry(ports)
    junction = StreamJunction("mix", ports, 1e-6, 1e6, 10.0, 1e5)
    junction.validate(reg)
    assert len(junction.residual_ids()) == len(junction.residual_scales()) == 2 * len(ports)
    assert junction.residual_scales().tolist() == [1e6, 1e6, 10.0, 1e5, 1e5, 1e5]
    y = reg.initial_state()
    assert junction(y, reg).shape == (6,)
    different_medium = ports[:2] + (FluidStreamPort(*(
        ports[2].p, ports[2].m_dot, ports[2].h_outflow, ports[2].h_instream),
        medium="water"),)
    with pytest.raises(ValueError, match="different fluid media"):
        StreamJunction("bad", different_medium, 1e-6, 1e6, 10, 1e5).validate(reg)
    duplicate_pressure = (ports[0], FluidStreamPort(ports[0].p, ports[1].m_dot,
                                                     ports[1].h_outflow, ports[1].h_instream))
    with pytest.raises(ValueError, match="duplicate p"):
        StreamJunction("bad", duplicate_pressure, 1e-6, 1e6, 10, 1e5).validate(reg)


def test_missing_or_array_port_key_fails_at_assembly_construction():
    class BadPort(Component):
        def register(self, registry):
            registry.register("a.p", kind="algebraic", initial=1.0, unit="Pa")
            registry.register("a.m", kind="algebraic", initial=1.0, unit="kg/s")
            registry.register("a.h_out", kind="algebraic", initial=1.0, unit="J/kg",
                              shape=(2,))

        def residual(self, t, y, ydot, context):
            return np.zeros(4)

        def ports(self):
            return {"out": FluidStreamPort("a.p", "a.m", "a.h_out", "a.h_in")}

    with pytest.raises(ValueError, match="must be scalar"):
        SystemAssembly([BadPort()])


def test_flat_and_composite_dae_match_residual_jacobian_and_dassl():
    flat, composite = make_flat(), make_composite()
    assert flat.registry.names == composite.registry.names
    y = np.array([310.0, 315.0, 25.0])
    yp = np.array([0.1, -0.2, 0.0])
    y_copy, yp_copy = y.copy(), yp.copy()
    assert flat.residual_ids(0, y, yp) == composite.residual_ids(0, y, yp)
    assert np.array_equal(flat.residual_scales(0, y, yp), composite.residual_scales(0, y, yp))
    assert np.allclose(flat.residual(0, y, yp), composite.residual(0, y, yp))
    assert np.array_equal(y, y_copy) and np.array_equal(yp, yp_copy)
    displaced = y + np.array([3.0, -2.0, 1.0])
    composite.residual(0.5, displaced, yp)
    assert np.allclose(composite.residual(0, y, yp), flat.residual(0, y, yp))
    j_flat = finite_difference_jacobian(lambda z: flat.residual(0, z, yp), y)
    j_composite = finite_difference_jacobian(lambda z: composite.residual(0, z, yp), y)
    assert np.allclose(j_flat, j_composite, atol=1e-7)
    diag = composite.diagnose_jacobian(0, y, yp, cj=1.0)
    assert diag.rank == diag.size == 3
    assert np.isfinite(diag.condition)
    options = dict(dt=0.1, rtol=1e-6, atol=1e-8, max_order=3)
    def solve(system):
        return DASSL(residual_scale=system.residual_scales(), **options).solve(
            system.residual, t_span=(0, 2), y0=system.registry.initial_state(),
            differential_mask=system.registry.differential_mask(),
        )
    result_flat, result_composite = solve(flat), solve(composite)
    assert np.array_equal(result_flat.t, result_composite.t)
    assert np.allclose(result_flat.y, result_composite.y, atol=1e-7)
    assert np.allclose(result_flat.ydot, result_composite.ydot, atol=1e-7)


def test_core_can_register_stateful_power_with_its_own_residual():
    class StatefulPower:
        def register(self, registry):
            registry.register("test.power", kind="algebraic", unit="W",
                              initial=50.0, scale=50.0)

        def evaluate(self, t, y, context):
            return np.array([context.get(y, "test.power")])

        def residual(self, t, y, ydot, context):
            return np.array([context.get(y, "test.power") - 50.0])

        def residual_ids(self):
            return ("test.power_constraint",)

        def residual_scales(self):
            return np.array([50.0])

    core = Core(fluid=SyntheticFluid(), solid=SyntheticSolid(),
                power=StatefulPower(), coupling=SyntheticCoupling())
    system = SystemAssembly([core, SyntheticHeatConstraint()])
    assert system.registry.names == ("test.cold_T", "test.hot_T", "test.power", "test.q_dot")
    assert system.residual_ids() == ("test.cold_energy", "test.hot_energy",
                                     "test.power_constraint", "test.heat_link")
    assert system.diagnose_jacobian(0, system.registry.initial_state(),
                                    np.zeros(4), cj=1.0).rank == 4


def test_derived_heat_rate_matches_algebraic_temperature_trajectory():
    algebraic = make_composite()
    derived = make_derived()
    def solve(system):
        return DASSL(dt=0.1, rtol=1e-6, atol=1e-8, max_order=3,
                     residual_scale=system.residual_scales()).solve(
            system.residual, t_span=(0, 2), y0=system.registry.initial_state(),
            differential_mask=system.registry.differential_mask(),
        )
    a, b = solve(algebraic), solve(derived)
    assert np.allclose(np.interp(a.t, b.t, b.y[:, 0]), a.y[:, 0], atol=0.005)
    assert np.allclose(np.interp(a.t, b.t, b.y[:, 1]), a.y[:, 1], atol=0.005)


def test_thermal_link_conserves_different_grids_and_rejects_bad_transfer():
    transfer = ThermalLink("core_to_vessel", lambda t, y, context: ([-2.0, -3.0], [5.0]))
    heat = transfer.evaluate(0, np.zeros(1), None)
    assert np.sum(heat.to_a) + np.sum(heat.to_b) == 0
    assert not heat.to_a.flags.writeable
    with pytest.raises(ValueError, match="not conservative"):
        ThermalLink("bad", lambda t, y, c: ([-2.0], [3.0])).evaluate(0, [], None)


def test_square_system_can_still_have_singular_scaled_jacobian():
    class Redundant(Component):
        def register(self, registry):
            registry.register("test.a", kind="algebraic", unit="1", initial=1, scale=1)
            registry.register("test.b", kind="algebraic", unit="1", initial=1, scale=1)

        def residual(self, t, y, ydot, context):
            return np.array([y[0] + y[1], 2 * y[0] + 2 * y[1]])

    system = SystemAssembly([Redundant()])
    diagnosis = system.diagnose_jacobian(0, system.registry.initial_state(), np.zeros(2))
    assert diagnosis.rank == 1
    assert diagnosis.size == 2


def test_closed_loop_transient_and_steady_inventory_closure_ledgers():
    transient = make_loop()
    steady = make_loop(steady_inventory=True)
    y = np.array([10.0, 10.0, 8.0, 2.0, 2.0])
    ydot = np.zeros(5)
    assert transient.registry.names == steady.registry.names
    assert transient.diagnose_jacobian(0, y, ydot, cj=1.0).rank == 5
    assert transient.diagnose_jacobian(0, y, ydot, cj=0.0).rank == 4
    assert steady.diagnose_jacobian(0, y, ydot, cj=0.0).rank == 5
    assert np.allclose(steady.residual(0, y, ydot), 0.0)
    with pytest.raises(ValueError, match="not square"):
        make_loop(duplicate_mass=True).residual(0, y, ydot)
