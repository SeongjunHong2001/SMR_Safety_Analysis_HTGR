from functools import partial

import numpy as np
import pytest
from interfaces.component import Component
from interfaces.state_registry import StateRegistry
from interfaces.ports import FluidPort, FluidStreamPort, ThermalPort
from solvers import DASSL, solve_steady_state
from system.assembly import SystemAssembly
from system.connections import (
    ideal_fluid_two_port_residual,
    ideal_fluid_stream_junction_residual,
    ideal_thermal_two_port_residual,
)


def test_variable_registry_indices_and_duplicate_rejection():
    reg = StateRegistry()
    reg.register('core.T', kind='differential', initial=300., unit='K', shape=(3,))
    reg.register('link.m_dot', kind='algebraic', initial=0., unit='kg/s')
    assert reg.size == 4
    assert np.allclose(reg.get(reg.initial_state(), 'core.T'), [300., 300., 300.])
    assert reg.differential_mask().tolist() == [True, True, True, False]
    with pytest.raises(ValueError):
        reg.register('core.T', kind='differential', initial=301., unit='K')
    reg.freeze()
    with pytest.raises(RuntimeError):
        reg.register('new', kind='algebraic', initial=0., unit='-')


def test_two_port_ideal_connection_constraint():
    reg = StateRegistry()
    vals = [('a.p', 5000.), ('a.h', 20.), ('a.m', 2.),
            ('b.p', 5000.), ('b.h', 20.), ('b.m', -2.)]
    for name, value in vals:
        reg.register(name, kind='algebraic', initial=value, unit='-')
    a = FluidPort('a.p', 'a.h', 'a.m')
    b = FluidPort('b.p', 'b.h', 'b.m')
    assert np.allclose(ideal_fluid_two_port_residual(reg.initial_state(), reg, a, b), 0.)


def test_thermal_connection_closes_global_residual():
    class Boundary(Component):
        def __init__(self, name, prescribed):
            self.port = ThermalPort(f"{name}.T", f"{name}.Q")
            self.prescribed = prescribed

        def register(self, registry):
            registry.register(self.port.T, kind="algebraic", unit="K", initial=310.0)
            registry.register(self.port.Q_dot, kind="algebraic", unit="W", initial=0.0)

        def residual(self, t, y, ydot, registry):
            name, value = self.prescribed
            return np.array([registry.get(y, getattr(self.port, name)) - value])

    hot = Boundary("hot", ("T", 330.0))
    cold = Boundary("cold", ("Q_dot", 20.0))
    connection = partial(ideal_thermal_two_port_residual, a=hot.port, b=cold.port)
    system = SystemAssembly([hot, cold], connections=[connection])
    y, yp = solve_steady_state(system.residual, 0.0, system.registry.initial_state())
    assert np.allclose(y, [330.0, -20.0, 330.0, 20.0])
    assert np.array_equal(yp, np.zeros_like(y))
    result = DASSL().solve(system.residual, t_span=(0.0, 1.0), y0=y, ydot0=yp)
    assert np.allclose(result.y, y)


def test_assembly_rejects_nonvector_residual():
    class InvalidComponent(Component):
        def register(self, registry):
            registry.register("invalid.y", kind="algebraic", unit="1", initial=0.0)

        def residual(self, t, y, ydot, context):
            return np.array([[0.0]])

    system = SystemAssembly([InvalidComponent()])
    with pytest.raises(ValueError, match="1-D vector"):
        system.residual(0.0, system.registry.initial_state(), np.zeros(1))


def test_stream_junction_mixes_two_inflows_and_handles_reversal_and_stagnation():
    registry = StateRegistry()
    ports = tuple(
        FluidStreamPort(f"{name}.p", f"{name}.m", f"{name}.h_out", f"{name}.h_in")
        for name in ("a", "b", "c")
    )
    for port in ports:
        for name, unit in (
            (port.p, "Pa"), (port.m_dot, "kg/s"),
            (port.h_outflow, "J/kg"), (port.h_instream, "J/kg"),
        ):
            registry.register(name, kind="algebraic", initial=0.0, unit=unit)

    def state(flows, outflows, instreams):
        values = registry.initial_state()
        for port, flow, outflow, instream in zip(ports, flows, outflows, instreams):
            for name, value in (
                (port.p, 5e6), (port.m_dot, flow),
                (port.h_outflow, outflow), (port.h_instream, instream),
            ):
                values[registry.spec(name).indices] = value
        return values

    def connection(flows, instreams):
        return ideal_fluid_stream_junction_residual(
            state(flows, [400.0, 300.0, 350.0], instreams),
            registry, ports, mass_flow_epsilon=1e-8,
        )

    forward = connection([-2.0, -1.0, 3.0], [300.0, 400.0, 1100.0 / 3.0])
    assert np.allclose(forward, 0.0, atol=1e-10)

    reverse = connection([2.0, 1.0, -3.0], [350.0, 350.0, 0.0])
    assert np.allclose(reverse[:4], 0.0, atol=1e-10)
    assert np.isfinite(reverse).all()

    stagnant = connection([0.0, 0.0, 0.0], [325.0, 375.0, 350.0])
    assert np.allclose(stagnant, 0.0, atol=1e-10)
    with pytest.raises(ValueError, match="mass_flow_epsilon"):
        ideal_fluid_stream_junction_residual(
            registry.initial_state(), registry, ports, mass_flow_epsilon=0.0,
        )


def test_stream_junction_closes_two_component_steady_state():
    class Endpoint(Component):
        def __init__(self, name, prescribed):
            self.port = FluidStreamPort(
                f"{name}.p", f"{name}.m", f"{name}.h_out", f"{name}.h_in",
            )
            self.prescribed = prescribed

        def register(self, registry):
            for name, unit in (
                (self.port.p, "Pa"), (self.port.m_dot, "kg/s"),
                (self.port.h_outflow, "J/kg"), (self.port.h_instream, "J/kg"),
            ):
                registry.register(name, kind="algebraic", initial=1.0, unit=unit)

        def residual(self, t, y, ydot, registry):
            return np.array([
                registry.get(y, getattr(self.port, name)) - value
                for name, value in self.prescribed.items()
            ])

    source = Endpoint("source", {"p": 5e6, "m_dot": -2.0, "h_outflow": 400.0})
    sink = Endpoint("sink", {"h_outflow": 300.0})
    connection = partial(
        ideal_fluid_stream_junction_residual,
        ports=(source.port, sink.port),
        mass_flow_epsilon=1e-6,
    )
    system = SystemAssembly([source, sink], connections=[connection])
    y, ydot = solve_steady_state(
        system.residual, 0.0, system.registry.initial_state(),
    )
    assert np.array_equal(ydot, np.zeros_like(y))
    assert np.max(np.abs(system.residual(0.0, y, ydot))) < 1e-8
    assert system.registry.get(y, source.port.h_instream) == pytest.approx(300.0)
    assert system.registry.get(y, sink.port.h_instream) == pytest.approx(400.0)
    result = DASSL().solve(system.residual, t_span=(0.0, 1.0), y0=y, ydot0=ydot)
    assert np.allclose(result.y, y)
