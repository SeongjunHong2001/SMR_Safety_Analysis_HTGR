from functools import partial

import numpy as np
import pytest
from interfaces.component import Component
from interfaces.state_registry import StateRegistry
from interfaces.ports import FluidPort, ThermalPort
from solvers import DASSL, solve_steady_state
from system.assembly import SystemAssembly
from system.connections import ideal_fluid_two_port_residual, ideal_thermal_two_port_residual


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
    system.context = system.registry
    y, yp = solve_steady_state(system.residual, 0.0, system.registry.initial_state())
    assert np.allclose(y, [330.0, -20.0, 330.0, 20.0])
    assert np.array_equal(yp, np.zeros_like(y))
    result = DASSL().solve(system.residual, t_span=(0.0, 1.0), y0=y, ydot0=yp)
    assert np.allclose(result.y, y)
