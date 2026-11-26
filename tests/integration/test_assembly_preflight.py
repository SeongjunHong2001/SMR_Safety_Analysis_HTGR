from functools import partial

import numpy as np
import pytest
from scipy import sparse

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from system.assembly import SystemAssembly
from system.connections import StreamJunction, ideal_fluid_stream_junction_residual


class Endpoint(Component):
    def __init__(self, name):
        self.port = FluidStreamPort(f"{name}.p", f"{name}.m", f"{name}.h_out",
                                    f"{name}.h_in")

    def register(self, registry):
        for key, unit in ((self.port.p, "Pa"), (self.port.m_dot, "kg/s"),
                          (self.port.h_outflow, "J/kg"),
                          (self.port.h_instream, "J/kg")):
            registry.register(key, kind="algebraic", initial=1.0, unit=unit)

    def residual(self, t, y, ydot, context):
        return np.empty(0)

    def ports(self):
        return {"fluid": self.port}


def test_same_port_cannot_enter_two_named_junctions():
    a, b = Endpoint("a"), Endpoint("b")
    junctions = [StreamJunction(name, (a.port, b.port), 1e-6, 1, 1, 1)
                 for name in ("first", "second")]
    with pytest.raises(ValueError, match=r"Port .*connection\[0\].*connection\[1\]"):
        SystemAssembly([a, b], connections=junctions)


def test_port_use_guard_includes_legacy_partial_connections():
    a, b = Endpoint("a"), Endpoint("b")
    junction = StreamJunction("first", (a.port, b.port), 1e-6, 1, 1, 1)
    legacy = partial(ideal_fluid_stream_junction_residual,
                     ports=(a.port, b.port), mass_flow_epsilon=1e-6)
    with pytest.raises(ValueError, match="used by connection"):
        SystemAssembly([a, b], connections=[junction, legacy])


def test_opaque_connection_cannot_bypass_port_ownership_check():
    with pytest.raises(TypeError, match="declare ports"):
        SystemAssembly([], connections=[lambda y, registry: np.empty(0)])


def test_preflight_rejects_square_but_rank_deficient_system():
    class Dependent(Component):
        def register(self, registry):
            for name in ("x", "y"):
                registry.register(name, kind="algebraic", initial=0.0, unit="1")

        def residual(self, t, y, ydot, context):
            return np.array([y[0] + y[1], 2 * (y[0] + y[1])])

    system = SystemAssembly([Dependent()])
    y = system.registry.initial_state()
    assert system.residual(0.0, y, np.zeros_like(y)).size == y.size
    with pytest.raises(ValueError, match="rank 1/2"):
        system.validate_for_solve(0.0, y, np.zeros_like(y), cj=1.0)
    with pytest.raises(ValueError, match="numerically singular"):
        system.validate_for_solve(0.0, y, np.zeros_like(y), cj=1.0,
                                  jacobian_sparsity=np.ones((2, 2), dtype=bool))


def test_sparse_preflight_accepts_full_rank_system():
    class Independent(Component):
        def register(self, registry):
            for name in ("x", "y"):
                registry.register(name, kind="algebraic", initial=0.0, unit="1")

        def residual(self, t, y, ydot, context):
            return y - np.array([1.0, 2.0])

    system = SystemAssembly([Independent()])
    y = system.registry.initial_state()
    diagnostic = system.validate_for_solve(
        0.0, y, np.zeros_like(y), cj=1.0,
        jacobian_sparsity=sparse.eye(2, format="csc"),
    )
    assert diagnostic.rank == diagnostic.size == 2
    assert diagnostic.condition is None
