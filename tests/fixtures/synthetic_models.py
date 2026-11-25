"""Two-temperature synthetic DAE for contract and solver integration tests."""
import numpy as np

from interfaces.component import Component
from models.core import Core
from system.assembly import SystemAssembly


class SyntheticFluid(Component):
    def register(self, registry):
        registry.register("test.cold_T", kind="differential", unit="K",
                          initial=300.0, scale=300.0)

    def residual(self, t, y, ydot, context, *, heat_rate):
        cold = context.get(y, "test.cold_T")
        rate = context.get(ydot, "test.cold_T")
        return np.array([40.0 * rate - heat_rate[0] + 3.0 * (cold - 300.0)])

    def residual_ids(self):
        return ("test.cold_energy",)

    def residual_scales(self):
        return np.array([50.0])


class SyntheticSolid(Component):
    def register(self, registry):
        registry.register("test.hot_T", kind="differential", unit="K",
                          initial=320.0, scale=300.0)

    def residual(self, t, y, ydot, context, *, heat_rate):
        return np.array([20.0 * context.get(ydot, "test.hot_T") - heat_rate[0]])

    def residual_ids(self):
        return ("test.hot_energy",)

    def residual_scales(self):
        return np.array([50.0])


class SyntheticPower:
    def evaluate(self, t, y, context):
        return np.array([50.0])


class SyntheticCoupling:
    def evaluate(self, t, y, context):
        q = context.get(y, "test.q_dot")
        return np.array([q]), np.array([q])


class SyntheticHeatConstraint(Component):
    def register(self, registry):
        registry.register("test.q_dot", kind="algebraic", unit="W",
                          initial=0.0, scale=50.0)

    def residual(self, t, y, ydot, context):
        q = context.get(y, "test.q_dot")
        hot = context.get(y, "test.hot_T")
        cold = context.get(y, "test.cold_T")
        return np.array([q - 5.0 * (hot - cold)])

    def residual_ids(self):
        return ("test.heat_link",)

    def residual_scales(self):
        return np.array([50.0])


class FlatSynthetic(Component):
    def register(self, registry):
        SyntheticFluid().register(registry)
        SyntheticSolid().register(registry)
        SyntheticHeatConstraint().register(registry)

    def residual(self, t, y, ydot, context):
        q = context.get(y, "test.q_dot")
        return np.concatenate((
            SyntheticFluid().residual(t, y, ydot, context, heat_rate=np.array([q])),
            SyntheticSolid().residual(t, y, ydot, context, heat_rate=np.array([50.0-q])),
            SyntheticHeatConstraint().residual(t, y, ydot, context),
        ))

    def residual_ids(self):
        return ("test.cold_energy", "test.hot_energy", "test.heat_link")

    def residual_scales(self):
        return np.array([50.0, 50.0, 50.0])


class DerivedSynthetic(Component):
    """Same temperatures with Q derived rather than registered as a state."""

    def register(self, registry):
        SyntheticFluid().register(registry)
        SyntheticSolid().register(registry)

    def residual(self, t, y, ydot, context):
        q = 5.0 * (context.get(y, "test.hot_T") - context.get(y, "test.cold_T"))
        return np.concatenate((
            SyntheticFluid().residual(t, y, ydot, context, heat_rate=np.array([q])),
            SyntheticSolid().residual(t, y, ydot, context, heat_rate=np.array([50.0-q])),
        ))

    def residual_ids(self):
        return ("test.cold_energy", "test.hot_energy")

    def residual_scales(self):
        return np.array([50.0, 50.0])


def make_composite():
    core = Core(fluid=SyntheticFluid(), solid=SyntheticSolid(),
                power=SyntheticPower(), coupling=SyntheticCoupling())
    return SystemAssembly([core, SyntheticHeatConstraint()])


def make_flat():
    return SystemAssembly([FlatSynthetic()])


def make_derived():
    return SystemAssembly([DerivedSynthetic()])
