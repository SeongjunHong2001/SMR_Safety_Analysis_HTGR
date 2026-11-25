"""Three-capacity SG fixture; every number is synthetic test data."""
import numpy as np

from interfaces.component import Component
from interfaces.ports import FluidStreamPort
from models.heat_exchanger import HeatExchanger
from system.assembly import SystemAssembly


class SyntheticSide:
    def __init__(self, label, medium, capacity, temperature):
        self.label = label
        self.medium = medium
        self.capacity = capacity
        self.temperature = temperature
        self.port_names = (f"{label}_inlet", f"{label}_outlet")

    def _key(self, suffix):
        return f"sg.{self.label}.{suffix}"

    def _port_fields(self):
        for port in self.port_names:
            for field, unit, target in (("p", "Pa", 10.0),
                                        ("m", "kg/s", 2.0),
                                        ("h_out", "J/kg", 100.0),
                                        ("h_in", "J/kg", 100.0)):
                yield self._key(f"{port}.{field}"), unit, target

    def register(self, registry):
        registry.register(self._key("T"), kind="differential", unit="K",
                          initial=self.temperature, scale=300.0)
        for key, unit, initial in self._port_fields():
            registry.register(key, kind="algebraic", unit=unit,
                              initial=initial, scale=max(initial, 1.0))

    def residual(self, t, y, ydot, context, *, heat_rate):
        energy = self.capacity * context.get(ydot, self._key("T")) - heat_rate[0]
        return np.array([energy, *self.port_residuals(y, context)])

    def port_residuals(self, y, context):
        return [context.get(y, key) - target
                for key, _, target in self._port_fields()]

    def ports(self):
        return {
            port: FluidStreamPort(
                self._key(f"{port}.p"), self._key(f"{port}.m"),
                self._key(f"{port}.h_out"), self._key(f"{port}.h_in"),
                medium=self.medium,
            ) for port in self.port_names
        }

    def residual_ids(self):
        return (self._key("energy"),
                *(f"{key}.fixture_constraint" for key, _, _ in self._port_fields()))

    def residual_scales(self):
        return np.array([100.0, *[max(value, 1.0) for _, _, value
                                 in self._port_fields()]])

    def outputs(self, t, y, context):
        return {self._key("temperature"): context.get(y, self._key("T"))}


class SyntheticWall:
    capacity = 3.0

    def register(self, registry):
        registry.register("sg.wall.T", kind="differential", unit="K",
                          initial=320.0, scale=300.0)

    def residual(self, t, y, ydot, context, *, heat_rate):
        return np.array([self.capacity * context.get(ydot, "sg.wall.T")
                         - heat_rate[0]])

    def residual_ids(self):
        return ("sg.wall.energy",)

    def residual_scales(self):
        return np.array([100.0])

    def outputs(self, t, y, context):
        return {"sg.wall.temperature": context.get(y, "sg.wall.T")}


class SyntheticSGCoupling:
    def evaluate(self, t, y, context):
        he = context.get(y, "sg.helium.T")
        wall = context.get(y, "sg.wall.T")
        secondary = context.get(y, "sg.secondary.T")
        q_he_wall = 2.0 * (he - wall)
        q_wall_secondary = 1.5 * (wall - secondary)
        return (np.array([-q_he_wall]), np.array([q_he_wall]),
                np.array([-q_wall_secondary]), np.array([q_wall_secondary]))


def make_parts():
    return (SyntheticSide("helium", "helium", 2.0, 350.0), SyntheticWall(),
            SyntheticSide("secondary", "water", 4.0, 300.0))


class FlatSyntheticSG(Component):
    """Independent explicit energy equations in the composite's state order."""

    def __init__(self):
        self.helium, self.wall, self.secondary = make_parts()

    def register(self, registry):
        for part in (self.helium, self.wall, self.secondary):
            part.register(registry)

    def residual(self, t, y, ydot, context):
        he = context.get(y, "sg.helium.T")
        wall = context.get(y, "sg.wall.T")
        secondary = context.get(y, "sg.secondary.T")
        q_he_wall = 2.0 * (he - wall)
        q_wall_secondary = 1.5 * (wall - secondary)
        return np.array([
            2.0 * context.get(ydot, "sg.helium.T") + q_he_wall,
            *self.helium.port_residuals(y, context),
            3.0 * context.get(ydot, "sg.wall.T") - q_he_wall + q_wall_secondary,
            4.0 * context.get(ydot, "sg.secondary.T") - q_wall_secondary,
            *self.secondary.port_residuals(y, context),
        ])

    def ports(self):
        return {**self.helium.ports(), **self.secondary.ports()}

    def residual_ids(self):
        return (self.helium.residual_ids() + self.wall.residual_ids()
                + self.secondary.residual_ids())

    def residual_scales(self):
        return np.concatenate((self.helium.residual_scales(),
                               self.wall.residual_scales(),
                               self.secondary.residual_scales()))


def make_composite():
    helium, wall, secondary = make_parts()
    model = HeatExchanger(helium_side=helium, wall=wall,
                          secondary_side=secondary,
                          coupling=SyntheticSGCoupling(), name="sg")
    return SystemAssembly([model])


def make_flat():
    return SystemAssembly([FlatSyntheticSG()])
