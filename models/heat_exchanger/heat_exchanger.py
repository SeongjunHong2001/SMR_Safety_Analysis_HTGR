"""SG composite boundary: injected submodels share one Newton trial.

The injected coupling owns the physical heat-transfer calculation. Production
couplings must use ``physics.heat_transfer``; this class only maps the resulting
W arrays to the three energy residuals.
"""
from __future__ import annotations

import numpy as np

from interfaces.component import Component


class HeatExchanger(Component):
    """Assemble helium side, wall and secondary side without registering twice.

    ``coupling.evaluate(t, y, context)`` returns four heat-received arrays in W:
    (to_helium, to_wall_from_helium, to_wall_from_secondary, to_secondary).
    Each interface sums to zero, and both wall arrays use the wall grid.
    """

    def __init__(self, config=None, properties=None, *, helium_side=None,
                 wall=None, secondary_side=None, coupling=None,
                 name="heat_exchanger"):
        self.name = name
        self.config = {} if config is None else config
        self.properties = properties
        self.helium_side = helium_side
        self.wall = wall
        self.secondary_side = secondary_side
        self.coupling = coupling

    def _parts(self):
        missing = [path for path, part in (
            ("models/heat_exchanger/helium_side.py", self.helium_side),
            ("models/heat_exchanger/wall.py", self.wall),
            ("models/heat_exchanger/secondary_side.py", self.secondary_side),
            ("models/heat_exchanger/heat_exchanger.py:coupling", self.coupling),
        ) if part is None]
        if missing:
            raise NotImplementedError("HeatExchanger requires injected implementations: "
                                      + ", ".join(missing))
        parts = (self.helium_side, self.wall, self.secondary_side)
        if len({id(part) for part in parts}) != len(parts):
            raise ValueError("HeatExchanger submodels must be distinct objects")
        return parts

    def register(self, registry):
        for part in self._parts():
            part.register(registry)

    def residual(self, t, y, ydot, context):
        parts = self._parts()
        rates = self.coupling.evaluate(t, y, context)
        if not isinstance(rates, (tuple, list)) or len(rates) != 4:
            raise ValueError("HeatExchanger coupling must return four W arrays")
        q_he, q_wall_he, q_wall_secondary, q_secondary = (
            np.asarray(rate, dtype=float) for rate in rates
        )
        for label, rate in zip(("helium", "wall from helium", "wall from secondary",
                                "secondary"),
                               (q_he, q_wall_he, q_wall_secondary, q_secondary)):
            if rate.ndim != 1 or not rate.size or not np.isfinite(rate).all():
                raise ValueError(f"HeatExchanger {label} heat rate must be a nonempty finite 1-D W array")
        if q_wall_he.shape != q_wall_secondary.shape:
            raise ValueError("HeatExchanger wall heat-rate grids differ")
        if not np.isclose(q_he.sum(), -q_wall_he.sum(), rtol=1e-10, atol=1e-9):
            raise ValueError("HeatExchanger helium/wall heat rate is not conservative")
        if not np.isclose(q_wall_secondary.sum(), -q_secondary.sum(),
                          rtol=1e-10, atol=1e-9):
            raise ValueError("HeatExchanger wall/secondary heat rate is not conservative")

        heat_rates = (q_he, q_wall_he + q_wall_secondary, q_secondary)
        residuals = []
        for part, heat_rate in zip(parts, heat_rates):
            value = np.asarray(part.residual(t, y, ydot, context,
                                             heat_rate=heat_rate), dtype=float)
            if value.ndim != 1 or not np.isfinite(value).all():
                raise ValueError("HeatExchanger submodel residual must be a finite 1-D array")
            residuals.append(value)
        return np.concatenate(residuals)

    def ports(self):
        self._parts()
        helium = dict(self.helium_side.ports())
        secondary = dict(self.secondary_side.ports())
        if not {"helium_inlet", "helium_outlet"} <= helium.keys():
            raise ValueError("HeatExchanger helium side must expose inlet and outlet")
        secondary_ports = {"secondary_inlet", "secondary_outlet"} & secondary.keys()
        if len(secondary_ports) == 1:
            raise ValueError("HeatExchanger secondary inlet/outlet must appear together")
        if helium.keys() & secondary.keys():
            raise ValueError("HeatExchanger submodel port names overlap")
        return {**helium, **secondary}

    def residual_ids(self):
        return self._metadata("residual_ids")

    def residual_scales(self):
        return np.concatenate(self._metadata("residual_scales"))

    def _metadata(self, method):
        result = []
        for part in self._parts():
            getter = getattr(part, method, None)
            if not callable(getter):
                raise TypeError(f"HeatExchanger submodels must define {method}()")
            value = getter()
            if value is None:
                raise TypeError(f"HeatExchanger submodels must define {method}()")
            result.append(value)
        if method == "residual_ids":
            return tuple(item for group in result for item in group)
        return result

    def outputs(self, t, y, context):
        values = {}
        for part in self._parts():
            getter = getattr(part, "outputs", None)
            if not callable(getter):
                continue
            part_values = getter(t, y, context)
            if values.keys() & part_values.keys():
                raise ValueError("HeatExchanger submodel output names overlap")
            values.update(part_values)
        return values
