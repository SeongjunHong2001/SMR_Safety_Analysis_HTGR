"""Core composite: one assembly-visible model, no embedded physical formula."""
from __future__ import annotations

import numpy as np

from interfaces.component import Component


class Core(Component):
    """Coordinate injected fluid, solid, power and conservative coupling.

    Legacy ``Core(config, properties)`` construction still works. Registering or
    evaluating it requires the actual submodels; none is silently substituted.
    """

    def __init__(self, config=None, properties=None, *, fluid=None, solid=None,
                 power=None, coupling=None, name="core"):
        self.name = name
        self.config = {} if config is None else config
        self.properties = properties
        self.fluid = fluid
        self.solid = solid
        self.power = power
        self.coupling = coupling

    def _require_parts(self):
        missing = [path for path, value in (
            ("models/core/fluid.py", self.fluid),
            ("models/core/solid.py", self.solid),
            ("models/core/power.py", self.power),
            ("models/core/coupling.py", self.coupling),
        ) if value is None]
        if missing:
            raise NotImplementedError("Core requires injected implementations: " + ", ".join(missing))

    def register(self, registry):
        self._require_parts()
        self.fluid.register(registry)
        self.solid.register(registry)
        has_register = callable(getattr(self.power, "register", None))
        has_residual = callable(getattr(self.power, "residual", None))
        if has_register != has_residual:
            raise TypeError("Stateful Core power requires both register and residual")
        if has_register:
            if not all(callable(getattr(self.power, name, None))
                       for name in ("residual_ids", "residual_scales")):
                raise TypeError("Stateful Core power requires residual IDs and scales")
            self.power.register(registry)

    def residual(self, t, y, ydot, context):
        self._require_parts()
        q_to_fluid, q_from_solid = self.coupling.evaluate(t, y, context)
        q_to_fluid = np.asarray(q_to_fluid, dtype=float)
        q_from_solid = np.asarray(q_from_solid, dtype=float)
        q_generated = np.asarray(self.power.evaluate(t, y, context), dtype=float)
        for label, rates in (("fluid", q_to_fluid), ("solid", q_from_solid),
                             ("power", q_generated)):
            if rates.ndim != 1 or not np.all(np.isfinite(rates)):
                raise ValueError(f"Core {label} heat rate must be a finite 1-D W array")
        if q_generated.shape != q_from_solid.shape:
            raise ValueError("Core power and solid heat-rate shapes differ")
        if not np.isclose(q_to_fluid.sum(), q_from_solid.sum(), rtol=1e-10, atol=1e-9):
            raise ValueError("Core fluid/solid heat-rate mapping does not conserve W")
        fluid_r = np.asarray(self.fluid.residual(t, y, ydot, context,
                                                 heat_rate=q_to_fluid), dtype=float)
        solid_r = np.asarray(self.solid.residual(t, y, ydot, context,
                                                 heat_rate=q_generated - q_from_solid), dtype=float)
        power_r = (np.asarray(self.power.residual(t, y, ydot, context), dtype=float)
                   if callable(getattr(self.power, "residual", None)) else np.empty(0))
        if fluid_r.ndim != 1 or solid_r.ndim != 1 or power_r.ndim != 1:
            raise ValueError("Core submodel residuals must be 1-D")
        return np.concatenate((fluid_r, solid_r, power_r))

    def ports(self):
        self._require_parts()
        return {**{f"fluid.{k}": v for k, v in self.fluid.ports().items()},
                **{f"solid.{k}": v for k, v in self.solid.ports().items()}}

    def residual_ids(self):
        if self.fluid is None or self.solid is None:
            return None
        providers = [self.fluid, self.solid]
        if callable(getattr(self.power, "residual", None)):
            providers.append(self.power)
        ids = [part.residual_ids() for part in providers]
        return None if any(item is None for item in ids) else tuple(value for item in ids for value in item)

    def residual_scales(self):
        if self.fluid is None or self.solid is None:
            return None
        providers = [self.fluid, self.solid]
        if callable(getattr(self.power, "residual", None)):
            providers.append(self.power)
        scales = [part.residual_scales() for part in providers]
        return None if any(item is None for item in scales) else np.concatenate(scales)

    def outputs(self, t, y, context):
        self._require_parts()
        power_outputs = (self.power.outputs(t, y, context)
                         if callable(getattr(self.power, "outputs", None)) else {})
        return {**self.fluid.outputs(t, y, context),
                **self.solid.outputs(t, y, context), **power_outputs}
