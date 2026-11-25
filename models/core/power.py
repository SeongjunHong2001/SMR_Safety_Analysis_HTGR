"""Role 3 Core power source contract; no default reactor power is assumed."""
from typing import Protocol


class CorePower(Protocol):
    def evaluate(self, t, y, context):
        """Return a finite 1-D W array on the solid grid for this trial."""
        ...


class StatefulCorePower(CorePower, Protocol):
    """Optional point-kinetics state owner, with its own residual ledger."""

    def register(self, registry): ...
    def residual(self, t, y, ydot, context): ...
    def residual_ids(self): ...
    def residual_scales(self): ...
