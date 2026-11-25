"""Role 4 secondary-side SG contract; implementation awaits water properties."""
from typing import Protocol


class SecondarySide(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D secondary-side residuals with current-trial W heat input."""
        ...
    def ports(self): ...
