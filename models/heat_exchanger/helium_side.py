"""Role 4 helium-side SG contract; implementation awaits geometry and EOS."""
from typing import Protocol


class HeliumSide(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D primary-side residuals with current-trial W heat input."""
        ...
    def ports(self): ...
