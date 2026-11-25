"""Role 4 SG wall contract; implementation awaits radial and axial mesh."""
from typing import Protocol


class SGWall(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D wall residuals with net current-trial W heat input."""
        ...
