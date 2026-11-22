"""물리모델 순회는 Residual 평가 순서일 뿐 시간전진 순서가 아님."""
import numpy as np
from interfaces.state_registry import StateRegistry

class SystemAssembly:
    def __init__(self, components, *, context=None, connections=()):
        """Assemble component residuals and connection callbacks (y, registry)."""
        self.components = tuple(components)
        self.connections = tuple(connections)
        self.context = context
        self.registry = StateRegistry()
        for component in self.components:
            component.register(self.registry)
        self.registry.freeze()

    def residual(self, t, y, ydot):
        y = np.asarray(y, dtype=float)
        ydot = np.asarray(ydot, dtype=float)
        n = self.registry.size
        if y.shape != (n,) or ydot.shape != (n,):
            raise ValueError(f"Expected y and ydot shape ({n},)")
        r_parts = [
            self._residual_part(component.residual(t, y, ydot, self.context), type(component).__name__)
            for component in self.components
        ]
        for index, connection in enumerate(self.connections):
            r_parts.append(self._residual_part(connection(y, self.registry), f"connection {index}"))
        r = np.concatenate(r_parts) if r_parts else np.empty(0)
        if r.size != n:
            raise ValueError(f"DAE is not square: {n} unknowns, {r.size} residuals")
        return r

    @staticmethod
    def _residual_part(values, source):
        values = np.asarray(values, dtype=float).ravel()
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"Non-finite residual from {source}")
        return values

    def outputs(self, t, y):
        out = {}
        for comp in self.components:
            out.update(comp.outputs(t, y, self.context))
        return out
