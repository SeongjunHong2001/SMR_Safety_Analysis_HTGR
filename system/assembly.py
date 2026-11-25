"""물리모델 순회는 Residual 평가 순서일 뿐 시간전진 순서가 아님."""
from dataclasses import dataclass
from functools import partial
import numpy as np
from interfaces.contracts import SimulationContext
from interfaces.state_registry import StateRegistry
from models.thermal_link import ThermalLink
from solvers.jacobian import finite_difference_jacobian
from system.connections import validate_port, validate_stream_junction


@dataclass(frozen=True)
class JacobianDiagnostic:
    rank: int
    size: int
    condition: float
    singular_values: np.ndarray
    residual_ids: tuple[str, ...]

class SystemAssembly:
    def __init__(self, components, *, context=None, connections=(), thermal_links=()):
        """Assemble components, fluid connections and derived thermal links."""
        self.components = tuple(components)
        self.connections = tuple(connections)
        self.thermal_links = tuple(thermal_links)
        if any(not isinstance(link, ThermalLink) for link in self.thermal_links):
            raise TypeError("thermal_links must contain ThermalLink instances")
        rate_keys = [key for link in self.thermal_links for key in link.rate_keys]
        if len(set(rate_keys)) != len(rate_keys):
            raise ValueError("Duplicate ThermalLink heat-rate keys")
        self.registry = StateRegistry()
        for component in self.components:
            component.register(self.registry)
        self.registry.freeze()
        if context is None:
            self.context = (SimulationContext(self.registry) if self.thermal_links
                            else self.registry)
        elif isinstance(context, SimulationContext):
            self.context = SimulationContext(
                self.registry, config=context.config, scenario=context.scenario,
                providers=context.providers, heat_rates=context.heat_rates,
            )
        else:
            self.context = context
        if self.thermal_links and not isinstance(self.context, SimulationContext):
            raise TypeError("ThermalLink assembly requires SimulationContext")
        if self.thermal_links and set(rate_keys) & set(self.context.heat_rates):
            raise ValueError("ThermalLink heat-rate key conflicts with context")
        for component in self.components:
            for port in component.ports().values():
                validate_port(port, self.registry)
        for connection in self.connections:
            if hasattr(connection, "validate"):
                connection.validate(self.registry)
            elif isinstance(connection, partial):
                arguments = connection.keywords or {}
                if "ports" in arguments:
                    validate_stream_junction(arguments["ports"], self.registry)
                for key in ("a", "b"):
                    if key in arguments:
                        validate_port(arguments[key], self.registry)

    def residual(self, t, y, ydot):
        return self.residual_details(t, y, ydot)[0]

    def residual_details(self, t, y, ydot):
        """Return residual, fixed scalar IDs, and physical residual scales."""
        y = np.asarray(y, dtype=float)
        ydot = np.asarray(ydot, dtype=float)
        n = self.registry.size
        if y.shape != (n,) or ydot.shape != (n,):
            raise ValueError(f"Expected y and ydot shape ({n},)")
        # Components must observe one immutable Newton trial. Views avoid copying
        # the global vectors on every residual/Jacobian evaluation.
        trial_y, trial_ydot = y.view(), ydot.view()
        trial_y.setflags(write=False)
        trial_ydot.setflags(write=False)
        trial_context = self._trial_context(t, trial_y)
        r_parts, ids, scales = [], [], []
        for index, component in enumerate(self.components):
            source = f"{component.name}[{index}]"
            part = self._residual_part(component.residual(t, trial_y, trial_ydot, trial_context), source)
            r_parts.append(part)
            part_ids, part_scales = self._metadata(component, source, part.size)
            ids.extend(part_ids)
            scales.extend(part_scales)
        for index, connection in enumerate(self.connections):
            source = f"connection[{index}]"
            part = self._residual_part(connection(trial_y, self.registry), source)
            r_parts.append(part)
            part_ids, part_scales = self._metadata(connection, source, part.size)
            ids.extend(part_ids)
            scales.extend(part_scales)
        r = np.concatenate(r_parts) if r_parts else np.empty(0)
        if r.size != n:
            raise ValueError(f"DAE is not square: {n} unknowns, {r.size} residuals")
        if len(set(ids)) != len(ids):
            raise ValueError("Duplicate residual IDs in global assembly")
        return r, tuple(ids), np.asarray(scales)

    def _trial_context(self, t, y):
        if not self.thermal_links:
            return self.context
        rates = dict(self.context.heat_rates)
        for link in self.thermal_links:
            exchange = link.evaluate(t, y, self.context)
            a_key, b_key = link.rate_keys
            rates[a_key], rates[b_key] = exchange.to_a, exchange.to_b
        return self.context.with_heat_rates(rates)

    @staticmethod
    def _metadata(source_object, source, length):
        ids_method = getattr(source_object, "residual_ids", None)
        scale_method = getattr(source_object, "residual_scales", None)
        declared_ids = ids_method() if callable(ids_method) else None
        declared_scales = scale_method() if callable(scale_method) else None
        ids = tuple(f"{source}.R[{i}]" for i in range(length)) if declared_ids is None else tuple(declared_ids)
        scales = np.ones(length) if declared_scales is None else np.asarray(declared_scales, dtype=float)
        if len(ids) != length or len(set(ids)) != len(ids) or any(not isinstance(i, str) or not i for i in ids):
            raise ValueError(f"{source}: residual IDs must be unique and match residual length")
        if scales.shape != (length,) or not np.all(np.isfinite(scales)) or np.any(scales <= 0):
            raise ValueError(f"{source}: residual scales must be positive and match residual length")
        return ids, scales

    def residual_ids(self, t=0.0, y=None, ydot=None):
        y = self.registry.initial_state() if y is None else y
        ydot = np.zeros(self.registry.size) if ydot is None else ydot
        return self.residual_details(t, y, ydot)[1]

    def residual_scales(self, t=0.0, y=None, ydot=None):
        y = self.registry.initial_state() if y is None else y
        ydot = np.zeros(self.registry.size) if ydot is None else ydot
        return self.residual_details(t, y, ydot)[2]

    def diagnose_jacobian(self, t, y, ydot, *, cj=0.0):
        """Rank/condition of diag(1/s_R)(F_y+cj F_ydot)diag(s_y)."""
        y = np.asarray(y, dtype=float)
        ydot = np.asarray(ydot, dtype=float)
        if not np.isfinite(cj):
            raise ValueError("cj must be finite")
        _, ids, residual_scales = self.residual_details(t, y, ydot)
        state_scales = self.registry.scale_vector()
        scaled_y = y / state_scales

        def scaled_residual(z):
            delta = (z - scaled_y) * state_scales
            return self.residual(t, y + delta, ydot + cj * delta) / residual_scales

        jacobian = finite_difference_jacobian(scaled_residual, scaled_y,
                                               step_scales=1.0)
        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        tolerance = max(jacobian.shape) * np.finfo(float).eps * singular_values[0]
        rank = int(np.count_nonzero(singular_values > tolerance))
        condition = (float(singular_values[0] / singular_values[-1])
                     if singular_values[-1] > 0 else float("inf"))
        return JacobianDiagnostic(rank, y.size, condition, singular_values, ids)

    @staticmethod
    def _residual_part(values, source):
        values = np.asarray(values, dtype=float)
        if values.ndim != 1:
            raise ValueError(f"{source} residual must be a 1-D vector")
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"Non-finite residual from {source}")
        return values

    def outputs(self, t, y):
        y = np.asarray(y, dtype=float)
        if y.shape != (self.registry.size,):
            raise ValueError(f"Expected y shape ({self.registry.size},)")
        trial_y = y.view()
        trial_y.setflags(write=False)
        trial_context = self._trial_context(t, trial_y)
        out = {}
        for comp in self.components:
            values = comp.outputs(t, trial_y, trial_context)
            if out.keys() & values.keys():
                raise ValueError("Duplicate output names in global assembly")
            out.update(values)
        return out
