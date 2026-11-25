"""Provider and trial-context contracts; no physical correlations live here."""
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Protocol, runtime_checkable

import numpy as np

from interfaces.state_registry import StateRegistry


def _readonly(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _readonly(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_readonly(item) for item in value)
    if isinstance(value, np.ndarray):
        array = value.copy()
        array.setflags(write=False)
        return array
    return value


@dataclass(frozen=True)
class SimulationContext:
    """Configuration and providers chosen before integration, plus one trial's rates.

    The default legacy context remains the bare StateRegistry. ``heat_rates`` is
    replaced for each trial and must never hold accepted-step history.
    """

    registry: StateRegistry
    config: Mapping = field(default_factory=dict)
    scenario: Mapping = field(default_factory=dict)
    providers: Mapping = field(default_factory=dict)
    heat_rates: Mapping = field(default_factory=dict)

    def __post_init__(self):
        for name in ("config", "scenario", "providers", "heat_rates"):
            object.__setattr__(self, name, _readonly(getattr(self, name)))

    def with_heat_rates(self, rates: Mapping) -> "SimulationContext":
        return replace(self, heat_rates=rates)

    def get(self, vector, name):
        """Legacy-compatible shorthand for ``registry.get``."""
        return self.registry.get(vector, name)


@runtime_checkable
class HeliumProperties(Protocol):
    """SI units, scalar/NumPy broadcasting; valid domain is provider-specific."""

    def temperature(self, p, h): ...
    def enthalpy(self, p, T): ...
    def density(self, p, h): ...
    def viscosity(self, p, T): ...
    def conductivity(self, p, T): ...
    def cp(self, p, T): ...
    def cv(self, p, T): ...
    def drho_dp(self, p, h): ...
    def drho_dh(self, p, h): ...


@runtime_checkable
class HeatTransferCorrelation(Protocol):
    """Current trial state and geometry -> h_conv in W/(m² K)."""

    def coefficient(self, state: Mapping, geometry: Mapping) -> np.ndarray: ...


@runtime_checkable
class HeliumInventory(Protocol):
    """One owner per helium volume; arrays use the component's cell order."""

    def helium_volumes(self) -> np.ndarray:
        """Return a finite positive 1-D array of owned control-volume m³."""
        ...

    def helium_inventory(self, t, y, context) -> np.ndarray:
        """Return the corresponding rho(p_sys, h) * V masses in kg."""
        ...
