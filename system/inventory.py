"""Shared closed-loop inventory sum; local mass equations stay with components."""
import numpy as np


def total_helium_inventory(components, t, y, context) -> float:
    """Sum each selected helium volume exactly once, returning kg.

    The builder supplies inventory-bearing components only. Their density is
    evaluated at the single thermodynamic ``p_sys`` by the component's EOS.
    This derived value is never an extra transient mass residual.
    """
    components = tuple(components)
    if not components:
        raise ValueError("Helium inventory needs at least one component")
    if len({id(component) for component in components}) != len(components):
        raise ValueError("Duplicate helium inventory component")
    total = 0.0
    for component in components:
        volumes = np.asarray(component.helium_volumes(), dtype=float)
        masses = np.asarray(component.helium_inventory(t, y, context), dtype=float)
        if (volumes.ndim != 1 or not volumes.size or masses.shape != volumes.shape
                or not np.isfinite(volumes).all() or not np.isfinite(masses).all()
                or np.any(volumes <= 0) or np.any(masses <= 0)):
            raise ValueError("Helium inventory volumes and masses must be positive finite matching 1-D arrays")
        total += float(masses.sum())
    return total
