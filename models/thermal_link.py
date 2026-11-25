"""Derived heat exchange between two separate temperature nodes."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HeatExchange:
    to_a: np.ndarray
    to_b: np.ndarray


class ThermalLink:
    """Evaluate one injected transfer path once per trial, then share both signs."""

    def __init__(self, name, transfer):
        if not name or not callable(transfer):
            raise ValueError("ThermalLink needs a name and callable transfer")
        self.name = name
        self.transfer = transfer

    def evaluate(self, t, y, context):
        to_a, to_b = self.transfer(t, y, context)
        a = np.asarray(to_a, dtype=float).copy()
        b = np.asarray(to_b, dtype=float).copy()
        if a.ndim != 1 or b.ndim != 1 or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError(f"{self.name}: heat rates must be finite 1-D W arrays")
        if not np.isclose(a.sum() + b.sum(), 0.0, rtol=1e-10, atol=1e-9):
            raise ValueError(f"{self.name}: heat transfer is not conservative")
        a.setflags(write=False)
        b.setflags(write=False)
        return HeatExchange(a, b)
