"""Synthetic closed-loop pressure/flow ledger; not an HTGR flow model."""
import numpy as np

from interfaces.component import Component
from system.assembly import SystemAssembly


class SyntheticLoop(Component):
    def __init__(self, *, steady_inventory=False, duplicate_mass=False):
        self.steady_inventory = steady_inventory
        self.duplicate_mass = duplicate_mass

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=10.0,
                          scale=10.0, owner="system")
        for name in ("p_a", "p_b"):
            registry.register(f"loop.{name}", kind="algebraic", unit="Pa", initial=10.0,
                              scale=10.0)
        for name in ("m_ab", "m_ba"):
            registry.register(f"loop.{name}", kind="algebraic", unit="kg/s", initial=2.0,
                              scale=2.0)

    def residual(self, t, y, ydot, context):
        p_sys, p_a, p_b, m_ab, m_ba = y
        mass_or_inventory = p_sys - 10.0 if self.steady_inventory else ydot[0]
        values = [p_a - p_b - m_ab,
                  p_b - p_a + 6.0 - 2.0 * m_ba,
                  m_ab - m_ba,
                  p_a - p_sys,
                  mass_or_inventory]
        if self.duplicate_mass:
            values.append(ydot[0])
        return np.array(values)

    def residual_ids(self):
        ids = ("loop.momentum_ab", "loop.momentum_ba", "loop.cell_mass",
               "loop.pressure_reference",
               "loop.initial_inventory" if self.steady_inventory else "loop.mass_dynamics")
        return ids + (("loop.duplicate_global_mass",) if self.duplicate_mass else ())

    def residual_scales(self):
        values = [10.0, 10.0, 2.0, 10.0, 10.0]
        return np.array(values + ([10.0] if self.duplicate_mass else []))


def make_loop(**kwargs):
    return SystemAssembly([SyntheticLoop(**kwargs)])
