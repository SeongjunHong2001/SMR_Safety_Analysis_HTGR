"""Synthetic two-cell closed loop with a single thermodynamic pressure.

All numbers and the linear density law are test values, not helium properties.
"""
import numpy as np

from interfaces.component import Component
from system.assembly import SystemAssembly
from system.inventory import total_helium_inventory


class SyntheticLowMachLoop(Component):
    volumes = np.array([1.0, 1.0])
    initial_inventory = 22.2

    def __init__(self, *, steady=False, closure="inventory", duplicate_mass=False):
        if closure not in ("inventory", "pressure"):
            raise ValueError("closure must be inventory or pressure")
        self.steady = steady
        self.closure = closure
        self.duplicate_mass = duplicate_mass

    def register(self, registry):
        registry.register("p_sys", kind="differential", unit="Pa", initial=10.0,
                          scale=10.0, owner="primary_loop")
        registry.register("loop.h_he", kind="differential", unit="J/kg",
                          shape=(2,), initial=[12.0, 10.0], scale=10.0)
        registry.register("loop.p_fluid_cell", kind="algebraic", unit="Pa",
                          shape=(2,), initial=[12.0, 10.0], scale=10.0)
        registry.register("loop.m_dot_fluid_face", kind="algebraic", unit="kg/s",
                          shape=(2,), initial=[2.0, 2.0], scale=2.0)

    def helium_volumes(self):
        return self.volumes.copy()

    def density(self, y, context):
        # rho(p_sys,h): the local momentum pressures are intentionally unused.
        return context.get(y, "p_sys") + 0.1 * context.get(y, "loop.h_he")

    def helium_inventory(self, t, y, context):
        return self.density(y, context) * self.volumes

    def internal_energy(self, y, context):
        rho = self.density(y, context)
        h = context.get(y, "loop.h_he")
        p_sys = context.get(y, "p_sys")
        return float(np.sum(self.volumes * (rho * h - p_sys)))

    def residual(self, t, y, ydot, context):
        p_sys = context.get(y, "p_sys")
        h = context.get(y, "loop.h_he")
        p_local = context.get(y, "loop.p_fluid_cell")
        m = context.get(y, "loop.m_dot_fluid_face")
        pdot = context.get(ydot, "p_sys")
        hdot = context.get(ydot, "loop.h_he")
        rho = self.density(y, context)
        rho_dot = pdot + 0.1 * hdot
        mass = self.volumes * rho_dot + np.array([m[0] - m[1], m[1] - m[0]])
        q_hot = 4.0 + 0.5 * np.sin(t)
        q_sink = -4.0 - 0.5 * (h[1] - 10.0)
        energy = rho * self.volumes * hdot - self.volumes * pdot + np.array([
            m[1] * (h[0] - h[1]) - q_hot,
            m[0] * (h[1] - h[0]) - q_sink,
        ])
        if self.steady:
            mass[1] = (total_helium_inventory([self], t, y, context)
                       - self.initial_inventory if self.closure == "inventory"
                       else p_sys - 10.0)
        result = np.array([mass[0], mass[1], *energy,
                           p_local[0] - p_local[1] - m[0],
                           p_local[1] - p_local[0] + 4.0 - m[1],
                           p_local[1] - p_sys])
        if self.duplicate_mass:
            result = np.append(result, mass.sum())
        return result

    def residual_ids(self):
        second = ("loop.initial_inventory" if self.closure == "inventory"
                  else "loop.initial_pressure") if self.steady else "loop.mass[1]"
        ids = ("loop.mass[0]", second, "loop.energy[0]", "loop.energy[1]",
               "loop.momentum[0]", "loop.momentum[1]", "loop.pressure_reference")
        return ids + (("loop.duplicate_global_mass",) if self.duplicate_mass else ())

    def residual_scales(self):
        return np.array([1.0, 20.0 if self.steady else 1.0,
                         10.0, 10.0, 10.0, 10.0, 10.0]
                        + ([1.0] if self.duplicate_mass else []))

    def outputs(self, t, y, context):
        return {"loop.M_he_loop": total_helium_inventory([self], t, y, context),
                "loop.U_he_loop": self.internal_energy(y, context)}


def make_low_mach_loop(**kwargs):
    return SystemAssembly([SyntheticLowMachLoop(**kwargs)])
