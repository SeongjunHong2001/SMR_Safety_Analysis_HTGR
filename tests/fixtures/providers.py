"""Deliberately synthetic property/correlation/boundary providers for tests."""
import numpy as np


class SyntheticHeliumProperties:
    """Ideal-gas arithmetic for tests; these constants are not reactor data."""

    R = 2.0
    CP = 5.0

    def temperature(self, p, h):
        return np.asarray(h) / self.CP

    def enthalpy(self, p, T):
        return self.CP * np.asarray(T)

    def density(self, p, h):
        return np.asarray(p) * self.CP / (self.R * np.asarray(h))

    def viscosity(self, p, T):
        return np.broadcast_arrays(p, T)[0] * 0 + 1.0

    def conductivity(self, p, T):
        return np.broadcast_arrays(p, T)[0] * 0 + 1.0

    def cp(self, p, T):
        return np.broadcast_arrays(p, T)[0] * 0 + self.CP

    def cv(self, p, T):
        return np.broadcast_arrays(p, T)[0] * 0 + (self.CP - self.R)

    def drho_dp(self, p, h):
        return self.CP / (self.R * np.asarray(h)) + np.asarray(p) * 0

    def drho_dh(self, p, h):
        return -np.asarray(p) * self.CP / (self.R * np.asarray(h)**2)


class SyntheticHeatCorrelation:
    def coefficient(self, state, geometry):
        return np.broadcast_to(2.0, np.broadcast_shapes(np.shape(state["T_wall"]),
                                                     np.shape(geometry["A"])))


class SyntheticSGBoundary:
    def heat_sink(self, t, signal):
        return 10.0 * signal
