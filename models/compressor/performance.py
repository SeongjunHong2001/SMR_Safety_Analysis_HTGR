"""Role 2 compressor performance provider contract."""
from typing import Protocol


class CompressorPerformance(Protocol):
    def evaluate(self, *, speed, mass_flow, inlet_state, geometry):
        """Return pressure rise and shaft work using a validated map."""
        ...
