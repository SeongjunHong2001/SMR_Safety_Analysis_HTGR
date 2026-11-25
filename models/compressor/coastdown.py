"""Role 2 compressor coastdown provider contract."""
from typing import Protocol


class CompressorCoastdown(Protocol):
    def evaluate(self, t, speed, scenario):
        """Return externally defined shaft torque at the current trial."""
        ...
