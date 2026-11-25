"""Role 4 SG-specific correlation contract; no unvalidated default."""
from typing import Protocol
from interfaces.contracts import HeatTransferCorrelation


class SGHeliumCrossflow(HeatTransferCorrelation, Protocol):
    """coefficient(state, geometry) for validated SG helium crossflow."""
