"""Role 3 Core coupling contract; geometry mapping will live here."""
from typing import Protocol


class CoreCoupling(Protocol):
    def evaluate(self, t, y, context):
        """Return (Q_to_fluid, Q_from_solid) in W with equal sums."""
        ...
