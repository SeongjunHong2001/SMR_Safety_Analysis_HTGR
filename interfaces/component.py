"""전역 DAE에 참여하는 Component의 공통 계약."""
from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping
import numpy as np

class Component(ABC):
    """시간전진이 없는 순수 Residual 제공자.

    register()는 모듈이 소유한 **독립 변수만** 등록한다.
    다른 모듈 변수가 필요하면 전역 Registry에서 읽는다.
    """

    @property
    def name(self) -> str:
        return getattr(self, "_name", type(self).__name__.lower())

    @name.setter
    def name(self, value: str) -> None:
        if not isinstance(value, str) or not value:
            raise ValueError("Component name must be a nonempty string")
        self._name = value

    def key(self, variable: str) -> str:
        if not isinstance(variable, str) or not variable:
            raise ValueError("Variable name must be a nonempty string")
        return f"{self.name}.{variable}"

    @abstractmethod
    def register(self, registry):
        """Register owned unknown variables; no time integration."""

    @abstractmethod
    def residual(self, t: float, y: np.ndarray, ydot: np.ndarray, context) -> np.ndarray:
        """Return 1-D residual vector for the current global Newton guess."""

    def outputs(self, t: float, y: np.ndarray, context) -> dict:
        """Return postprocessing values; not used by the time integrator."""
        return {}

    def ports(self) -> Mapping[str, object]:
        """Named external ports. Resolve their Registry keys after registration."""
        return {}

    def residual_ids(self) -> tuple[str, ...] | None:
        """Fixed scalar order; None retains the legacy generated-ID behavior."""
        return None

    def residual_scales(self) -> np.ndarray | None:
        """Positive physical scale for each residual; distinct from state scale."""
        return None
