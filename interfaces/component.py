"""담당 팀원에게 제공할 물리 Component의 최소 인터페이스."""
from abc import ABC, abstractmethod
import numpy as np

class Component(ABC):
    """시간전진이 없는 순수 Residual 제공자.

    register()는 모듈이 소유한 **독립 변수만** 등록한다.
    다른 모듈 변수가 필요하면 전역 Registry에서 읽는다.
    """
    @abstractmethod
    def register(self, registry):
        """Register owned unknown variables; no time integration."""

    @abstractmethod
    def residual(self, t: float, y: np.ndarray, ydot: np.ndarray, context) -> np.ndarray:
        """Return 1-D residual vector for the current global Newton guess."""

    def outputs(self, t: float, y: np.ndarray, context) -> dict:
        """Return postprocessing values; not used by the time integrator."""
        return {}
