"""전역 상태벡터 y/ydot의 이름, 소유권, 순서 관리."""
from dataclasses import dataclass
from typing import Literal
import numpy as np

Kind = Literal["differential", "algebraic"]

@dataclass(frozen=True)
class VariableSpec:
    name: str
    indices: slice
    shape: tuple[int, ...]
    kind: Kind
    unit: str
    initial: np.ndarray
    scale: float

    @property
    def size(self) -> int:
        return self.indices.stop - self.indices.start

class StateRegistry:
    """등록 이후 freeze()로 잠가 새 변수 등록을 금지."""
    def __init__(self):
        self._variables: dict[str, VariableSpec] = {}
        self._size = 0
        self._frozen = False

    def register(self, name: str, *, kind: Kind, initial, unit: str,
                 shape: tuple[int, ...] = (), scale: float = 1.0) -> VariableSpec:
        if self._frozen:
            raise RuntimeError("Registry is frozen")
        if name in self._variables:
            raise ValueError(f"Duplicate state name: {name}")
        if kind not in ("differential", "algebraic"):
            raise ValueError(f"Invalid kind: {kind}")
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("scale must be positive and finite")
        shape = tuple(shape)
        if any(not isinstance(n, int) or n <= 0 for n in shape):
            raise ValueError("shape entries must be positive integers")
        arr = np.asarray(initial, dtype=float)
        if shape:
            if arr.size == 1:
                arr = np.full(shape, float(arr.item()))
            elif arr.shape != shape:
                raise ValueError(f"Initial {name}: expected shape={shape}, got {arr.shape}")
        else:
            if arr.size != 1:
                raise ValueError(f"Initial {name}: scalar expected")
            arr = arr.reshape(())
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"Initial {name} contains non-finite entries")
        sz = int(np.prod(shape)) if shape else 1
        spec = VariableSpec(name, slice(self._size, self._size+sz), shape,
                            kind, unit, arr.copy(), scale)
        self._variables[name] = spec
        self._size += sz
        return spec

    def freeze(self):
        self._frozen = True

    @property
    def size(self):
        return self._size

    @property
    def names(self):
        return tuple(self._variables)

    def spec(self, name: str) -> VariableSpec:
        return self._variables[name]

    def get(self, vector: np.ndarray, name: str):
        spec = self.spec(name)
        values = np.asarray(vector)[spec.indices]
        return values.reshape(spec.shape) if spec.shape else values[0]

    def initial_state(self):
        y = np.empty(self.size, dtype=float)
        for spec in self._variables.values():
            y[spec.indices] = spec.initial.ravel()
        return y

    def differential_mask(self):
        mask = np.zeros(self.size, dtype=bool)
        for spec in self._variables.values():
            if spec.kind == "differential":
                mask[spec.indices] = True
        return mask

    def scale_vector(self):
        scale = np.ones(self.size)
        for spec in self._variables.values():
            scale[spec.indices] = spec.scale
        return scale
