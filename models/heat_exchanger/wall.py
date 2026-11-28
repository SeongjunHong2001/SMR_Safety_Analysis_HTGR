"""Role 4 SG 관벽.

``SGWall`` 계약과 SG-0용 ``NoStorageWall``. SG-0은 관벽 축열을 풀지 않는다 (팀 결정 2026-11-26: 1차계통
해석 코드이므로 SG는 열이 빠져나가는 것만 모사). 결합 모델이 He → 벽, 벽 → 물 열률을 같은 크기로 주므로
벽에 남는 순열률은 0이고, 이 벽은 상태·잔차가 없다. 순열률이 0이 아니면 오류를 낸다 (에너지 누락 방지).
SG-1/2 관벽 축열(T_sg_wall)은 아직 구현하지 않았다.
"""
from typing import Protocol

import numpy as np


class SGWall(Protocol):
    def register(self, registry): ...
    def residual(self, t, y, ydot, context, *, heat_rate):
        """Return 1-D wall residuals with net current-trial W heat input."""
        ...
    def residual_ids(self): ...
    def residual_scales(self): ...
    def outputs(self, t, y, context): ...


class NoStorageWall:
    """SG-0 관벽: 상태 없음, 순열률 0 확인만 한다."""

    def __init__(self, name="sg.wall", n_cells=None):
        self.name = name
        self.n_cells = n_cells

    def register(self, registry):
        pass

    def residual(self, t, y, ydot, context, *, heat_rate):
        q = np.asarray(heat_rate, dtype=float)
        scale = max(float(np.max(np.abs(q))) if q.size else 0.0, 1.0)
        if self.n_cells is not None and q.shape != (self.n_cells,):
            raise ValueError(f"{self.name}: wall heat rate must have shape ({self.n_cells},)")
        if np.any(np.abs(q) > 1e-9 * scale + 1e-6):
            raise ValueError(f"{self.name}: SG-0 wall has no heat storage but received net heat {q.sum():.3e} W")
        return np.empty(0)

    def residual_ids(self):
        return ()

    def residual_scales(self):
        return np.empty(0)

    def outputs(self, t, y, context):
        return {}
