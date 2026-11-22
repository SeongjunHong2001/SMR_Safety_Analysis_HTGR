"""검증 데이터 비교 (동일한 시간격자에 맞춘 뒤 사용)."""
import numpy as np

def rmse(predicted, observed):
    p = np.asarray(predicted, dtype=float)
    o = np.asarray(observed, dtype=float)
    if p.shape != o.shape:
        raise ValueError("Predicted and observed series must have matching shapes")
    return float(np.sqrt(np.mean((p-o)**2)))
