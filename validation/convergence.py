"""시간/공간 수렴성 비교의 기본 지표."""
import numpy as np

def max_relative_difference(reference, candidate, floor=1e-12):
    a = np.asarray(reference, dtype=float)
    b = np.asarray(candidate, dtype=float)
    return float(np.max(np.abs(a-b)/np.maximum(np.abs(a),floor)))
