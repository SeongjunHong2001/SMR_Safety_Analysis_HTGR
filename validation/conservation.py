"""전역 보존 오차 평가 (단순 유틸리티)."""

def balance_residual(accumulation_rate, source_rate, sink_rate):
    """d(저장량)/dt - 유입원 + 유출원. 단위 일치가 선행되어야 함."""
    return accumulation_rate - source_rate + sink_rate
