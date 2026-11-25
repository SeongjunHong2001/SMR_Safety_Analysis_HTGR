"""RPV 전도·복사·외부 냉각 경로.

TODO: 담당자가 변수 및 공간 이산화 모델 확정 후 완성.
"""
from interfaces.component import Component

class Vessel(Component):
    def __init__(self, config=None, properties=None):
        self.config = config or {}
        self.properties = properties

    def register(self, registry):
        raise NotImplementedError("Vessel: 변수/방정식 확정 후 등록")

    def residual(self, t, y, ydot, context):
        raise NotImplementedError("Vessel: 공간 이산화 잔차 방정식 구현 필요")
