"""HX 1차 헬륨, 벽, 2차 경계/상태.

TODO: 담당자가 변수 및 공간 이산화 모델 확정 후 완성.
"""
from interfaces.component import Component

class HeatExchanger(Component):
    def __init__(self, config=None, properties=None):
        self.config = config or {}
        self.properties = properties

    def register(self, registry):
        raise NotImplementedError("HeatExchanger: 변수/방정식 확정 후 등록")

    def residual(self, t, y, ydot, context):
        raise NotImplementedError("HeatExchanger: 공간 이산화 잔차 방정식 구현 필요")
