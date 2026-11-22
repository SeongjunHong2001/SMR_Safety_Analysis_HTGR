"""노심의 1D 유동 및 고체 열전달·출력 연성.

TODO: 담당자가 변수 및 공간 이산화 모델 확정 후 완성.
"""
from interfaces.component import Component

class Core(Component):
    def __init__(self, config=None, properties=None):
        self.config = config or {}
        self.properties = properties

    def register(self, registry):
        raise NotImplementedError("Core: 변수/방정식 확정 후 등록")

    def residual(self, t, y, ydot, context):
        raise NotImplementedError("Core: 공간 이산화 잔차 방정식 구현 필요")
