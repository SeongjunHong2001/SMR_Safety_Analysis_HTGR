"""두 열용량 + 공유 열유량의 3변수 Index-1 DAE 예제.

q는 대수미지수, T_hot/T_cold는 미분미지수.
실제 HTGR의 Core/HX 모델이 아님을 분명히 한다.
"""
import numpy as np
from interfaces.component import Component
from system.assembly import SystemAssembly

class DemoHeater(Component):
    def __init__(self, heat_capacity=20.0, heat_input=50.0):
        self.C = heat_capacity
        self.Q = heat_input

    def register(self, registry):
        registry.register("demo.hot_T", kind="differential", unit="K", initial=320.0, scale=300.0)

    def residual(self, t, y, ydot, context):
        reg = context
        Thdot = reg.get(ydot, "demo.hot_T")
        q = reg.get(y, "demo.q_dot")
        return np.array([self.C * Thdot - self.Q + q])

    def outputs(self, t, y, context):
        return {"hot_temperature_K": float(context.get(y, "demo.hot_T"))}

class DemoCooler(Component):
    def __init__(self, heat_capacity=40.0, sink_ua=3.0, sink_temperature=300.0):
        self.C = heat_capacity
        self.sink_ua = sink_ua
        self.Tsink = sink_temperature

    def register(self, registry):
        registry.register("demo.cold_T", kind="differential", unit="K", initial=300.0, scale=300.0)

    def residual(self, t, y, ydot, context):
        reg = context
        Tc = reg.get(y, "demo.cold_T")
        Tcdot = reg.get(ydot, "demo.cold_T")
        q = reg.get(y, "demo.q_dot")
        return np.array([self.C * Tcdot - q + self.sink_ua*(Tc-self.Tsink)])

    def outputs(self, t, y, context):
        return {"cold_temperature_K": float(context.get(y, "demo.cold_T"))}

class DemoHeatLink(Component):
    """두 열용량을 연결하는 대수식 q-UA*(T_hot-T_cold)=0."""
    def __init__(self, ua=5.0):
        self.ua = ua

    def register(self, registry):
        registry.register("demo.q_dot", kind="algebraic", unit="W", initial=0.0, scale=50.0)

    def residual(self, t, y, ydot, context):
        reg = context
        Th = reg.get(y, "demo.hot_T")
        Tc = reg.get(y, "demo.cold_T")
        q = reg.get(y, "demo.q_dot")
        return np.array([q-self.ua*(Th-Tc)])

    def outputs(self, t, y, context):
        return {"heat_transfer_W": float(context.get(y, "demo.q_dot"))}

def make_demo_system():
    components = [DemoHeater(), DemoCooler(), DemoHeatLink()]
    sys = SystemAssembly(components)
    sys.context = sys.registry
    return sys
