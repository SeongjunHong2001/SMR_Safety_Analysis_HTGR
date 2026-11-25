"""RCCS component boundary pending Role 4 air-loop model and inputs."""
from interfaces.component import Component


class RCCS(Component):
    def __init__(self, config=None, air_properties=None):
        self.config = {} if config is None else config
        self.air_properties = air_properties

    def register(self, registry):
        raise NotImplementedError("RCCS requires models/rccs/rccs.py Role 4 implementation and air-path inputs")

    def residual(self, t, y, ydot, context):
        raise NotImplementedError("RCCS air-loop residual requires Role 4 implementation")
