"""계통 연결 이름과 방향을 저장. 흐름역전과 계산 순서는 독립."""
from dataclasses import dataclass

@dataclass(frozen=True)
class Edge:
    left: str
    right: str

class Topology:
    def __init__(self, edges=()):
        self.edges = tuple(edges)
        for e in self.edges:
            if not isinstance(e, Edge):
                raise TypeError("Topology edges must be Edge instances")
            if e.left == e.right:
                raise ValueError("Connection cannot be self-loop at same port")
