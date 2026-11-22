"""Port 심볼만 정의. 값은 항상 StateRegistry와 전역 y에서 읽는다.

두 FluidPort의 m_dot는 각 Component에 **들어가는 방향이 양(+)**이다.
압력 p는 절대압력[Pa], h는 비엔탈피[J/kg].
유동역전/다분기 Junction은 별도의 stream/mixing 모델이 필요하다.
"""
from dataclasses import dataclass

@dataclass(frozen=True)
class FluidPort:
    p: str       # StateRegistry에 등록된 전역 변수명
    h: str       # StateRegistry에 등록된 전역 변수명
    m_dot: str   # StateRegistry에 등록된 전역 변수명

@dataclass(frozen=True)
class ThermalPort:
    T: str
    Q_dot: str   # Component에 들어오는 열량을 양(+)으로
