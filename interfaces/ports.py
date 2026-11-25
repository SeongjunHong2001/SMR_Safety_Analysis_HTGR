"""Port 심볼만 정의. 값은 항상 StateRegistry와 전역 y에서 읽는다.

유량은 각 Component에 들어가는 방향이 양(+)이다. 압력 p는 절대압력[Pa],
엔탈피 h는 비엔탈피[J/kg]. FluidPort는 기존 두 포트 연결용이며,
FluidStreamPort는 유동 역전과 다분기 혼합용이다.
"""
from dataclasses import dataclass

@dataclass(frozen=True)
class FluidPort:
    p: str       # StateRegistry에 등록된 전역 변수명
    h: str       # StateRegistry에 등록된 전역 변수명
    m_dot: str   # StateRegistry에 등록된 전역 변수명


@dataclass(frozen=True)
class FluidStreamPort:
    """분기점에 연결할 포트. 유량은 구성요소로 유입할 때 양수다.

    h_outflow는 구성요소에서 분기점으로 나가는 유체의 상태이고,
    h_instream은 분기점에서 구성요소로 들어오는 유체의 상태다.
    네 이름은 Registry에 스칼라로 등록한다. h_outflow는 구성요소가
    소유한 상태를 참조할 수 있고, 연결식이 대수변수 h_instream을 결정한다.
    """

    p: str
    m_dot: str
    h_outflow: str
    h_instream: str
    medium: str = "helium"  # connection metadata, never a Registry state


@dataclass(frozen=True)
class ThermalPort:
    T: str
    Q_dot: str   # Component에 들어오는 열량을 양(+)으로
