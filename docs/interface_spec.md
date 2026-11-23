# Component Interface Specification (초안, v0.1)

## 당장 고정할 규약

- 각 모델은 `interfaces.component.Component`의 `register(registry)`와
  `residual(t, y, ydot, context)`를 구현한다.
- `residual`은 현재 전역 Newton 추정값에서 **1차원 NumPy 배열**을 반환한다.
- `residual` 안에서 `y`, `ydot`, 시뮬레이션 시각이나 모듈 내부 상태를 변경하지 않는다.
- 각 미지수는 딱 한 번 등록한다. 공유 변수는 상대방 Registry 이름을 참조한다.
- 내부 Cell/Face/FVM/FEM 이산화 및 해당 Residual 계산은 각 물리 담당자의 책임이다.
- 시간 이력 저장, BDF, Newton, 스텝 수용 및 거부는 Solver 담당자(Role 5)의 책임이다.
- 단위는 SI. `p`: 절대압력 Pa, `T`: K, `h`: J/kg, `m_dot`: kg/s.
- 유체 포트 질량유량은 **해당 Component로 유입하는 방향이 양(+)**이다.
- 역할 5가 조립하는 전체 잔차는 **전체 미지수 수와 길이가 같아야 한다.** 개수만 같아도
  독립성/일치 초기조건은 보장되지 않으므로 Jacobian 및 물리 보존 테스트가 필요하다.
- 포트 연결식은 `SystemAssembly(..., connections=...)`에 `(y, registry)`를 받는
  콜백으로 전달한다. 조립 순서는 Component 잔차, 연결 잔차 순이다. 예를 들어
  `functools.partial(ideal_thermal_two_port_residual, a=port_a, b=port_b)`를 사용할 수 있다.

## 유체 Stream 포트와 다분기 접속

`FluidStreamPort(p, m_dot, h_outflow, h_instream)`는 압력, 질량유량,
구성요소에서 나가는 유체의 엔탈피, 구성요소로 들어오는 유체의 엔탈피를
서로 다른 Registry 이름으로 참조한다. 각 이름은 스칼라로 등록한다.
`h_outflow`는 구성요소가 계산하는 상태이고, `h_instream`은 연결식이 결정하는
대수변수다.
구성요소는 에너지식에서 유입 유량에는 `h_instream`, 유출 유량에는
`h_outflow`를 사용한다. 기존 `FluidPort(p, h, m_dot)`는 두 포트 연결 예제에만 사용한다.

`ideal_fluid_stream_junction_residual`은 N개 포트에 대해 압력 동일 조건 N−1개,
유량 합 0 조건 1개, 포트별 유입 엔탈피 조건 N개를 반환한다. 각 포트로 들어갈
유체의 엔탈피는 **자기 포트를 제외한** 분기점 공급 포트들의 `h_outflow`를
질량유량으로 가중해 계산한다. 공급 포트는 `m_dot < 0`인 포트다.
이상 접속부 자체의 압력 손실은 0으로 둔다. 배관·장치의 손실은 해당
Component가 계산한다. 정체·역류 부근에는 양의 `mass_flow_epsilon`(kg/s)으로
유입 가중치를 매끄럽게 한다. 유량이 모두 0이면 다른 포트들의 출구
엔탈피 평균을 유입 엔탈피로 둔다. 이때 실제 이류 에너지 전달은 0이다.

```python
from functools import partial
from interfaces.ports import FluidStreamPort
from system.connections import ideal_fluid_stream_junction_residual

a = FluidStreamPort("a.p", "a.m_dot", "a.h_outflow", "a.h_instream")
b = FluidStreamPort("b.p", "b.m_dot", "b.h_outflow", "b.h_instream")
c = FluidStreamPort("c.p", "c.m_dot", "c.h_outflow", "c.h_instream")
connection = partial(
    ideal_fluid_stream_junction_residual,
    ports=(a, b, c), mass_flow_epsilon=1e-6,
)
```

각 포트 이름은 먼저 Component의 `register`에서 등록해야 한다. 연결 콜백은
`SystemAssembly(..., connections=[connection])`에 전달한다. 전역 미지수와
잔차의 수를 맞추려면 구성요소와 연결식의 방정식 수를 함께 계산한다.

## 아직 확정하지 않은 항목

- MHTGR-350 통합변수표에 대응하는 실제 사고/검증 데이터, Cell/Face 배치, 모델별 독립 미지수 목록
- Heat transfer model, 핵분열/붕괴열 모델, EOS/물성, Component별 유량역전 응답
- Junction 압력손실의 실제 모델과 `mass_flow_epsilon`의 계통별 선정값
- 전역 Residual의 물리 단위별 Scaling, 희소 Jacobian 정책

## 새 Component 제출 체크

1. `docs/variable_registry.csv`: 변수·단위·shape·초기값·소유자
2. `docs/equation_registry.csv`: 지배방정식·공간 이산화·잔차 수
3. Python Component의 `register`, `residual`, `outputs` 구현
4. 단독 테스트 및 열/질량 보존 확인

## 예제

`python main.py --demo`는 실제 HTGR이 아닌 2열용량+1대수 열유량 모형을 실행한다.
그 목적은 전역 y/ydot, 서로 공유하는 변수, DAE의 초기화와 시간적분을 독립적으로 검증하는 것이다.
