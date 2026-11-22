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
- FluidPort 질량유량은 **해당 Component로 유입하는 방향이 양(+)**이다.
- 역할 5가 조립하는 전체 잔차는 **전체 미지수 수와 길이가 같아야 한다.** 개수만 같아도
  독립성/일치 초기조건은 보장되지 않으므로 Jacobian 및 물리 보존 테스트가 필요하다.
- 포트 연결식은 `SystemAssembly(..., connections=...)`에 `(y, registry)`를 받는
  콜백으로 전달한다. 조립 순서는 Component 잔차, 연결 잔차 순이다. 예를 들어
  `functools.partial(ideal_thermal_two_port_residual, a=port_a, b=port_b)`를 사용할 수 있다.

## 아직 확정하지 않은 항목

- 대상 노형, 실제 사고/검증 데이터, Cell/Face 배치, 동적/대수 변수 목록
- Heat transfer model, 핵분열/붕괴열 모델, EOS/물성, 유량역전 모델
- Junction 압력손실 및 다분기 Stream Mixing 규약
- 전역 Residual의 물리 단위별 Scaling, 희소 Jacobian 정책

## 새 Component 제출 체크

1. `docs/variable_registry.csv`: 변수·단위·shape·초기값·소유자
2. `docs/equation_registry.csv`: 지배방정식·공간 이산화·잔차 수
3. Python Component의 `register`, `residual`, `outputs` 구현
4. 단독 테스트 및 열/질량 보존 확인

## 예제

`python main.py --demo`는 실제 HTGR이 아닌 2열용량+1대수 열유량 모형을 실행한다.
그 목적은 전역 y/ydot, 서로 공유하는 변수, DAE의 초기화와 시간적분을 독립적으로 검증하는 것이다.
