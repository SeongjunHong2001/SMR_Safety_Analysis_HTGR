# 공통 인터페이스 계약 v3 (Role 5, G0)

개발 기준 `main`: `d035aaed`. 계약 v3 구현은 [PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)에서 병합 SHA `52016a8`로 `main`에 반영됐다.
이 계약은 교육용 MHTGR-350 해석 코드의 **조립 경계**다. 실제 물성, 상관식,
형상, 사고 이력은 아직 구현 또는 검증되지 않았다. G0의 실행 증거는
`tests/contracts/`의 합성 모델과 기존 데모다. 이 문서에서 `계약`으로 적은
실제 장치 포트는 담당 물리 모델이 구현되기 전까지 구성요소에 존재하지 않는다.

## 1. 전역 DAE와 Component

`SystemAssembly`는 대표 Component와 연결식을 고정된 순서로 등록하고
`F(t, y, ydot)=0` 하나를 DASSL에 전달한다. Component는
`register(registry)`, `residual(t,y,ydot,context)`, `outputs(t,y,context)`를
제공한다. `name`, `key(variable)`, `ports()`, `residual_ids()`,
`residual_scales()`가 공통 API다. `key(variable)`은
`{instance}.{variable}`을 만든다. 배열 잔차는 각 스칼라에 고정 ID를
부여한다. 구형 Component가 ID 또는 scale을 제공하지 않으면
조립기가 순서 기반 ID와 1의 scale을 생성한다. 새 물리 모델은 둘 다 명시한다.

`residual`은 현재 trial의 전역 `y`, `ydot`만 읽고 1차원 유한 NumPy 배열을
반환한다. 상태 변경, 파일 쓰기, 내부 Newton/BDF, 이전 평가값 캐시를 하지
않는다. Component마다 상태 수와 잔차 수가 같을 필요는 없다. 전체만
`n_unknown = n_residual`이어야 하며, 이는 정칙성 보증이 아니다.
`SystemAssembly.residual_details()`는 잔차·ID·`s_R`를 같은 순서로 반환한다.
`diagnose_jacobian(t,y,ydot,cj)`는
`diag(1/s_R) (F_y+cj F_ydot) diag(s_y)`의 rank와 조건수를 보고한다.
`s_y`는 `StateRegistry.register(scale=...)`의 상태 척도이고 `s_R`은
방정식 단위별 잔차 척도다. 실제 계통의 scale은 각 담당자가 정해야 한다.

`StateRegistry`에는 독립 Differential/Algebraic만 한 번 등록한다.
Parameter, External Input, Derived, Alias, Residual, Diagnostic은 y에
추가하지 않는다. `registry.alias(name,target)`은 새 미지수를 만들지 않고
등록 완료 후 `freeze()`에서 참조를 검사한다. 포트는 반드시 등록된 **스칼라**
키를 참조한다. 유닛은 SI: 절대압 Pa, T K, h J/kg, 열률 W, 유량 kg/s.
전역 `p_sys`는 한 키만 사용하고 지역 포트/셀 압력은 참조하거나 별도
물리식이 있는 경우에만 독립 등록한다.

`SimulationContext(registry, config, scenario, providers, heat_rates)`는
읽기 전용 설정·시나리오·주입 provider를 보유한다. `heat_rates`는
`with_heat_rates()`로 현재 평가용 새 context에만 넣는다. 기존 데모와
콜백은 `context=registry`도 계속 사용한다. Assembly에
`SimulationContext` 템플릿을 전달하면 등록을 마친 Assembly Registry로
다시 바인딩한다. provider 선택은 builder/config
단계에서 완료하고 residual 도중 바꾸지 않는다. 실제 builder와 필수 입력
검증은 Role 1/물리 모델 연결 뒤 구현한다. fixture는 자동으로 선택하지 않는다.

## 2. 유체·열 접속

`FluidStreamPort(p, m_dot, h_outflow, h_instream, medium="helium")`의 앞
네 값은 Registry 스칼라 키이고 `medium`은 접속 검증용 메타데이터다.
포트 `m_dot > 0`은 **해당 Component로 유입**이다. `h_outflow`는
Component 출구 상태의 별칭 또는 정의 잔차가 있는 대수변수이고,
`h_instream`은 접속식이 정하는 대수변수다. 서로 다른 매질은 연결하지
않는다. 셀/포트 별칭을 새 독립 상태로 중복 등록하지 않는다.

N포트 이상 접속은 압력연속 `N-1`, 유량합 `1`, 각 포트의 유입 stream
혼합식 `N`개로 **2N 잔차**다. 각 `h_instream`은 자기 포트를 제외한
공급 포트(`m_dot < 0`)의 `h_outflow`로 계산한다. 정체·역류 근처는
양의 `m_dot_epsilon`으로 완화한다. 모든 유량이 0이면 다른 출구들의
평균을 택하지만 이류 열전달은 0이다. `StreamJunction`은 고정 잔차 ID,
압력·유량·엔탈피 별 scale을 받는다. `m_dot_epsilon`은 기존
`eps_mdot_rel`과 검증된 접속부 정격 유량으로 산출한다. 실제 정격 유량이
없으면 숫자를 정하지 않는다. 압력 손실은 접속식에 숨기지 않고 장치의
운동량식에 둔다. 기존 `FluidPort`의 두 포트 식은 데모 호환용이다.

`ThermalPort(T,Q_dot)`의 기존 이상 연결은 같은 온도 노드의 등식과
열률합 0을 사용한다. 온도가 다른 Core→RPV, RPV→RCCS 경계에는
이상 온도 동일식을 쓰지 않는다. `ThermalLink`는 주입한 순수
`transfer(t,y,context)`를 trial마다 한 번 호출해 `to_a`, `to_b`
열률 배열을 만들고 `sum(to_a)+sum(to_b)=0`을 검증한다. 열률은 기본
Derived이며 별도 대수 상태라면 독립 변수와 정의 잔차를 같이 추가한다.

### 대표 장치의 계약 포트

아래는 구현 담당자가 `ports()`로 공개할 이름과 매질/경계다. 아직
미구현인 대표 클래스는 `register()`에서 명시적 예외를 낸다.

| 대표 장치 | 공개 이름 | 형식/매질 | 조건·경계 | 담당 |
| --- | --- | --- | --- | --- |
| Core | `fluid.helium_inlet`, `fluid.helium_outlet` | FluidStreamPort/He | 필수, 두 포트 모두 유입 양수 | 2/5 |
| Core | `solid.outer_surface` | 온도 키 + ThermalLink/열 | RPV 경로 선택 시 필수; 고체 격자와 매핑 | 3/5 |
| HeatExchanger | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수 | 4 |
| HeatExchanger | `secondary_inlet`, `secondary_outlet` | FluidStreamPort/water | SG-2에서는 필수, SG-0 단순 열침원에서는 없음 | 4 |
| Compressor | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수; 축속도는 별도 상태/입력 | 2 |
| Piping | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수 | 2 |
| Plenum | `helium_ports[i]` | FluidStreamPort/He | 실제 접속구 수와 방향에 맞춰 2개 이상 | 2 |
| Vessel | `inner_surface`, `outer_surface` | 온도 키 + ThermalLink/열 | 구조 열 모델 선택 시 필수; 별도 노드/별칭 구분 | 4 |
| RCCS | `panel_surface` | 온도 키 + ThermalLink/열 | RCCS-1 열 경계 필수 | 4 |
| Boundary | `fluid` | FluidStreamPort/명시 매질 | 외부 유체 경계가 있을 때만 | 5 |

Helium과 water는 같은 FluidStreamPort 접속에 섞지 않는다. 각 포트의
압력·유량·stream은 해당 Component가 등록한 키만 참조한다. 실제
장치별 잔차 수·ID·scale·초기조건은 모델 PR에서 ledger로 검증한다.

## 3. 물성·유동·열전달 provider

`interfaces.contracts.HeliumProperties`는
`temperature(p,h)`, `enthalpy(p,T)`, `density(p,h)`,
`viscosity(p,T)`, `conductivity(p,T)`, `cp(p,T)`, `cv(p,T)`,
`drho_dp(p,h)`, `drho_dh(p,h)`를 요구한다. 입력은 SI, 스칼라 또는
NumPy broadcasting 가능한 배열이다. 출력은 broadcast shape의 값이다.
실제 EOS의 온도·압력 유효 범위 및 도함수 일관성은 Role 2가 정해 검증한다.
현재 생산용 구현은 없으며 시험용 물성만 `tests/fixtures/`에 있다.

Role 2 공통 flow/pressure_loss 계약은 셀량 shape `(N,)`, 면 유량
shape `(N+1,)`, 기하학적 축 방향 유량 양수, Component 포트의 유입
양수로 둔다. 경계 면은 포트 stream과 일치시킨다. 압력손실은 절대압
Pa에서 계산하며 마찰·형상·중력·가속도 항의 부호와 저유량/역류
적용범위를 모델 PR에 기록한다. 실제 kernel 시그니처와 수식은 Role 2
자료 확인 전까지 고정하지 않았다.

공통 열전달 API는 `physics.heat_transfer.base`와 동일 경로 패키지
`physics.heat_transfer`에서 제공한다.

| 함수 | 인자와 단위 | 반환/부호 |
| --- | --- | --- |
| `conduction(T_a,T_b,k,A,L)` | K,K,W/(m·K),m²,m | a→b 양수 W |
| `convection(T_a,T_b,h_conv,A)` | K,K,W/(m²·K),m² | a→b 양수 W |
| `radiation(T_a,T_b,epsilon,A,sigma)` | K,K,1,m²,W/(m²·K⁴) | a→b 양수 W; 절대 K |
| `heat_transfer_coefficient(Nu,k,D_h)` | 1,W/(m·K),m | W/(m²·K) |

상관식 provider는 `coefficient(state,geometry)->h_conv`를 구현한다.
`state`에는 현재 trial의 `T_wall`, `T_bulk`, `p`, `rho`, `mu`, `k`,
`cp`, `m_dot`을, `geometry`에는 `D_h`, `A`, `L`, 방향을 명시한다.
벽/벌크의 역할은 상관식 입력이고 일반 열률 함수의 a→b 부호와 별개다.
Re/Pr는 provider가 이 상태·기하에서 계산한다. 유효 Re/Pr/온도 범위,
층류·천이·자연대류·저유량·역류 처리와 shape는 각 Role 3/4 구현이
검증하고 문서화한다. 기본 난류 상관식이나 유효 범위를 넘는 값은
자동 선택하지 않는다. `mode_nu`, `include_nc`, `mode_friction`,
`mode_he_eos`, `mode_sg_he_crossflow` 등 기존 mode를 선택에 쓴다.
현재 base 함수는 구현 파일 경로를 명시한 `NotImplementedError`를 낸다.

## 4. 복합 Core와 열률 단일 경로

`Core(config=None,properties=None,*,fluid,solid,power,coupling,name="core")`
는 기존 생성자 인자와 import를 보존한다. 실제 사용에는 네 구현을
생성자에 주입한다. Core만 SystemAssembly에 등록되고, 등록 순서는
fluid→solid→상태가 있는 power다. power가 독립 점동특성 상태를
소유하면 `register`와 `residual`을 함께 제공한다. 상태 없는 지정
출력 provider는 둘 다 제공하지 않는다. `coupling.evaluate(t,y,context)`는 같은 trial에서
`(Q_to_fluid,Q_from_solid)`를 W 배열로 반환한다. 길이가 달라도 합이
같아야 한다. `power.evaluate(t,y,context)`는 고체 격자 W 배열을 반환한다.
Core는 fluid에 `heat_rate=Q_to_fluid`, solid에
`heat_rate=power-Q_from_solid`를 전달한다. 하위 모델은 이 열률로
각자 잔차를 구성하며 전역 캐시와 호출 순서에 의존하지 않는다.
상관식은 `physics/heat_transfer/` 한 곳에서 계산하고
`models/core/coupling.py`에는 격자·면적·열률 매핑만 둔다. SG도 같은
단일 경로를 사용한다.

## 5. 폐루프 압력·재고량 ledger

실제 `p_sys`의 방정식, 국부압 기준점, `M_he_loop_initial`은 Role 2/5
공동 결정 사항이다. 셀 질량식을 합한 `R_he_loop_mass`는 진단값이며
독립 잔차가 아니다. 정상 초기화에서 시간미분이 0이 되어 종속되는
질량식이 확인되면, 그 **한 식을** 기준압 또는 검증된 재고량 조건으로
대체한다. 물리식에 임의 방정식을 더하거나 빼지 않는다.

`tests/fixtures/closed_loop.py`의 5미지수 합성 루프는 다음을 시험한다.

| 단계 | 독립 미지수 | 잔차 | 결합 Jacobian rank |
| --- | ---: | --- | ---: |
| 과도 `cj=1` | 5 | 운동량 2 + 질량 연속 1 + 국부 기준 1 + p_sys 질량동역학 1 | 5 |
| 정상 `cj=0`, 동역학 식 유지 | 5 | 같은 5식이지만 시간미분 항 0 | 4 |
| 정상, 종속 식을 재고량 조건으로 대체 | 5 | 운동량 2 + 질량 연속 1 + 국부 기준 1 + 재고량 1 | 5 |
| 전역 질량식을 중복 추가 | 5 | 6식 | 비정방 |

이 숫자는 fixture의 증거일 뿐 실제 MHTGR-350의 자유도·index-1
판정이 아니다. 실제 셀/포트 수, 압축성, EOS, 경계와 정상 초기화를
정한 뒤 동일 ledger와 scaled Jacobian을 다시 검증해야 한다.

## 6. 설정과 완료 범위

파일별 키·단위·shape·필수 조건은 [설정 스키마](config_schema.md)에,
변수와 식의 변경안은 `docs/registry_patches/role5/`에 둔다.
기본 YAML의 `null`은 미정이지 0이 아니다. 실제 물성/형상/사고
이력이 없으면 경로와 항목을 알려주는 오류로 중단해야 한다.
현재 G0는 합성 결합 시험과 기존 데모까지만 통과했으며 실제 HTGR
정상상태 또는 LOFC 결과를 제공하지 않는다.
