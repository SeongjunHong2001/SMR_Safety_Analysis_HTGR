# 공통 인터페이스 계약 v3 (Role 5, G0)

`d035aaed`는 계약 작성 당시의 기준이다. 계약 v3 구현은 [PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)에서 병합 SHA `52016a8`로 `main`에 반영됐다. 이후 개발은 최신 `main`을 확인한다.
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
등록 완료 후 `freeze()`에서 참조를 검사한다. 유체 포트와 `ThermalPort`는
등록된 **스칼라** 키를 참조한다. `ThermalSurfacePort`는 배열 상태의
검증된 한 셀도 `index`로 가리킬 수 있다. 유닛은 SI: 절대압 Pa, T K,
h J/kg, 열률 W, 유량 kg/s.
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

`SystemAssembly(components, *, context=None, connections=(), thermal_links=())`에서
`thermal_links`는 `ThermalLink` 객체의 열이다. 링크가 있으면 템플릿
context는 `SimulationContext`여야 한다(`None`이면 조립기가 생성한다).
각 잔차 평가와 `outputs()`에서 링크의 `transfer(t,y,context)`를 링크별로
정확히 한 번 호출한다. 현재 trial의 열률은 `context.heat_rates`의
`{link.name}.to_a`, `{link.name}.to_b`에 읽기 전용 W 배열로 전달한다.
이 키는 다른 링크 또는 템플릿 context와 중복될 수 없고, 다음 평가로
이전 열률을 넘기지 않는다. 각 장치는 해당 배열을 자기 에너지 잔차에
한 번만 반영하며 배열의 격자 매핑은 장치 담당자가 소유한다.

## 2. 유체·열 접속

`FluidStreamPort(p, m_dot, h_outflow, h_instream, medium="helium")`의 앞
네 값은 Registry 스칼라 키이고 `medium`은 접속 검증용 메타데이터다.
포트 `m_dot > 0`은 **해당 Component로 유입**이다. `h_outflow`는
Component 출구 상태의 별칭 또는 정의 잔차가 있는 대수변수이고,
`h_instream`은 접속식이 정하는 대수변수다. 서로 다른 매질은 연결하지
않는다. 같은 접속부의 어느 `h_outflow`와도 같은 상태를 가리키면 안 된다.
셀/포트 별칭을 새 독립 상태로 중복 등록하지 않는다.

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
열률합 0을 사용한다. 이 연결에 쓰는 `Q_dot`은 Registry에 등록된
스칼라 **대수변수**여야 하며 해당 Component가 정의 잔차를 제공한다.
열률이 Derived인 경로는 `ThermalPort`에 미등록 키를 넣지 않고
`ThermalSurfacePort(T,index=None)`로 온도 경계만 공개하고
`ThermalLink`로 계산한다. `T`는 등록된 K 상태다. 스칼라 상태는
`index=None`, 배열 상태의 한 셀은 `index` 정수 또는 인덱스 튜플로
지정한다. `port.temperature(y,registry)`가 그 스칼라 K 값을 읽는다.
온도가 다른 Core→RPV, RPV→RCCS 경계에는
이상 온도 동일식을 쓰지 않는다. `ThermalLink`는 주입한 순수
`transfer(t,y,context)`를 trial마다 한 번 호출해 `to_a`, `to_b`
비어 있지 않은 1-D 열률 배열을 만들고 `sum(to_a)+sum(to_b)=0`을 검증한다. 열률은 기본
Derived이며 별도 대수 상태라면 독립 변수와 정의 잔차를 같이 추가한다.
변수표의 `Q_dot_thermal_port` Derived 분류는 기본 열경로를 가리킨다.

### 대표 장치의 계약 포트

아래는 구현 담당자가 `ports()`로 공개할 이름과 매질/경계다. 아직
미구현인 대표 클래스는 `register()`에서 명시적 예외를 낸다.

| 대표 장치 | 공개 이름 | 형식/매질 | 조건·경계 | 담당 |
| --- | --- | --- | --- | --- |
| Core | `fluid.helium_inlet`, `fluid.helium_outlet` | FluidStreamPort/He | 필수, 두 포트 모두 유입 양수 | 2/5 |
| Core | `solid.outer_surface` 또는 `solid.outer_surface[i]` | ThermalSurfacePort/K + ThermalLink/열 | RPV 경로 선택 시 필수; 고체 격자와 매핑 | 3/5 |
| HeatExchanger | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수; 상위 Composite가 하위 포트를 그대로 공개 | 4 |
| HeatExchanger | `secondary_inlet`, `secondary_outlet` | FluidStreamPort/water | SG-2에서는 두 포트 모두 필수, SG-0/1 열침원에서는 없음 | 4 |
| Compressor | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수; 축속도는 별도 상태/입력 | 2 |
| Piping | `helium_inlet`, `helium_outlet` | FluidStreamPort/He | 필수 | 2 |
| Plenum | `helium_ports[i]` | FluidStreamPort/He | 실제 접속구 수와 방향에 맞춰 2개 이상 | 2 |
| Vessel | `inner_surface`, `outer_surface` | ThermalSurfacePort/K + ThermalLink/열 | 구조 열 모델 선택 시 필수; 별도 노드/별칭 구분 | 4 |
| RCCS | `panel_surface` | ThermalSurfacePort/K + ThermalLink/열 | RCCS-1 열 경계 필수 | 4 |
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

`HeatExchanger(config=None,properties=None,*,helium_side=None,wall=None,secondary_side=None,coupling=None,name="heat_exchanger")`는 구현된 상위 조립기다. 실제 구현 네 개를 생성자에 명시적으로 주입하고 상위 객체만 Assembly에 등록한다. 자식의 `register`, `residual`, `residual_ids`, `residual_scales` 순서는 helium_side→wall→secondary_side이며 `outputs`를 병합한다. `helium_side`와 `secondary_side`의 `ports()`에서 외부 포트를 공개한다. 상세 모델이 없으면 `register`에서 명시적으로 실패하며 fixture를 자동 선택하지 않는다. Role 4는 세 자식과 coupling의 실제 물리 구현을 맡는다.

`coupling.evaluate(t,y,context)`는 같은 trial에서 `(to_helium,to_wall_from_helium,to_wall_from_secondary,to_secondary)`의 비어 있지 않은 유한 1-D W 배열을 반환한다. 첫째·둘째 배열의 합, 셋째·넷째 배열의 합은 각각 0이고 두 벽 배열의 shape는 같다. 상위는 helium_side에 첫째, wall에 둘째+셋째, secondary_side에 넷째를 `heat_rate=`로 전달한다. 자식은 받은 열률을 에너지 잔차에 한 번만 넣는다. Role 4 coupling은 실제 열전달 계산에 `physics.heat_transfer` 공개 API를 사용하고 장치별 격자·면적 매핑만 소유한다. 기본 열전달 커널과 상세 물리식은 아직 미구현이다. 합성 조립·flat/Composite·DASSL 증거는 `tests/contracts/test_heat_exchanger_composite.py`에 있다.
SG-0/1의 단순 2차 열침원도 `secondary_side` 인스턴스를 명시적으로 주입한다. 이 경계 구현은 물 포트와 독립 상태·잔차가 없을 수 있지만 `register`, `residual(...,heat_rate=...)`, ID·scale·`ports`·`outputs` 계약을 지켜 외부로 빠지는 열률을 기록한다.

## 5. 폐루프 압력·재고량 ledger

기본 저마하수 계약의 전역 `p_sys`는 **Role 2가 구현할 지정 lower-plenum 인스턴스만** `register(registry)`에서 `name="p_sys", kind="differential", owner="primary_loop"`로 한 번 등록한다. Role 5의 실제 builder는 해당 인스턴스를 하나만 지정하고 `p_sys` 등록 수·owner·kind를 검사한다. 다른 헬륨 장치는 이 키를 읽되 다시 등록하지 않는다. `h_he[i]`는 셀별 Differential, `p_fluid_cell[i]`와 관성항을 생략한 기본 모델의 `m_dot_fluid_face[j]`는 Algebraic이다. 후자의 모델에서 운동량 관성항을 채택하면 해당 state kind와 식을 별도 patch로 변경해야 한다. `rho_he[i]=rho(p_sys,h_he[i])`와 `M_he_loop`는 Derived이며 독립 상태가 아니다.

국부 운동량식은 `p_fluid_cell`의 압력차를 사용한다. 지정 lower-plenum 셀 **한 곳의 국부 압력 기준식** `R_p_ref=p_fluid_cell[lower]-p_sys=0`을 대수 잔차로 등록한다. 이 식은 정상 재고량 조건과 별개다. 기본 EOS는 국부 정압이 아니라 `p_sys`로 평가한다. 각 고정 체적 `V_i`의 국부 질량식은 다음과 같다. 면 유량의 유입 부호는 장치 안에서 일관되게 정한다.

```text
rho_i = rho(p_sys, h_i)
R_M,i = V_i * [drho_dp(p_sys,h_i) * p_sys_dot
             + drho_dh(p_sys,h_i) * h_i_dot]
        + sum(outgoing m_dot) - sum(incoming m_dot) = 0
M_he_loop = sum_i rho_i * V_i
```

헬륨 체적을 소유한 모든 Component의 담당자(Role 2의 Core/배관/플레넘/순환기와 Role 4의 SG 헬륨 측)는 `interfaces.contracts.HeliumInventory`의 `helium_volumes()->1-D m³`와 `helium_inventory(t,y,context)->1-D kg`를 같은 셀 순서로 제공한다. 후자는 현재 trial의 `rho(p_sys,h)*V`다. Role 5는 `system.inventory.total_helium_inventory(components,t,y,context)`로 플레넘·Core·배관·SG의 각 체적을 합산한다. 헬퍼는 동일 Component 중복을 거부하며 실제 서로 다른 장치 간 체적 중복은 builder의 기하 ledger에서 검사한다. `sum_i R_M,i = dM_he_loop/dt`는 닫힌 경계에서 진단값 `R_he_loop_mass`이며 추가 독립 잔차가 아니다.

과도 잔차는 모든 국부 질량식, 각 셀 엔탈피형 에너지식, 운동량식, 국부 압력 기준식을 유지한다. 정상 `ydot=0`에서는 폐루프 국부 질량식 한 개가 종속된다. 확인된 **한 ID**를 정상 초기화 전용 잔차에서 `M_he_loop-M_he_loop_initial=0`으로 교체한다. 검증된 초기 재고량이 없고 기준압만 있으면 그 한 자리에 `p_sys-p_sys_initial=0`을 사용하고 정상해에서 재고량을 유도한다. 미정값을 0으로 대입하지 않는다. 두 선택은 배타적이고 둘 다 국부 기준식 `R_p_ref`와 별개다. 정상해는 교체 전 과도 잔차에 `ydot=0`을 넣어 다시 확인한다. 정상용 제약을 과도 잔차에 유지하거나 전역 질량식을 추가하지 않는다.

기본 엔탈피형 에너지 잔차는 `rho_i*V_i*h_i_dot - V_i*p_sys_dot + sum_in[m_in*(h_i-h_upstream)] - Q_i=0`이다. 위 질량식과 함께 `U_i=V_i*(rho_i*h_i-p_sys)`의 보존형 내부에너지 수지와 동치인지 검산한다. 전체 셀에서 내부 면의 엔탈피·질량 전달을 상쇄하고 외부 열·일과 `d(sum U_i)/dt`를 비교한다. 이동 경계나 축일은 해당 물리 모델의 식과 ledger에 별도로 명시한다.

`tests/fixtures/low_mach_loop.py`의 2셀 합성 폐루프는 7미지수(`p_sys` 1, `h` 2, 국부압 2, 면유량 2)를 사용한다. 과도는 질량 2+에너지 2+운동량 2+국부 기준 1=7식, `cj=1` 결합 Jacobian rank 7이다. 정상에서 원래 식 그대로면 `cj=0` rank 6이고, 한 국부 질량식만 재고량 또는 기준압 제약으로 교체하면 7식 rank 7이다. 전역 질량식을 더하면 8식으로 비정방이다. 정상해·일치 초기화·질량/에너지 수지·DASSL 적분은 `tests/contracts/test_low_mach_closure.py`가 검사한다. 기존 `tests/fixtures/closed_loop.py`의 5미지수 시험도 회귀 증거다. 이는 **공통 수학 구조의 합성 검증**이며 실제 MHTGR-350 EOS·기하·경계의 Jacobian rank와 DASSL 적분은 아직 검증하지 않았다.

## 6. 설정과 완료 범위

파일별 키·단위·shape·필수 조건은 [설정 스키마](config_schema.md)에,
변수와 식의 변경안은 `docs/registry_patches/role5/`에 둔다.
기본 YAML의 `null`은 미정이지 0이 아니다. 실제 물성/형상/사고
이력이 없으면 경로와 항목을 알려주는 오류로 중단해야 한다.
현재 G0는 합성 결합 시험과 기존 데모까지만 통과했으며 실제 HTGR
정상상태 또는 LOFC 결과를 제공하지 않는다.
