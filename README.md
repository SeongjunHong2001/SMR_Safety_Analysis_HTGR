# SMR Safety Analysis · HTGR

MHTGR-350의 정상상태와 가압 냉각재 유량 상실(P-LOFC/PCC)을
해석하기 위한 Python 프로젝트다. 구성요소의 상태와 보존식을 하나의
미분대수방정식(DAE) `F(t, y, ydot) = 0`으로 조립한다.

**현재 실행 가능한 계산은 두 열용량 데모와 합성 계약 시험이다.** 실제
원자로 모델에 필요한 물성, 장치별 잔차식, 형상과 사고 이력은 아직
확정되지 않았다. [구현 현황](docs/implementation-status.md)에서 기능별
범위를 확인할 수 있다.

## 설치와 첫 실행

Python 3.10 이상이 필요하다. 저장소를 내려받은 뒤 프로젝트 최상위 폴더에서 실행한다.

```bash
git clone https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR.git
cd SMR_Safety_Analysis_HTGR
python -m pip install -r requirements.txt
python main.py --demo
```

`main.py`의 기본 실행과 `--demo`는 같은 **합성 데모**를 실행한다.
고온부·저온부 온도는 미분변수, 열유량은 대수변수다. 기본 설정은
0초부터 5초까지 적분하며 수용·거절 스텝과 사용한 최고 차수를 출력한다.
초기 시점과 수용 시점의 결과는 `output/results/demo_dae.csv`에 저장된다.

| CSV 열 | 의미 |
| --- | --- |
| `time_s` | 시뮬레이션 시간(s) |
| `demo.hot_T_K`, `demo.cold_T_K` | 고온부와 저온부 온도(K) |
| `demo.q_dot_W` | 고온부에서 저온부로 전달되는 열유량(W) |

적분 구간, 첫 스텝 후보 간격(`dt`), 오차 허용값(`rtol`, `atol`)은
[솔버 설정](config/solver.yaml)에 있다. 다른 설정은
`python main.py --demo --solver-config 경로/설정.yaml`로 지정한다.

```bash
python -m pytest -q
python -m pytest -q tests/contracts
python scripts/compare_variable_registry.py
```

2026-11-26 기준 전체 55개 테스트가 통과했다(기존 48개 포함). 계약 시험은
Core·HeatExchanger 합성 flat/composite DAE의 잔차·Jacobian·DASSL 적분,
포트 검증, 저마하수 폐루프 정상/과도 closure와 변수표 동기화를 확인한다.
합성 수치는 실제 HTGR 검증값이 아니다.

## 현재 제공 범위

| 부분 | 상태 | 사용 범위 |
| --- | --- | --- |
| 상태·포트·전역 DAE 계약 | 구현 | `StateRegistry`, `Component`, `SystemAssembly`, `StreamJunction`, `ThermalSurfacePort`, `ThermalLink` |
| 시간적분 | 구현 | DASSL, 일치 초기화, 정상상태 풀이, 데모 계산 |
| Core·HeatExchanger 상위 및 장치 패키지 | Composite 조립 구현 | 실제 fluid/solid/power/SG 하위/압축기 등 물리 구현 필요 |
| 헬륨 물성·열전달 상관식 | API만 제공 | 실제 provider와 적용범위 검증 필요 |
| HTGR 정상상태·LOFC | 실행 불가 | 형상·운전점·사고 이력·검증자료 필요 |

기본 [형상](config/geometry.yaml)·[운전](config/operation.yaml) 설정의
`null`은 미정 값이다. 실제 모델 실행에 시험용
`tests/fixtures/` provider를 자동 주입하지 않는다.

## Python에서 솔버 사용

`DASSL`은 지수 1 DAE를 적분하는 Python 솔버다. 1~5차 후진차분공식(BDF)과 뉴턴법을 사용하며, 가중 RMS 오차에 따라 스텝 간격과 차수를 조정한다.

```python
from examples.demo_system import make_demo_system
from solvers import DASSL, solve_steady_state

system = make_demo_system()
registry = system.registry
result = DASSL(
    dt=0.1, rtol=1e-6, atol=1e-8,
    residual_scale=system.residual_scales(),
).solve(
    system.residual,
    t_span=(0.0, 5.0),
    y0=registry.initial_state(),
    differential_mask=registry.differential_mask(),
)

y_equilibrium, ydot_equilibrium = solve_steady_state(
    system.residual, 0.0, registry.initial_state()
)
```

`differential_mask`를 주면 미분변수의 초깃값을 유지하면서 잔차식을 만족하는 대수변수와 미분변수의 초기 변화율을 구한다. 잔차식이 대수변수의 시간미분에 의존하지 않을 때 사용할 수 있다. 일치하는 초기조건을 알고 있다면 마스크를 생략하고 `y0`와 `ydot0`를 함께 전달한다. `solve_steady_state`는 모든 변화율을 0으로 두고 평형 상태를 구한다.

`result.t`, `result.y`, `result.ydot`에는 초기점을 포함한 시각, 상태, 변화율이 담긴다. `result.orders[i]`, `result.error_norms[i]`, `result.newton_iterations[i]`는 `t[i]`에서 `t[i+1]`까지 수용된 스텝의 정보다. 해석적 야코비안을 사용할 때는 `∂F/∂y + cj·∂F/∂ydot`를 반환하는 `jacobian(t, y, ydot, cj)`를 `solve`에 전달한다. 기본값은 유한차분이다.

## 모델 개발과 문서

담당 모델은 [공통 계약 v3](docs/interface_spec.md)에 따라 독립
미지수만 Registry에 등록하고, 같은 trial의 상태에서 계산한 1차원
잔차와 고정된 잔차 ID·scale을 제공한다. Core와 HeatExchanger의 실제
하위 모델·coupling 구현은 생성자에 명시적으로 주입한다. 접속 포트는
모든 구성요소를 등록한 뒤 검증한다. 각 장치의 실제 입력·방정식이
없으면 명시적 예외가 발생한다.

| 문서 | 내용 |
| --- | --- |
| [인터페이스 계약 v3](docs/interface_spec.md) | 상태·별칭, 포트 부호, provider API, Core 열률 매핑, 폐루프 자유도 |
| [설정 스키마](docs/config_schema.md) | YAML 키, 단위·shape, 조건별 필수 입력과 `null`의 의미 |
| [아키텍처](docs/architecture.md) | `physics/`·`models/`·`interfaces/`·`system/`·`solvers/`의 책임 |
| [구현 현황](docs/implementation-status.md) | 실행 가능한 범위와 테스트 증거 |
| [통합 계획](docs/integration_plan.md) | 역할별 소유권과 G0~G3 연결 순서 |
| [변수표](docs/variable_registry.csv)·[방정식표](docs/equation_registry.csv) | 모델 변수·식의 master 목록 |
| [변경안 제출법](docs/registry_patches/README.md)·[Role 5 인수인계](docs/handoffs/role5.md) | 담당자별 patch와 시험·미정 항목 |
| [물리 가정](docs/assumptions.md) | 문헌 기준값의 차이와 아직 정하지 않은 모델 선택 |
| [담당자별 지시서](docs/instructions/) | Role 1~5의 독립 개발 범위와 현재 코드 API |

master 변수표는 기존 361개 항목, 공통 계약 후보 11개,
기준설계 데이터표[1] 대조 신규 88개를 합친 **460개**이며 [Excel 변수표](docs/HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx)와
12개 대응 필드의 의미상 불일치를 0개로 맞췄다. CSV가 코드 개발용
단일 기준이고 Role 5가 patch를 검토해 Excel을 동기화한다. 기존 62셀
차이와 근거는 [변경 이력](docs/registry_patches/step01c_variable_sync_history.json)에
보존했다. 새 기준값은 출처의 설계·검산 조건으로 기록했으며 실제 초기값은 추정해 채우지 않았다.

## MHTGR-350 SG·RCCS 조사 입력

팀의 「MHTGR-350 SG/RCCS 보완본[4]」(2026-11, pp. 4–24)을 대조해 원래 통합안 288개에 SG/RCCS 관련 변수 73개를 추가했다. 아래 수치는 설계·검산 자료이며 실행 가능한 모델 입력 세트를 뜻하지 않는다. 현재 [열교환기](models/heat_exchanger/)의 상위 Composite만 구현됐고 하위 물리 모델과 [압력용기](models/vessel/)에는 실제 변수 등록과 잔차식이 필요하다.

원 설계서 DOE-HTGR-86-129 Rev.1/HFD-42102 Rev.2에 따른 SG 기준은 연속 전열관 439개, 외경 22.225 mm, 벽두께 3.302 mm, 외면 기준 총 전열면적 약 4026.539 m²다. 관벽 모델은 He 쪽과 물·증기 쪽 반벽의 전도와 축열을 다룬다. He 관군 마찰차압 9.651 kPa와 SG 전체 양단 차압 25.5 kPa는 측정 경계가 다르다. SG-2의 관내 비등·2상 차압 상관식은 보완본의 예비 기준이며 실제 고압 나선관 검증이 필요하다.

RCCS 상승관은 227개이고 각 관의 내형은 40.48 × 240.48 mm다. 보완본의 총 내부 유로면적 2.20976 m²와 수력직경 0.06930 m를 기록했다. 전체 자연순환 유량 계산에는 상승관 외에도 하강통로, 루버, 플레넘, 배기부의 손실과 외기 경계가 필요하다. 설계값·기하 계산값·기존 검산값의 차이와 미정 항목은 [물리 가정](docs/assumptions.md)에 기록했다.
