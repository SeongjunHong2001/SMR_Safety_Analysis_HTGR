# SMR Safety Analysis · HTGR

MHTGR-350의 정상상태와 가압 냉각재 유량 상실(P-LOFC/PCC)을
해석하기 위한 Python 프로젝트다. 구성요소의 상태와 보존식을 하나의
미분대수방정식(DAE) `F(t, y, ydot) = 0`으로 조립한다.

**현재 S0 정격 정상상태를 실제 Role 1–4 초기 모델로 조립해 풀 수 있다.**
- **노심 단독:** OECD/NEA 벤치마크 조건(`core_benchmark`)
- **1차계통 폐루프:** 순환기·덕트·SG-0 열침원 포함(`system_loop`)

사고 과도(LOFC)는 아직 실행할 수 없다. 시나리오 처리, 점동특성·붕괴열, 순환기 감속 모델이 없다.
기능별 범위와 시험 증거는 [구현 현황](docs/implementation-status.md)에 있다.

## 설치와 첫 실행

Python 3.10 이상이 필요하다(3.13에서 시험함). 저장소를 내려받은 뒤 프로젝트 최상위 폴더에서 실행한다.

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

## S0 정상상태 실행

```bash
python main.py --scenario steady_state                       # 노심 단독 벤치마크 (기본 조건 core_benchmark)
python main.py --scenario steady_state --condition system_loop \
    --calibration config/calibration/s0_pressure_loss_K.yaml # 1차계통 폐루프 (PSID 압력손실 K 보정 적용)
```

config(`geometry.yaml`, `operation.yaml`, `scenarios/steady_state.yaml`)만 읽어 실제 모델을 조립한다
([`system/builder.py`](system/builder.py)). 실행 순서는 다음과 같다 ([`system/simulation.py`](system/simulation.py)).

1. 정상해를 구한다.
2. 그 해를 과도 잔차식에 다시 넣어 확인한다.
3. 일치 초기화 후 외란 없이 10 s 적분해 정상이 유지되는지 본다.
4. 전역 질량·에너지 수지를 진단한다.
5. 결과 파일을 쓴다.

한 번에 약 10–20 s 걸린다.

| 결과 파일 (`output/results/`) | 내용 |
| --- | --- |
| `steady_state_<조건>_states.csv` | `time_s`와 모든 상태 변수 (변수표 이름, 배열은 `이름[i]`) |
| `steady_state_<조건>_outputs.csv` | 장치별 출력 (온도·압력·열률 등) |
| `steady_state_<조건>_meta.json` | 단위, 입력 파일·해시, 가정·보정값, 코드 버전, 솔버 통계, 보존 진단, 참고값 비교·판정 |

| 옵션 | 뜻 |
| --- | --- |
| `--condition 이름` | `conditions` 중 하나 (기본값은 시나리오의 `default_condition`) |
| `--calibration 파일` | 보정 결과(연결부 K) 덮어쓰기. 보정한 조건에만 쓸 수 있다 |
| `--calibrate` | `system_loop` 정상해로 K를 다시 보정해 `--calibration-out`(기본 `config/calibration/s0_pressure_loss_K.yaml`)에 쓴다 |
| `--assumptions 파일` | config의 `null` 입력을 [가정]값으로 채운다 ([형식](config/assumptions/README.md)) |
| `--no-hold` | 무외란 유지 적분을 건너뛴다 |
| `--config-dir`, `--output-dir` | 입력·출력 폴더 바꾸기 |

config에 필요한 값이 `null`이면 계산하지 않고, 빠진 파일·키·조건 목록을 출력한 뒤 종료 코드 2로 끝난다.
아직 지원하지 않는 시나리오(`loss_of_flow`)도 같은 방식으로 멈춘다.

| S0 주요 결과 (초기 모델) | core_benchmark | system_loop | 참고값 (판정 아님) |
| --- | --- | --- | --- |
| 노심 출구 혼합 온도 | 960.5 K | 960.3 K | R1 960.15 K |
| 노심 압력강하 | 20.7 kPa | 34.7 kPa (K 보정) | R3 19.95 / PSID 34.5 kPa |
| RPV 최고 온도 | 462.6 K | 462.5 K | 판정 기준 643.15 K 미만 → 통과 |
| SG 제거열 / 순환기 승압 | — | 352.7 MW / 91.9 kPa | PSID 352.2 MW / 91 kPa |

결과는 균일 출력분포 [가정], R1 고정 물성, 일부 도면해석 형상 위의 초기 모델 값이다.
근거와 한계는 [구현 현황](docs/implementation-status.md)과 [Role 5 인수인계](docs/handoffs/role5.md)에 있다.

## 시험

```bash
python -m pytest -q
python -m pytest -q tests/contracts
python scripts/compare_variable_registry.py
```

2026-11-27 기준 전체 218개 테스트가 통과한다(약 2분).

| 폴더 | 시험 대상 |
| --- | --- |
| `tests/contracts` | 공통 계약: 합성 DAE의 잔차·Jacobian·DASSL, 포트 검증, 희소 Jacobian, 저마하수 폐루프, 변수표 동기화 |
| `tests/inputs` | Role 1 입력 일관성 |
| `tests/physics` | 물성·FVM 커널·압력손실·열전달 |
| `tests/components` | 장치 단독: 노심 유체·고체, 플레넘·배관·순환기, RPV·RCCS-0, SG-0 |
| `tests/integration` | builder, S0 두 조건 실행, K 보정, 데모 |

## 현재 제공 범위

| 부분 | 상태 | 사용 범위 |
| --- | --- | --- |
| 상태·포트·전역 DAE 계약 | 구현 | `StateRegistry`, `Component`, `SystemAssembly`, `StreamJunction`, `ThermalSurfacePort`, `ThermalLink` |
| 시간적분 | 구현 | DASSL, 일치 초기화, 정상상태 풀이, 희소 Jacobian |
| 헬륨 물성·열전달 (Role 2·3) | 초기 모델 | R1 헬륨 물성, 1D FVM·압력손실, McEligot 대류, 전도·복사, R1 고체 물성(Ex. 2b) |
| 노심 (Role 2·3) | 초기 모델 | 병렬 3유로 유체(평균·고온·우회), 단위격자 고체, 지정 출력(균일 [가정]) |
| 플레넘·배관·순환기 (Role 2) | 초기 모델 | 0D 플레넘, 1D 배관, 유량 지정 순환기. MLSV 밸브·순환기 감속 없음 |
| RPV·RCCS·SG (Role 4) | 초기 모델 | RPV 측벽, RCCS-0 고정 온도 경계, SG-0 일정 온도 물 열침원. SG-1/2, RCCS-1, 헤드 열구조체 없음 |
| builder·실행기·K 보정 (Role 5) | 구현 | config → 계통 조립, 정상해·검증·유지 적분·출력, PSID 설계손실 K 보정 |
| HTGR 정상상태 (S0) | 실행 가능 | `core_benchmark`, `system_loop` |
| 사고 과도 (LOFC) | 실행 불가 | 시나리오 처리, 점동특성·붕괴열, 순환기 성능곡선·관성(M09) 필요 |

[형상](config/geometry.yaml)·[운전](config/operation.yaml) 설정의 `null`은 미정 값이다(0이 아니다).
값마다 근거를 [원문]·[계산]·[도면해석]·[가정]·[제안]으로 표시했다.
builder는 `null`을 기본값으로 채우지 않고, 시험용 `tests/fixtures/` provider를 실제 실행에 주입하지 않는다.

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
    preflight=system.validate_for_solve,
)

y_equilibrium, ydot_equilibrium = solve_steady_state(
    system.residual, 0.0, registry.initial_state()
)
```

`differential_mask`를 주면 미분변수의 초깃값을 유지하면서 잔차식을 만족하는 대수변수와 미분변수의 초기 변화율을 구한다. 잔차식이 대수변수의 시간미분에 의존하지 않을 때 사용할 수 있다. 일치하는 초기조건을 알고 있다면 마스크를 생략하고 `y0`와 `ydot0`를 함께 전달한다. `solve_steady_state`는 모든 변화율을 0으로 두고 평형 상태를 구한다.

`result.t`, `result.y`, `result.ydot`에는 초기점을 포함한 시각, 상태, 변화율이 담긴다. `result.orders[i]`, `result.error_norms[i]`, `result.newton_iterations[i]`는 `t[i]`에서 `t[i+1]`까지 수용된 스텝의 정보다. `preflight=system.validate_for_solve`는 일치 초기화 뒤 첫 스텝의 결합 Jacobian 정칙성을 확인한다. 해석적 야코비안을 사용할 때는 `∂F/∂y + cj·∂F/∂ydot`를 반환하는 `jacobian(t, y, ydot, cj)`를 `solve`에 전달한다. 기본값은 밀집 유한차분이며, `jacobian_sparsity`에 `∂F/∂y`와 `∂F/∂ydot`의 가능한 비영 성분을 모두 포함한 패턴을 주면 묶음 유한차분과 SciPy 희소 선형 풀이를 사용한다. 패턴을 누락하면 Jacobian이 잘못 계산되므로 모델의 결합 구조에서 검증해야 한다.

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
| [변경안 제출법](docs/registry_patches/README.md)·[담당자별 인수인계](docs/handoffs/) | 담당자별 patch, 구현 범위, 시험, 미정 항목 (role1–5) |
| [가정 입력 형식](config/assumptions/README.md)·[보정 결과](config/calibration/) | `null` 입력을 채우는 [가정] 파일 규칙, 연결부 K 보정값 |
| [물리 가정](docs/assumptions.md) | 문헌 기준값의 차이와 아직 정하지 않은 모델 선택 |
| [담당자별 지시서](docs/instructions/) | Role 1~5의 독립 개발 범위와 현재 코드 API |

master 변수표는 기존 361개 항목, 공통 계약 후보 11개,
기준설계 데이터표[1] 대조 신규 88개를 합친 **460개**이며 [Excel 변수표](docs/HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx)와
12개 대응 필드의 의미상 불일치를 0개로 맞췄다. CSV가 코드 개발용
단일 기준이고 Role 5가 patch를 검토해 Excel을 동기화한다. 기존 62셀
차이와 근거는 [변경 이력](docs/registry_patches/step01c_variable_sync_history.json)에
보존했다. 새 기준값은 출처의 설계·검산 조건으로 기록했으며 실제 초기값은 추정해 채우지 않았다.

## MHTGR-350 SG·RCCS 조사 입력

팀의 「MHTGR-350 SG/RCCS 보완본[4]」(2026-11, pp. 4–24)을 대조해 원래 통합안 288개에 SG/RCCS 관련 변수 73개를 추가했다. 아래 수치는 설계·검산 자료이며 실행 가능한 모델 입력 세트를 뜻하지 않는다.

실행 입력(`config/geometry.yaml`)의 SG 형상은 PSID(DOE-HTGR-86-024) Amd.3–6 기준이다 ([제안] 2026-11-26, 팀 합의 대기 — [Role 1 인수인계](docs/handoffs/role1.md)).
- **관:** 350개, 평균 유효 길이 164.5 m, 외면 전열면적 4020 m²
- **헬륨 측:** 최소 자유유로면적·수력직경을 Table 5.3-8 치수로 계산했다.

아래 보완본[4] Rev.1(439관) 값은 판본 비교용이다.
SG는 [SG-0 일정 온도 물 열침원](models/heat_exchanger/), RPV 측벽은 [압력용기](models/vessel/) 모델로 들어가 있다.

원 설계서 DOE-HTGR-86-129 Rev.1/HFD-42102 Rev.2에 따른 SG 기준은 연속 전열관 439개, 외경 22.225 mm, 벽두께 3.302 mm, 외면 기준 총 전열면적 약 4026.539 m²다. 관벽 모델은 He 쪽과 물·증기 쪽 반벽의 전도와 축열을 다룬다. He 관군 마찰차압 9.651 kPa와 SG 전체 양단 차압 25.5 kPa는 측정 경계가 다르다. SG-2의 관내 비등·2상 차압 상관식은 보완본의 예비 기준이며 실제 고압 나선관 검증이 필요하다.

RCCS 상승관은 227개이고 각 관의 내형은 40.48 × 240.48 mm다. 보완본의 총 내부 유로면적 2.20976 m²와 수력직경 0.06930 m를 기록했다. 전체 자연순환 유량 계산에는 상승관 외에도 하강통로, 루버, 플레넘, 배기부의 손실과 외기 경계가 필요하다. 설계값·기하 계산값·기존 검산값의 차이와 미정 항목은 [물리 가정](docs/assumptions.md)에 기록했다.
