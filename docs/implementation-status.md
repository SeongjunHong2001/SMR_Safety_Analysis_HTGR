# 개발 구현 현황

2026년 11월 25일 공통 계약 v3의 G0 골격이 [PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)로 `main`에 병합됐다(병합 SHA `52016a8`). STEP 01B에서 읽기 전용 trial과 열경계 연결을, STEP 01C에서 HeatExchanger 상위 Composite·변수표 동기화·저마하수 폐루프 수학 계약을 보완했다.

2026-11-27 기준으로 Role 1–4 초기 모델과 Role 5 builder·실행기가 들어와 **S0 정상상태 두 조건(core_benchmark, system_loop)을 config에서 조립해 푼다.** 다만 일부 입력은 원문 값이 아니라 [가정]이고(`config/assumptions/provisional_s0.yaml`), 압력손실 K 보정 전이며, LOFC(G4)는 시작하지 않았다.

## 게이트 (Role 5 지시서 STEP D-4)

| 게이트 | 상태 | 증거 | 차단 항목 (소유자) |
| --- | --- | --- | --- |
| G0 공통 계약 | pass | PR #5, `tests/contracts` | — |
| G1 역할별 단독 | pass (초기 모델) | `tests/physics` 58개, `tests/components` 34개: 물성·압력손실·FVM 커널, 노심 유체, 플레넘·배관·순환기, 노심 고체·출력·결합, RPV·RCCS-0, SG-0 단독 보존·범위 시험 | MLSV 밸브, 순환기 성능곡선·관성 M09 (2), 실제 출력분포·점동특성·붕괴열 (3), SG-1/2·RCCS-1·헤드 HS-09/10 (4) |
| G2 실제 접속 | pass (가정 입력 2개 포함) | `system/builder.py`가 연결부 표로 조립: system_loop 481 미지수 = 잔차 481, p_sys 단일 소유(하부 플레넘), 정상 폐합 1개(`lower_plenum.mass` → p_sys − p0), 유체 포트 전부 연결, Core→RPV→RCCS-0 열링크 3개. 희소 Jacobian 구조 rank 정칙(`validate_for_solve`), 일치 초기화 성공 | 가정 입력: SG 용기 상부 영역을 110·350에 나누는 비율(M10 남은 부분) (1). 우회 D_h·SG He 유로(M11)·300·330·340 체적은 R19·R1 치수로 채움(#31) |
| G3 정상상태 | partial | 정상해·과도 잔차 재확인 ≤ 1.04e-7(scaled), 무외란 DASSL 10 s 표류 3.6e-8, 전역 질량 오차 4e-16, 에너지 오차 7e-10(상대). 열수지: SG 352.67 MW = 350 + 3.1 − 0.43 MW. R1 노심 출구 960.5 K(참고 960.15 K), RPV 최고 463 K < 643.15 K(판정 통과). K 보정(`config/calibration/s0_pressure_loss_K.yaml`) 후 PSID 구간 설계손실 9개 일치(1 Pa 이내), 순환기 승압 91.93 kPa(PSID 91.0), 지점 압력 ±3.5 kPa. He 재고량 2,558.5 kg(R22 2,585, −1.0 %). core_benchmark는 가정 없이 노심 Δp 20.7 kPa(R3 19.95) | 110·350 분배 가정 (1), Ex. 2 출력분포 파일 없음(균일 가정) (1/3), pseudo-transient fallback 미구현 (5) |
| G4 과도/LOFC | not started | 출력 +10 % 램프 3 s 과도에서 질량 보존 2.5e-10, 에너지 수지 끝값 8e-5 (합성 아닌 실제 모델 과도의 보존 시험일 뿐 LOFC 아님) | `system/scenario.py`(STEP C), loss_of_flow.yaml, 순환기 감속·붕괴열 입력 (1/2/3/5) |

## 기능별 상태

| 기능 | 상태 | 증거·남은 일 |
| --- | --- | --- |
| [Component/Registry 계약](../interfaces/) | G0 구현 | 독립 상태·별칭·owner·shape·단위, 공개 name/key/ports/잔차 ID·scale, 읽기 전용 context와 HeliumProperties Protocol |
| [접속·전역 조립](../system/assembly.py) | G0+01C 구현 | N포트 2N 잔차, 읽기 전용 trial, 열링크 전달, 조립 순서·ID·scale, 전역 포트 중복 검사와 결합 Jacobian 사전 검증. `HeliumInventory` 합산 |
| [builder](../system/builder.py) | STEP B 구현 | config(geometry·operation·scenario) + 선택적 [가정] 파일 → SystemAssembly. 노드 ID→장치, `geometry.junctions`→StreamJunction(병렬 묶음·내부 생략), 개방 영역 경계, 폐루프 기준압 폐합, 열링크. 미정 입력은 파일·키·시나리오·조건을 모아 한 번에 실패(`MissingInputsError`). 가정 파일은 null만 채우고 원문 값 덮어쓰기는 오류 |
| [실행기](../system/simulation.py)·[main](../main.py) | STEP B 구현 | 정상해(희소 FD Newton) → 과도 잔차 재확인 → 일치 초기화 → DASSL 무외란 유지 → 보존 진단 → CSV(Registry 이름/shape 평탄화, 중복 검사)·메타 JSON(SI 단위, 시나리오, 가정, 코드 버전, 입력 해시, 솔버 통계). `python main.py --scenario steady_state [--condition ...] [--assumptions ...]` |
| [압력손실 K 보정](../system/calibration.py) | 구현 | PSID 구간 설계손실(비가역 = 마찰 + 형상)에 맞춰 구간별 연결부 K_form 증분을 식 7로 반복 보정, 보정 파일(`meta.kind: calibration`)은 기준 조건·가정 파일과 함께만 쓴다. `python main.py --calibrate` |
| [전역 보존 진단](../validation/conservation.py) | STEP B 구현 | `M_he_loop`, `U_he_loop`(Σ M h − p_sys V), `U_solid_total`(Σ ∫C(T)dT), 경계 열률(출력·순환기·SG·RCCS·개방 경계 엔탈피), 정상 수지와 과도 이력 수지 |
| [DASSL](../solvers/) | 기존 구현 | 일치 초기화와 가변 BDF 1~5차, 희소 색칠 FD. `consistent_initialize`에 패턴을 줄 때는 열 순서가 [대수 y, 미분 ẏ]임에 주의 |
| [Core](../models/core/) | Role 2·3 초기 모델 | 병렬 3유로 유체, 단위격자 고체(링 포함), 지정 출력, McEligot 결합 |
| [HeatExchanger](../models/heat_exchanger/) | SG-0 | He 측 9셀, 축열 없는 벽, 일정 온도 물 열침원(U 보정/고정). SG-1/2 미구현 |
| [압축기·배관·플레넘](../models/) | Role 2 초기 모델 | 순환기 유량 지정·성능곡선(곡선 미확보), 1D 배관, 0D 플레넘. 밸브 모델 없음(MLSV는 열린 플레넘 [가정]) |
| [용기·RCCS](../models/) | Role 4 초기 모델 | RPV 측벽 축방향 셀, RCCS-0 고정 온도 경계. 헤드 HS-09/10, RCCS-1 미구현 |
| [공통 열전달](../physics/heat_transfer/) | Role 3 구현 | 전도·대류·복사 기본 함수, McEligot·자연대류, R1 고체 물성(Ex. 2b) |
| [변수표](variable_registry.csv)·[Excel](HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx) | 동기화 | 460개, 12개 대응 필드 의미상 불일치 0개. 역할별 patch는 `docs/registry_patches/role1–4` (마스터 병합 대기) |
| [실제 설정](../config/)·[검증자료](../data/) | Role 1 입력 | 형상·운전점·S0 조건 확보. M11·우회 D_h·하부 블록 셀별 D_h·M10 대부분 보완(#31). 미확보: M04, M07, M09, M10 상부 영역 분배 |

## 검증 범위

`python -m pytest -q`: **218개 통과** (2026-11-27, 약 120 s).
`python scripts/compare_variable_registry.py`: CSV/Excel 각각 460개, 12개 대응 필드 의미상 불일치 0개.
`python main.py --demo`: 기본 0–5초 DAE 데모 완료, 최종 상태 `[316.19074506, 306.64335141, 47.73696821]` (K,K,W).
`python main.py --scenario steady_state --condition system_loop`: 가정 파일 없이는 미정 입력 8개를 적고 종료 코드 2.
`--assumptions config/assumptions/provisional_s0.yaml`을 주면 정상해·검증·유지 적분·출력까지 완료.

| S0 결과 (2026-11-27 #31 형상 보완 후) | core_benchmark | system_loop | 참고값 (판정 아님) |
| --- | --- | --- | --- |
| 입력 | config만 (가정 없음) | config + 가정 2개 + 보정 K 10개 | — |
| 미지수 | 347 | 481 (정상) / 480 (과도) | — |
| 노심 출구 혼합 | 960.5 K | 960.3 K | R1 960.15 K |
| 연료 최고 (컴팩트 중심) | 1124.7 K | 1119.2 K | R3 1261 K (다른 출력분포) |
| 노심 Δp (150→230) | 20.7 kPa | 34.7 kPa (K 보정) | R3 19.95 kPa / PSID 34.5 kPa |
| RPV 최고 | 462.6 K | 462.5 K | 판정 643.15 K 미만 → 통과 |
| RCCS 제거열 | 426 kW | 425 kW | — |
| SG 제거열 / U | — | 352.67 MW / 474 W/(m²·K) | PSID 352.2 MW |
| 순환기 승압 | — | 91.93 kPa (K 보정 전 약 37 kPa) | PSID 91 kPa |
| He 재고량 | 809 kg | 2,314.9 kg (+ 정체 243.6 = 2,558.5) | R22 2,585 kg |

이 결과는 초기 모델과 가정 입력의 정상해다. 문헌값과의 정량 비교는 같은 경계·운전점일 때만 의미가 있다.
실제 G1/G2/G3 조건과 담당자별 접속 순서는 [통합 계획](integration_plan.md)에, 물리적 미정과 문헌값의 차이는 [가정](assumptions.md)에 기록한다.
