# 개발 구현 현황

2026년 11월 25일 공통 계약 v3의 G0 골격이 [PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)로 `main`에 병합됐다(병합 SHA `52016a8`). STEP 01B에서 읽기 전용 trial과 열경계 연결을, STEP 01C에서 HeatExchanger 상위 Composite·변수표 동기화·저마하수 폐루프 수학 계약을 보완했다. 실제 MHTGR-350 정상상태 또는
LOFC 해석을 실행한 결과는 없다. 현재 실행 가능한 수치는 두 열용량
데모와 `tests/fixtures/`의 명시적 합성 모델에 한정된다.

| 기능 | 상태 | 증거·남은 일 |
| --- | --- | --- |
| [Component/Registry 계약](../interfaces/) | G0 구현 | 독립 상태·별칭·owner·shape·단위, 공개 name/key/ports/잔차 ID·scale, 읽기 전용 context와 HeliumProperties Protocol |
| [접속·전역 조립](../system/) | G0+01C 구현 | N포트 2N 잔차, 읽기 전용 trial, 열링크 전달, 조립 순서·ID·scale, 전역 포트 중복 검사와 결합 Jacobian 사전 검증. `HeliumInventory` 합산 계약 구현. 실제 압력손실과 HTGR 폐루프 물리 검증 대기 |
| [DASSL](../solvers/) | 기존 구현 확장 | 일치 초기화와 가변 BDF 1~5차. 비영 패턴을 지정하면 색칠 유한차분·희소 선형 풀이; 기본 밀집 경로 유지. 합성 flat/composite DAE 적분 결과 일치 |
| [Core](../models/core/) | 상위 조립 계약 | fluid/solid/power/coupling 생성자 주입과 현재 trial의 보존적 열률 전달. 실제 하위 물리식은 Roles 2/3 대기 |
| [HeatExchanger](../models/heat_exchanger/) | 상위 Composite 구현 | helium_side→wall→secondary_side DI 조립·ID·scale·외부 포트·보존 열률·출력, 합성 flat/Composite·DASSL 검증. 실제 하위 SG 물리식과 상관식은 Role 4 대기 |
| [압축기, 배관, 플레넘, 용기, RCCS](../models/) | 패키지·대표 클래스 | 실제 register/residual은 명시적 미구현. 압축기 provider는 계약만 존재 |
| [공통 열전달](../physics/heat_transfer/) | API 동결 | base 함수 서명과 Role 3/4 파일 경계. 실제 식은 명시적 미구현 |
| [외부 열경계](../models/thermal_link.py) | G0+01B 구현 | `ThermalSurfacePort`의 스칼라/배열 셀 온도, `ThermalLink`의 보존 열률을 `SimulationContext.heat_rates`에 trial별로 전달. 실제 Core/RPV/RCCS 열물리는 미구현 |
| [변수표](variable_registry.csv)·[Excel](HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx) | 01C 동기화 + 기준설계 대조 | CSV 단일 기준, 460개 고유 항목·12개 대응 필드 의미상 불일치 0개. 2026-11-26 기준설계 데이터표[1]에서 88개를 추가하고 기존 40개를 보완. STEP 01C 당시 372개·62셀 차이 근거는 [이력](registry_patches/step01c_variable_sync_history.json)에 보존. 미확정 초기 물리값은 미정 유지 |
| [방정식표](equation_registry.csv) | 저마하수 수학 계약 | `p_sys` 소유·국부압 기준·EOS·셀 질량/에너지·정상 재고량 교체·과도 구성 명시. 합성 7미지수 정상/과도 rank·일치 초기화·DASSL 검증, 실제 HTGR 검증 대기 |
| [실제 설정](../config/)·[검증자료](../data/) | 미정 | 형상, 물성, 운전점, 사고 이력·벤치마크 자료가 필요. `null`은 0이 아님 |

## 검증 범위

`python -m pytest -q`: **65개 통과**(2026-11-26 확인; 기존 55개 회귀 통과).
`python scripts/compare_variable_registry.py`: CSV/Excel 각각 460개,
12개 대응 필드 의미상 불일치 0개.
`python main.py --demo`: 기본 0–5초 DAE 데모 완료, 최종 상태
`[316.19074506, 306.64335141, 47.73696821]` (K,K,W).
합성 flat/composite 잔차·Jacobian·적분 일치, derived 열률의 온도 궤적,
매질/포트 검증, 외부 열링크의 trial별 전달, HeatExchanger 합성 조립·DASSL,
저마하수 합성 폐루프의 정상/과도 rank·일치 초기화·질량/에너지 수지를 시험했다.
이 결과는 실제 HTGR 물리 정확도 증거가 아니다.

실제 G1/G2/G3 조건과 담당자별 접속 순서는 [통합 계획](integration_plan.md)에,
물리적 미정과 문헌값의 차이는 [가정](assumptions.md)에 기록한다.
