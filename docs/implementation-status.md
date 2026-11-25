# 개발 구현 현황

2026년 11월 25일, 기준 `main` SHA `d035aaed1f99de838e3b0f2498b34777cbed5ed7`
위에 공통 계약 v3의 G0 골격을 작성했다. 실제 MHTGR-350 정상상태 또는
LOFC 해석을 실행한 결과는 없다. 현재 실행 가능한 수치는 두 열용량
데모와 `tests/fixtures/`의 명시적 합성 모델에 한정된다.

| 기능 | 상태 | 증거·남은 일 |
| --- | --- | --- |
| [Component/Registry 계약](../interfaces/) | G0 구현 | 독립 상태·별칭·owner·shape·단위, 공개 name/key/ports/잔차 ID·scale, 읽기 전용 context와 HeliumProperties Protocol |
| [접속·전역 조립](../system/) | G0 구현 | N포트 2N 잔차, 매질/스칼라/키 검증, 조립 순서·ID·scale, scaled 결합 Jacobian rank/condition 진단. 실제 압력손실·closure는 미정 |
| [DASSL](../solvers/) | 기존 구현 유지 | 일치 초기화와 가변 BDF 1~5차. 합성 flat/composite DAE 적분 결과 일치 |
| [Core](../models/core/) | 상위 조립 계약 | fluid/solid/power/coupling 생성자 주입과 현재 trial의 보존적 열률 전달. 실제 하위 물리식은 Roles 2/3 대기 |
| [SG, 압축기, 배관, 플레넘, 용기, RCCS](../models/) | 패키지·대표 클래스 | 기존 생성자/import 유지. 실제 register/residual은 명시적 미구현. SG 하위·압축기 provider는 계약만 존재 |
| [공통 열전달](../physics/heat_transfer/) | API 동결 | base 함수 서명과 Role 3/4 파일 경계. 실제 식은 명시적 미구현 |
| [변수표](variable_registry.csv)·[Excel](HTGR_전체시스템_통합변수표_SG_RCCS_보완반영.xlsx) | 동기화 | 기존 361개 + 신규 11개 = 372개. 5개 기존 Residual 분류 보정. 새 수치는 미정 유지 |
| [방정식표](equation_registry.csv) | 계약 ledger | 합성 DAE 3식, stream 접속 3종, 실제 폐루프 진단·closure 후보. 실제 장치 식 수·scale은 역할별 모델 PR에서 채움 |
| [실제 설정](../config/)·[검증자료](../data/) | 미정 | 형상, 물성, 운전점, 사고 이력·벤치마크 자료가 필요. `null`은 0이 아님 |

## 검증 범위

`python -m pytest -q`: **43개 통과**. 기존 32개와 공통 계약 11개다.
`python main.py --demo`: 기본 0–5초 DAE 데모 완료, 최종 상태
`[316.19074506, 306.64335141, 47.73696821]` (K,K,W).
합성 flat/composite 잔차·Jacobian·적분 일치, derived 열률의 온도 궤적,
매질/포트 검증, 폐루프 fixture의 정상/과도 rank를 시험했다.
이 결과는 실제 HTGR 물리 정확도 증거가 아니다.

실제 G1/G2/G3 조건과 담당자별 접속 순서는 [통합 계획](integration_plan.md)에,
물리적 미정과 문헌값의 차이는 [가정](assumptions.md)에 기록한다.
