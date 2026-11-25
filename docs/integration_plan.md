# Role 5 통합 계획 (계약 v3)

계약 v3는 기준 `d035aaed`에서 작성해 [PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)로 `main`에 병합됐다(병합 SHA `52016a8`). 이후 담당자는 작업 시작 시 최신 `main`을 가져와 자기 브랜치를 만들고 [인터페이스 계약](interface_spec.md)과 `tests/contracts/`를 따른다.

| 게이트 | 요구 증거 | 현재 상태 |
| --- | --- | --- |
| G0 공통 골격 | 기존 데모/pytest, 공개 import, Registry·연결·scale·ID, Core/HeatExchanger 상위 Composite, 합성 flat/composite DAE와 DASSL, 합성 저마하수 폐루프 정상/과도 rank·질량/에너지. STEP 01B의 읽기 전용 trial·열링크 전달과 STEP 01C의 미해결 계약 보완 | PR #5 병합, 01B/01C 로컬 보완 |
| G1 담당 모델 단독 | Role 1 입력 출처, Role 2 물성/유동, Role 3 노심 고체/열, Role 4 SG/RPV/RCCS 각각의 단독 보존·범위 시험 | 미구현 |
| G2 실제 접속 | 모두 같은 실제 입력과 EOS를 사용한 Core→SG→압축기 폐루프, RPV→RCCS 열경로, 정칙성·일치 초기화 | 미구현 |
| G3 정상상태 | 실제 운전점의 일치 초기화·보존 검산과 같은 경계의 문헌값 비교 | 자료·물리 모델 미확보 |
| G4 과도·검증 | 실제 부분 유량 감소/LOFC 이력과 같은 조건의 비교 및 적용 범위 | 자료·물리 모델 미확보 |

## 역할과 접속 순서

Role 1은 실제 `config/geometry.yaml`, `operation.yaml`, `scenarios/`,
`data/validation/`의 출처를 확정한다. Role 2는 헬륨 물성, 유동·압력손실,
Core fluid, piping/plenum/compressor를 맡는다. Role 3은 Core solid/power/
coupling과 공통 열전달 상관식을 맡는다. Role 4는 SG, vessel, RCCS와
SG/RCCS 전용 상관식을 맡는다. Role 5는 `interfaces/`, `system/`,
`solvers/`, `main.py`, Core·HeatExchanger 상위, ThermalLink, boundary, 공개 export,
공유 문서와 master 표를 통합한다. 물리 담당자는 자기 변수/식 변경안을
`docs/registry_patches/roleN/`에 제출하고 master는 Role 5가 병합한다.
CSV가 코드 개발용 기준이며 Excel은 동기화된 배포본이다. 병합마다
`python scripts/compare_variable_registry.py`로 모든 대응 필드를 확인한다.

1. 담당 모델은 `ports()`, `residual_ids()`, `residual_scales()`와
   Register의 독립 상태 목록을 먼저 확정한다. 미정 데이터는 `null`로 둔다.
2. fixture로 모델 단독의 shape, 유량 역전/정체, 물성 적용범위,
   열·질량 보존을 시험한다. fixture는 생산 builder에서 사용하지 않는다.
3. Role 2/5가 확정된 `p_sys` 단일 소유·국부압 기준식·정상 재고량 제약의 자유도
   ledger를 작성한다. 종속 질량식을 중복 등록하지 않는다.
4. 각 모델이 같은 trial을 읽도록 Core/SG 열률 매핑과 접속을 묶는다.
   동일 열경계는 한 계산 경로에서 양쪽에 반대 부호로 전달한다.
5. 전역 `n_unknown=n_residual`, scaled Jacobian rank/condition,
   index-1 구조, 일치 정상/과도 초기조건을 검증한다. 이후 DASSL로
   실제 해석을 수행한다.

담당 모델 PR은 최신 `main`을 반영하고 자기 테스트와 기존 pytest,
`python main.py --demo`를 실행한다. 공통 계약 변경은 Role 5의 작은
별도 PR로 검토한다. 실제 모델이 아직 미완성이라도 독립 단독 PR은
제출할 수 있다. 실제 접속 G2는 필요한 모든 실제 모델이 구현돼야 한다.
