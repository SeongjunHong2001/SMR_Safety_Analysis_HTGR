# Role 5 통합 계획 (계약 v3)

기준 `main` SHA: `d035aaed1f99de838e3b0f2498b34777cbed5ed7`.
작업 브랜치: `feat/role5-contract-v3`. 이 기준에서 각 담당자가
자기 브랜치를 만들고 `docs/interface_spec.md` 및 계약 테스트를 따른다.

| 게이트 | 요구 증거 | 현재 상태 |
| --- | --- | --- |
| G0 공통 골격 | 기존 데모/pytest, 공개 import, Registry·연결·scale·ID, 합성 flat/composite DAE, 합성 폐루프 rank | 계약 코드·fixture에 구현 |
| G1 담당 모델 단독 | Role 1 입력 출처, Role 2 물성/유동, Role 3 노심 고체/열, Role 4 SG/RPV/RCCS 각각의 단독 보존·범위 시험 | 미구현 |
| G2 실제 접속 | 모두 같은 실제 입력과 EOS를 사용한 Core→SG→압축기 폐루프, RPV→RCCS 열경로, 정칙성·일치 초기화 | 미구현 |
| G3 검증 | 문헌 운전점·LOFC 이력에 대한 정량 비교와 불확도/범위 | 자료 미확보 |

## 역할과 접속 순서

Role 1은 실제 `config/geometry.yaml`, `operation.yaml`, `scenarios/`,
`data/validation/`의 출처를 확정한다. Role 2는 헬륨 물성, 유동·압력손실,
Core fluid, piping/plenum/compressor를 맡는다. Role 3은 Core solid/power/
coupling과 공통 열전달 상관식을 맡는다. Role 4는 SG, vessel, RCCS와
SG/RCCS 전용 상관식을 맡는다. Role 5는 `interfaces/`, `system/`,
`solvers/`, `main.py`, Core 상위, ThermalLink, boundary, 공개 export,
공유 문서와 master 표를 통합한다. 물리 담당자는 자기 변수/식 변경안을
`docs/registry_patches/roleN/`에 제출하고 master는 Role 5가 병합한다.

1. 담당 모델은 `ports()`, `residual_ids()`, `residual_scales()`와
   Register의 독립 상태 목록을 먼저 확정한다. 미정 데이터는 `null`로 둔다.
2. fixture로 모델 단독의 shape, 유량 역전/정체, 물성 적용범위,
   열·질량 보존을 시험한다. fixture는 생산 builder에서 사용하지 않는다.
3. Role 2/5가 실제 p_sys, 국부압 기준점, 정상 재고량 제약의 자유도
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
