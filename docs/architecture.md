# 공통 해석 구조 (Role 5 계약 v3)

```text
config/scenario + 명시적으로 선택한 provider
                 ↓
Core · HeatExchanger · Compressor · Piping · Plenum · Vessel · RCCS
                 ↓ 현재 trial y, ydot
interfaces/StateRegistry + system/connections + ThermalLink
                 ↓ 전역 F(t,y,ydot)=0, residual IDs, s_R
solvers/DASSL       ↓
결과·보존 진단·검증
```

`physics/`에는 현재 입력에서 값을 계산하는 물성·유동·압력손실·공통
열전달 함수와 상관식만 둔다. Component 상태, 시나리오, 시간적분을
소유하지 않는다. `physics/heat_transfer/`가 공통 열률 API의 유일한
위치이고 `models/core/coupling.py` 및 SG 내부 coupling은 격자·면적
매핑만 맡는다. 현재 물리식은 미구현이며 기본값을 반환하지 않는다.

`models/`에는 장치 상태·기하·포트·잔차를 둔다. Core는 fluid/solid/
power/coupling을 생성자에서 주입받아 한 대표 Component로 조립한다.
HeatExchanger도 Role 4가 구현할 helium-side/wall/secondary-side와
coupling을 주입받아 한 대표 모델로 조립한다. 하위 물리식은 아직 없다.
단순 장치는 필요 없는 하위 계층을 만들지
않는다. ThermalSurfacePort는 등록된 스칼라 온도 또는 배열의 한 셀을
공개한다. SystemAssembly는 ThermalLink를 매 잔차 평가에서 한 번 계산해
같은 trial의 열률 배열을 양쪽 Component에 전달한다.

`interfaces/`는 상태·별칭·포트·provider 타입과 읽기 전용 context를,
`system/`은 연결, 등록 완료 후 포트 검증, 잔차 조립과 결합 Jacobian
진단을 맡는다. 동일 유체·열률 포트의 중복 연결을 거부하고, 해석 전
결합 Jacobian의 국소 정칙성을 검사할 수 있다. `solvers/`의 DASSL은 index-1 DAE 일치 초기화와
가변 스텝 BDF 1~5차를 사용한다. Python에 Modelica의 자동 flattening
또는 index reduction은 없다. 모델 개발자가 미지수·잔차 ledger와
독립성을 직접 검증해야 한다. 확인된 Jacobian 비영 패턴을 주면
색칠 유한차분과 희소 선형 풀이를 사용하며, 기본 경로는 밀집 행렬이다.

`tests/fixtures/`의 숫자와 provider는 합성 시험 전용이다. 실제 HTGR
설정에서 자동으로 선택되지 않는다. 실제 형상·운전점·물성·시나리오가
없어 현재 `main.py --demo`만 계산 가능하다. 모델 구현 게이트는
[통합 계획](integration_plan.md), 실제 파일별 상태는
[구현 현황](implementation-status.md)을 따른다.
