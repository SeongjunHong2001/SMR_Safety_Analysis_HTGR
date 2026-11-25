# Role 5 공통 골격 인수인계

계약 버전은 [v3](../interface_spec.md)이다. `d035aaed`는 최초 작성 기준이며,
계약 구현 commit은 `f214c8b`, `main` 병합 commit은 `52016a8`이다.
[PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)에 변경 범위와 시험 결과가 있다. 후속 작업은 최신 `main`에서 시작한다.
STEP 01B는 `ThermalSurfacePort`와 Assembly의 trial별 열링크 전달을 보완했다.
STEP 01C는 HeatExchanger 상위 DI 조립, `HeliumInventory` 합산 계약,
372행 CSV/Excel의 전체 대응 필드 동기화, 합성 저마하수 폐루프의
정상/과도 방정식 수·rank·DASSL 시험을 추가했다.

## 공개 API와 소유권

- `StateRegistry.register/alias/resolve/freeze`, `Component.name/key/ports/
  residual_ids/residual_scales`, `SimulationContext`, `HeliumProperties`가
  [interfaces](../../interfaces/)의 공통 API다.
- [SystemAssembly](../../system/assembly.py)는 대표 Component와
  연결식을 조립해 잔차·ID·scale을 같은 순서로 반환한다.
  `diagnose_jacobian`은 상태/잔차 척도를 적용한 결합 행렬을 보고한다.
- [StreamJunction](../../system/connections.py)은 압력 N−1, 유량 1,
  stream N의 2N 잔차를 가진다. 매질·스칼라 키·별칭 중복을 검증한다.
- `from models.core import Core` 및 기존 대표 장치 import가 유지된다.
  Core는 fluid/solid/power/coupling을 명시적으로 주입받고 Core만
  Assembly에 등록한다. 실제 submodel 공식은 Roles 2/3 소유다.
- `from models.heat_exchanger import HeatExchanger`는 네 구현을 생성자에
  명시적으로 주입하는 상위 Composite다. 하위 등록·잔차·ID·scale·포트·출력
  조립은 완료됐고 SG 물리 모델과 coupling은 Role 4 소유다.
- `interfaces.contracts.HeliumInventory`와
  `system.inventory.total_helium_inventory`는 Component별 체적·재고량을
  현재 trial에서 합산한다. Role 2와 SG 헬륨 측의 Role 4가 자기 체적에 구현한다.
- [ThermalLink](../../models/thermal_link.py)는 다른 온도 노드에
  동일 열경계의 양쪽 부호를 전달한다. `SystemAssembly(..., thermal_links=...)`가
  평가마다 한 번 실행하고 `{link.name}.to_a`/`.to_b`의 W 배열을
  `SimulationContext.heat_rates`에 넣는다. [ThermalSurfacePort](../../interfaces/ports.py)는
  상태 온도 또는 배열 셀을 `T,index`로 참조한다. SG·RPV·RCCS 물리는 Role 4다.
- `physics.heat_transfer` 패키지는 `conduction`, `convection`,
  `radiation`, `heat_transfer_coefficient`의 고정 서명을 제공한다.
  실제 식과 범위 검증은 Role 3/4가 구현한다.

## 시험 사용법

```bash
python -m pytest -q
python -m pytest -q tests/contracts
python scripts/compare_variable_registry.py
python main.py --demo
```

`tests/fixtures/synthetic_models.py`는 두 온도+대수 열률의 flat과
composite DAE를 같은 변수 순서로 만든다. 같은 residual/Jacobian/
적분 결과를 검사하고, derived 열률 모델과 온도 궤적을 비교한다.
`tests/fixtures/closed_loop.py`는 합성 정상/과도 압력 ledger에서
rank 5/4/5와 중복 질량식의 비정방 오류를 검증한다. 테스트의 수치와
provider는 실제 모델에 주입하면 안 된다.
`tests/fixtures/synthetic_heat_exchanger.py`는 SG 상위의 flat/Composite
잔차·Jacobian·DASSL·에너지 시험에만, `low_mach_loop.py`는 7미지수
저마하수 폐루프의 정상/과도 rank·질량·에너지 시험에만 사용한다.

## 제한과 다음 연결

실제 형상, 헬륨 EOS, 상관식 적용범위, SG/압축기/배관/플레넘/
RCCS 모델, 사고 이력, 실험 검증 데이터가 미정이다. 기본 YAML의
`null`을 사용 가능한 값으로 간주하지 않는다. 생산용 Core는 필수
하위 모델이 없으면 경로를 명시해 실패하고, 열전달 API도 실제 식이
없으면 `NotImplementedError`를 낸다.

`p_sys`의 지정 lower-plenum 단일 Differential 소유, 국부 기준식
`p_fluid_cell[lower]-p_sys=0`, `rho(p_sys,h_he)`, 셀 질량식의
`rho_p*p_sys_dot+rho_h*h_dot`, 정상 종속 질량식 한 개의 재고량 또는
기준압 제약 교체는 [인터페이스 계약](../interface_spec.md)에 확정했다.
Role 2/5는 이를 실제 셀·경계에 적용하고 독립 rank와 일치 초기조건을
재검증해야 한다. Role 3/4는 각 열경계의
면적·격자 매핑 및 반대 부호 전달, ID·scale·shape를 제출한다.
Role 1의 실제 데이터가 들어오면 G1 단독 시험 후 G2 전역 조립,
일치 초기화, scaled Jacobian, 실제 정상/LOFC 검증을 진행한다.
