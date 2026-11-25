# Role 5 공통 골격 인수인계

계약 버전은 [v3](../interface_spec.md)이다. 기준 `main`은 `d035aaed`,
계약 구현 commit은 `f214c8b`, `main` 병합 commit은 `52016a8`이다.
[PR #5](https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR/pull/5)에 변경 범위와 시험 결과가 있다. 후속 작업은 최신 `main`에서 시작한다.

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
- [ThermalLink](../../models/thermal_link.py)는 다른 온도 노드에
  동일 열경계의 양쪽 부호를 전달한다. SG·RPV·RCCS 물리는 Role 4다.
- `physics.heat_transfer` 패키지는 `conduction`, `convection`,
  `radiation`, `heat_transfer_coefficient`의 고정 서명을 제공한다.
  실제 식과 범위 검증은 Role 3/4가 구현한다.

## 시험 사용법

```bash
python -m pytest -q
python -m pytest -q tests/contracts
python main.py --demo
```

`tests/fixtures/synthetic_models.py`는 두 온도+대수 열률의 flat과
composite DAE를 같은 변수 순서로 만든다. 같은 residual/Jacobian/
적분 결과를 검사하고, derived 열률 모델과 온도 궤적을 비교한다.
`tests/fixtures/closed_loop.py`는 합성 정상/과도 압력 ledger에서
rank 5/4/5와 중복 질량식의 비정방 오류를 검증한다. 테스트의 수치와
provider는 실제 모델에 주입하면 안 된다.

## 제한과 다음 연결

실제 형상, 헬륨 EOS, 상관식 적용범위, SG/압축기/배관/플레넘/
RCCS 모델, 사고 이력, 실험 검증 데이터가 미정이다. 기본 YAML의
`null`을 사용 가능한 값으로 간주하지 않는다. 생산용 Core는 필수
하위 모델이 없으면 경로를 명시해 실패하고, 열전달 API도 실제 식이
없으면 `NotImplementedError`를 낸다.

Role 2/5는 `p_sys`, 국부압 기준점, 초기 재고량 식을 실제 셀
질량식과 함께 자유도 ledger로 결정해야 한다. Role 3/4는 각 열경계의
면적·격자 매핑 및 반대 부호 전달, ID·scale·shape를 제출한다.
Role 1의 실제 데이터가 들어오면 G1 단독 시험 후 G2 전역 조립,
일치 초기화, scaled Jacobian, 실제 정상/LOFC 검증을 진행한다.
