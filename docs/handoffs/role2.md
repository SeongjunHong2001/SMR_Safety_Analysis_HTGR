# Role 2 1차 헬륨 물성·열수력 인수인계

브랜치는 `role2/helium-properties`이고 `role1/mhtgr350-inputs`(`abd6d48`, `main` `8b98291` 기반) 위에 있다.
Role 1 입력(`data/properties/r1_benchmark_properties.yaml`, `helium.csv`, `config/operation.yaml`)을 시험에 쓰므로
Role 1 PR이 먼저 병합돼야 한다.

## 진행 상태

| 단계 | 파일 | 상태 |
| --- | --- | --- |
| 1. 헬륨 물성 provider | `physics/properties.py` | 구현 (이 문서) |
| 2. 1D FVM 보존 커널 | `physics/flow.py` | 미구현 |
| 3. 압력손실 | `physics/pressure_loss.py` | 미구현 |
| 4. 노심 유체 | `models/core/fluid.py` | 미구현 (Protocol만) |
| 5. 플레넘·배관·순환기 | `models/{plenum,piping,compressor}/` | 미구현 |

## 헬륨 물성 provider

`physics.properties.R1HeliumProperties`는 `interfaces.contracts.HeliumProperties`의 9개 메서드를 구현한다.
계수는 R1 Table AIV.13을 옮긴 `data/properties/r1_benchmark_properties.yaml`의 `helium` 절에서 읽는다.
코드 안에 물성 숫자를 따로 두지 않는다.

```python
from physics.properties import R1HeliumProperties, make_helium_properties

he = make_helium_properties("ideal_gas")          # mode_he_eos. 다른 값은 ValueError
he = R1HeliumProperties.from_yaml(T_ref=0.0)      # 같은 결과. 파일 경로를 줄 수도 있다
rho = he.density(p_sys, h_he)                     # SI, 스칼라·배열 broadcasting
```

| 메서드 | 식 | 단위 |
| --- | --- | --- |
| `temperature(p, h)` | T = T_ref + h / c_p | K |
| `enthalpy(p, T)` | h = c_p (T − T_ref) | J/kg |
| `density(p, h)` | ρ = p / (R T) | kg/m³ |
| `viscosity(p, T)` | μ = 3.674×10⁻⁷ T^0.7 | Pa·s |
| `conductivity(p, T)` | k = 2.682×10⁻³ (1 + 1.123×10⁻⁸ p) T^(0.71 (1 − 2×10⁻⁹ p)) | W/(m·K) |
| `cp(p, T)`, `cv(p, T)` | c_p = 5195, c_v = c_p − R = 3117.54 | J/(kg·K) |
| `drho_dp(p, h)` | ∂ρ/∂p│h = 1 / (R T) | kg/(m³·Pa) |
| `drho_dh(p, h)` | ∂ρ/∂h│p = −p / (R T² c_p) | kg²/(m³·J) |

- **기체상수:** R = 2077.46 J/(kg·K) (= 8.314 / 0.004002). 팀 결정(2026-11-26)대로 R1 세트를 전 계통에 쓴다.
- **엔탈피 기준:** `T_ref` 기본값은 변수표 `T_he_enthalpy_ref` 제안값 0 K다. 기준을 바꾸면 h 값만 평행이동하고 같은 (p, T)의 밀도·도함수는 같다.
- **출력 형태:** 입력의 broadcast shape이다. 스칼라 입력은 NumPy 스칼라로 돌려준다.
- **비물리 입력:** p ≤ 0, T ≤ 0(h로 계산한 T 포함), 비유한 값은 그 원소만 NaN이다. 예외나 경고를 내지 않는다.
  Newton 선탐색(`solvers/newton.py`)이 비유한 잔차를 거부하므로 시도 상태가 잠시 벗어나도 적분이 멈추지 않는다.
- **적용 범위 [가정]:** T 300–1900 K (보고서[1] §6.1), p 0.05–8 MPa (시나리오 범위 0.1 MPa 감압 ~ 7.18 MPa 안전밸브 설정압 포함).
  범위 밖 양수 입력도 계산은 한다. 확인은 `in_valid_range(p, T)` (bool 배열), `check_valid_range(p, T, where=...)` (ValueError)로 한다.
  초기화와 결과 진단에서 호출하고 잔차 계산 안에서는 호출하지 않는다.
- **이상기체 가정:** R1 벤치마크 정의를 따른다. 실기체 보정(KTA 등)은 `mode_he_eos`에 새 모드로 추가할 수 있게 팩토리를 두었다.

## 다른 담당이 쓰는 법

- **Role 2 유체 모델(노심 채널, 플레넘, 배관, SG He 측):** EOS는 전역 `p_sys`와 셀 `h_he`로 `density`, `drho_dp`, `drho_dh`를 부른다.
  질량식의 `V (∂ρ/∂p · ṗ_sys + ∂ρ/∂h · ḣ)` 항에 그대로 쓴다.
- **Role 3 대류 상관식:** 벌크 온도에서 `viscosity`, `conductivity`, `cp`로 Re·Pr을 만든다. R1 Nu 식(McEligot)은 Role 3 범위다.
- **Role 5 builder:** `operation.yaml`/`meta.helium_property_set`과 변수표 `mode_he_eos`로 `make_helium_properties(mode)`를 한 번 만들어 모든 Component에 같은 객체를 주입한다.

## 변수표 patch

`docs/registry_patches/role2/variable_patch.csv` 6행 (모두 update): `T_he`, `T_he_enthalpy_ref`, `rho_he`, `drho_dp_he`, `drho_dh_he`, `k_he`.
값은 바꾸지 않고 코드 참조와 식을 연결한다. `k_he`의 상태 "KTA 채택 제안"은 R1 식 채택으로 바꾼다.
`c_p_he`, `c_v_he`, `R_he`, `mu_he`, `mode_he_eos`는 role1 patch가 이미 R1 값과 물성 파일로 연결하므로 중복 행을 만들지 않았다.
`docs/registry_patches/role2/equation_patch.csv`는 만들지 않았다. 새 방정식이 없다 (물성은 대수 함수).

## 시험

```bash
python -m pytest -q tests/physics     # 16 passed
python -m pytest -q                   # 119 passed (main 65 + Role 1 입력 38 + Role 2 물성 16)
python main.py --demo                 # 최종 상태 [316.19074506, 306.64335141, 47.73696821]
```

`tests/physics/test_helium_properties.py`가 확인하는 것:
Protocol 만족, 파일 계수와 R = R_u/M·c_v = c_p − R 일관성, Role 1 `helium.csv` 51행 일치, R1 원식(bar 단위) 재현,
엔탈피 왕복과 이상기체 관계, 엔탈피 기준 변경 불변성, 밀도 도함수와 중앙차분 일치, broadcasting·스칼라 출력,
비물리 입력의 NaN 처리(경고 없음), 적용 범위 판정, `operation.yaml` 정격 분포 10개 지점의 범위·밀도 일치, 팩토리 모드 선택, 잘못된 계수 거부.
실제 HTGR 열수력 정확도 시험은 아니다.
