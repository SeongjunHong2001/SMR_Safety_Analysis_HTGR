# Role 3 노심 고체·출력·열전달 인수인계

브랜치는 `role3/heat-transfer-base`이고 `role2/helium-properties`(`95ddb10`) 위에 있다.
시험에서 Role 2 헬륨 물성과 Role 1 입력을 쓰므로 Role 1 → Role 2 → Role 3 순서로 병합한다.

## 진행 상태

| 단계 | 파일 | 상태 |
| --- | --- | --- |
| 공통 열전달 커널 | `physics/heat_transfer/{conduction,convection,radiation}.py` | 구현 (이 문서) |
| 헬륨 대류 상관식 | `physics/heat_transfer/correlations.py` | 구현: `mode_nu = "mceligot"` (이 문서) |
| 노심 고체 | `models/core/solid.py` | 미구현 (Protocol만) |
| 노심 출력 | `models/core/power.py`, `physics/neutronics.py` | 미구현 |
| 고체–유체 결합 | `models/core/coupling.py` | 미구현 |

## 공통 열전달 커널

`physics.heat_transfer`의 고정 서명(`base.py`)을 그대로 쓴다. 열률은 a → b 방향이 양수 [W], 온도는 절대 K, 입력은 NumPy broadcasting.

| 함수 | 식 |
| --- | --- |
| `conduction(T_a, T_b, k, A, L)` | k A (T_a − T_b) / L (원통 등은 호출자가 등가 A/L을 넘김) |
| `convection(T_a, T_b, h_conv, A)` | h_conv A (T_a − T_b) |
| `radiation(T_a, T_b, epsilon, A, sigma)` | ε σ A (T_a⁴ − T_b⁴), ε는 두 면의 유효 방사율 |
| `heat_transfer_coefficient(Nu, k, D_h)` | Nu k / D_h |

- **보조:** `radiation.STEFAN_BOLTZMANN` = 5.670374419×10⁻⁸ (CODATA 2018)
- **보조:** `radiation.effective_emissivity_gray_enclosure(eps_a, eps_b, A_a/A_b)` = 1 / (1/ε_a + (A_a/A_b)(1/ε_b − 1))
  (동심 원통·구, 평행평판). RPV ε 0.85, RCCS ε 0.80, r 3.408 / 4.633 m이면 0.735로 변수표 `epsilon_rccs_effective` 0.74와 맞는다.
- **입력 검사:**
  - 면적·길이·D_h·σ가 음수 또는 비유한이면 ValueError (설정 오류). 면적 0은 허용.
  - 계수(k, h, Nu, ε)가 유한한 음수면 ValueError. ε는 1을 넘으면 ValueError.
  - 계수가 NaN이면 그대로 NaN을 전달한다 (물성 provider가 비물리 상태에 NaN을 주는 정책과 같음).
  - 온도가 0 이하이거나 비유한이면 해당 원소만 NaN. Newton 선탐색이 거부한다.

## 헬륨 대류 상관식 (`McEligotConvection`)

팀 보고서 [2] 「헬륨 유동 및 열전달」 파트 C §3의 결정 트리를 구현했다. 매개변수 이름은 변수표와 같다.

```python
from physics.heat_transfer.correlations import make_convection_correlation
corr = make_convection_correlation("mceligot", include_nc=True)      # 다른 mode_nu는 ValueError
h = corr.coefficient(state, geometry)                                  # W/(m²·K)
diag = corr.evaluate(state, geometry)   # Re, Re_safe, Pr, w_tr, Nu_lam, Nu_turb, Nu_fc, Gr, Ra, Nu_nc, Ri, Nu, h_conv
```

| 입력 | 키 | 단위 |
| --- | --- | --- |
| state | `m_dot` (부호 허용), `T_bulk`, `T_wall`, `mu`, `k`, `cp`, `rho` (include_nc일 때) | kg/s, K, K, Pa·s, W/(m·K), J/(kg·K), kg/m³ |
| geometry | `D_h`, `A_flow` (m_dot이 지나는 면적), `shape` = `circular`(기본) \| `parallel_plates` | m, m² |

- **난류:** Nu_turb = 0.021 Re_safe^0.8 Pr^0.4 r_eff^−0.5 (R1 Table AIV.13). 냉각 구간(T_w < T_b)은 r_eff가 1로 매끄럽게 수렴한다 (`eps_T_ratio_clamp` 0.01).
- **층류 바닥:** 원관 4.364, 평행평판 8.235.
- **천이:** w_tr = ½[1 + tanh((Re_safe − 3150)/350)], Re_safe = √(Re² + 100²).
  보고서는 tanh 안에 Re_b를 쓰지만, 역류 대칭과 매끄러움을 위해 Re_safe를 쓴다. 그래서 Re = 3150에서 w_tr이 0.502다 (보고서 0.500).
- **자연대류:** 수직 유로 Churchill–Chu, 대표길이 D_h, β = 1/T_b, |ΔT|_sm = √(ΔT² + 0.1²). 결합은 Nu = (Nu_fc⁶ + Nu_nc⁶)^(1/6).
- **비보수 주의:** 관내 자연순환 유량은 1D 운동량식이 이미 계산하므로, P-LOFC 안전해석은 `include_nc=False`(층류 바닥만) 민감도를 함께 낸다.
- **적용 표시:** `evaluate(...)['mceligot_nominal']`는 McEligot 원 적용 조건(|Re| ≥ 10⁴, T_w/T_b ≥ 1)일 때 True다. 그 밖은 천이·층류·자연대류 결합 결과다.
- **미구현:** Taylor·Gnielinski·Dittus–Boelter 민감도 모드, 고온 벽면 층류화(부력 열화)는 2차 과제다 (보고서 [2] C).
  마찰계수(Churchill)는 Role 2 `physics/pressure_loss.py` 범위라 여기 두지 않았다.

## 공유 파일 변경 (Role 5 확인 필요)

1. **`physics/heat_transfer/__init__.py` 버그 수정:** 공개 함수 이름(`conduction`, `convection`, `radiation`)이 같은 이름의
   하위 모듈과 겹친다. 하위 모듈이 처음 import되면 패키지 속성이 함수에서 모듈로 바뀌어, 두 번째 호출부터
   `'module' object is not callable` 오류가 났다 (커널 구현 후 처음 드러남). 하위 모듈을 먼저 불러 둔 뒤 base 함수를
   바인딩하도록 두 줄을 추가했다. 공개 이름과 서명은 그대로다.
2. **`tests/contracts/test_contract_v3.py`:** "미구현이면 NotImplementedError" 검사에서 열전달 4개 함수를 빼고,
   유한한 양의 값을 두 번 연속 돌려주는지 검사로 바꿨다. 다른 Component의 미구현 검사는 그대로다.

## 변수표 patch

`docs/registry_patches/role3/variable_patch.csv` 14행 (update 13, add 1).
- **update:** `h_conv_core_wall`, `Nu_fluid`, `Re_fluid`, `Pr_he`, `Ri_fluid`(코드 참조), `Re_min_safe`, `Re_tr_mid`, `Re_tr_width`,
  `dT_eps_nc`, `eps_T_ratio_clamp`, `m_smooth_max`, `include_nc`(값 동일, 구현 기본값), `sigma_SB`(5.6704e-8 → CODATA 정확값)
- **add:** `Nu_lam_floor` (4.364 / 8.235)
- `mode_nu`는 role1 patch가 이미 `mceligot`로 연결하므로 중복 행을 만들지 않았다.
`equation_patch.csv`는 만들지 않았다. 커널·상관식은 대수 함수다.

## 시험

```bash
python -m pytest -q tests/physics      # 33 passed (물성 16 + 열전달 17)
python -m pytest -q                    # 136 passed
python main.py --demo                  # 최종 상태 [316.19074506, 306.64335141, 47.73696821]
```

`tests/physics/test_heat_transfer.py`가 확인하는 것:
- **커널:** 식, a → b 부호, 반복 호출 회귀, σ와 변수표 일치, broadcasting, 잘못된 기하·계수 거부, NaN 전달, 유효 방사율 0.735
- **보고서 [2] C 검산값:**
  - Re 5×10⁴, T_w/T_b 1.15·1.30에서 Nu 95.2·89.6
  - 천이 가중치 0.0077 / 0.50 / 0.9923
  - P-LOFC 6.39 MPa, 900 K, ΔT 20·150 K에서 Ra 3.6×10³·2.7×10⁴, Nu_nc ≈ 6.7
- **R1 연료 냉각채널 정격 근처:** Re 3×10⁴–10⁵, Pr 0.65–0.68, h 1500–3000 W/(m²·K)
- **연속성·특수 조건:** 영유량·역류 대칭, 층류 바닥과 형상, 천이·냉각 경계 연속성, 입력 오류

실제 HTGR 열전달 정확도 시험은 아니다.
