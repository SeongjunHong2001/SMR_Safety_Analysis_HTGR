# Role 3 노심 고체·출력·열전달 인수인계

1단계(열전달 커널·상관식)는 `role3/heat-transfer-base`(`ee4b007`, Role 2 물성 위),
2단계(노심 고체·출력·결합)는 `role3/core-solid-power-coupling`(Role 2 유동·장치 `d7130ec` 위)이다.
병합 순서는 Role 1 → Role 2 물성 → Role 3 열전달 → Role 2 유동·장치 → Role 3 노심 고체다.

## 진행 상태

| 단계 | 파일 | 상태 |
| --- | --- | --- |
| 공통 열전달 커널 | `physics/heat_transfer/{conduction,convection,radiation}.py` | 구현 (이 문서) |
| 헬륨 대류 상관식 | `physics/heat_transfer/correlations.py` | 구현: `mode_nu = "mceligot"` (이 문서) |
| 고체 물성 | `physics/heat_transfer/solid_properties.py` | 구현: R1 Ex. 2b 고정 열전도도 (아래) |
| 노심 고체 | `models/core/solid.py` | 구현: `UnitCellCoreSolid` (아래) |
| 노심 출력 | `models/core/power.py` | 구현: `SpecifiedCorePower` (지정 총출력). 점동특성·붕괴열(`physics/neutronics.py`)은 미구현 |
| 고체–유체 결합 | `models/core/coupling.py` | 구현: `UnitCellCoreCoupling` (아래) |

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

## 노심 고체·출력·결합 (2단계)

### 고체 물성 `R1SolidProperties` (`physics/heat_transfer/solid_properties.py`)

- **열전도도** (`mode="ex2b_fixed"`, R1 Table AIV.1): 컴팩트 20, 균질 연료블록 37, H-451 66, 2020 35, Alloy 800H 17.8 W/(m·K).
  온도·조사량 의존 식(Ex. 2c)은 참가자 한정 조사량 맵이 없어 미구현이다.
- **비열:** 흑연 cp = 4184(c0 + c1T + c2/T + …) (1000 K에서 1759 J/(kg·K), 보고서 [3]과 같음), Alloy 800H 515 J/(kg·K).
- **밀도:** H-451 1850, 2020 1780, Alloy 800H 7800 kg/m³.
  [가정] 컴팩트·균질 연료블록의 밀도·비열은 R1 파일에 없어 H-451 값을 쓴다. 과도 열용량에만 영향이 있다.

### 노심 고체 `UnitCellCoreSolid` (`models/core/solid.py`)

```python
solid = UnitCellCoreSolid(geometry, R1SolidProperties(), helium=he, initial=dict(T=profile),
                          channels=("average", "hot"), rings=True, outer_link=None)
```
- **가열 기둥:** average 65기둥(표준 53 + RSC 12), hot 1기둥(표준). 축 17셀 중 연료 셀(10개)은 컴팩트·블록 흑연 2노드, 나머지는 블록 1노드다.
  셀 재료: CRE는 Alloy 800H, 상·하부 반사체·전이는 H-451, 유량분배·포스트는 2020 [가정, HS-08].
- **반경 링** (HS-03~HS-06): 내측 반사체 r 0–0.824, 외측 교체 반사체 1.743–2.413, 영구 반사체 2.413–2.970, 배럴 2.973–3.049 m.
  평균 기둥 블록 노드가 연료 영역(0.824–1.743 m) 역할로 내·외측 반사체와 원통 껍질 전도를 한다. 영구 반사체–배럴은 3.5 mm 헬륨 간극(전도 + 복사)이다.
  고온 기둥은 반경 연결이 없다 (HS-02 outer: null). 링·블록 모두 축방향으로 전도한다.
- **단위격자 저항** (냉각공 1개당 단위 길이; 표준 블록 연료공 210 : 냉각공 108):
  - 컴팩트 R_c = (N_c/N_f)/(8πk_c)
  - 간극 R_gap = (N_c/N_f)/(2π r_compact h_gap), h_gap = k_He/δ + 4σε_eff T̄³ (δ = 0.1 mm)
  - 흑연 web R_web = 1/(6 k_H451 S_pair), S_pair = 2π/acosh((p² − r_f² − r_c²)/(2 r_f r_c)) (두 평행 원통 형상계수, 이웃 연료공 6개 [가정])
  - 블록 노드는 web 가운데에 둔다: G_fb = N_c Δz/(R_c + R_gap + R_web/2), 벽 쪽 G_bw = N_c Δz/(R_web/2).
- **상태:** `{name}.{유로}.T_fuel` (10,), `{name}.{유로}.T_block` (17,), `{name}.{링}.T` (17,). 모두 differential K이고, 고체 격자 순서는 `grid_order`·`slices`다.
- **에너지식:** ρc_p(T)V dT/dt − heat_rate − 내부 전도 − 배럴 외부 링크 = 0. heat_rate는 Core가 주는 power − Q_from_solid다.
- **외면:** `ports()` = `outer_surface[i]` (배럴 온도 ThermalSurfacePort 17개).
  `outer_link="이름"`이면 `context.heat_rates["이름.to_a"]` (17,) W를 배럴에 더한다 (Role 4 RPV·riser 링크용).
- **출력:** `{name}.{유로}.T_compact_center` (컴팩트 중심 = T_fuel + q'/(8πk_c)), `{name}.T_fuel_max`, 링 온도, 전체 열용량.
  T_fuel_max는 컴팩트 중심 기준이고 TRISO 커널 온도가 아니다.

- **S0 core_benchmark 전체 조립에서 (Role 4 단계):** 배럴 외면은 `outer_link="barrel_riser"`로 riser 헬륨과 연결된다.
  배럴 평균 569 K (R3 595 K), 연료 최고 1132 K. `R1SolidProperties`에 RPV 재료 `pressure_vessel`과 `fixed_conductivity('air')`를 더했다.

### 노심 출력 `SpecifiedCorePower` (`models/core/power.py`)

- **배분:** P_i(t) = P_total(t) × channel_fraction × axial_fraction을 컴팩트 노드에만 넣는다 (operation.yaml `reactor.power_distribution`, 균일 분포 [가정]).
- **입력 확인:** 분율 합이 1이 아니거나, 연료 없는 유로·셀에 분율이 있으면 오류다.
- **지원 모드:** `total_power`는 숫자 또는 t → W 함수다. `mode_core_power="fission_plus_decay"`는 붕괴열(D7)·점동특성 계수가 없어 `ModelInputMissing`을 낸다.

### 고체–유체 결합 `UnitCellCoreCoupling` (`models/core/coupling.py`)

- **가열 기둥:** Q = (T_block − T_he) G_bw G_h/(G_bw + G_h), G_h = h A_w, A_w = 4A/D_h × Δz.
  h는 `McEligotConvection`이고 벽온도를 두 번 평가한다 (T_w = T_b로 h₀ → 직렬 분배 T_w → h). 캐시는 쓰지 않는다.
- **우회:** Type I 영역 분율로 우회 젖은 면적을 나눈다 [가정]: 노심내 → 평균 연료 블록, 내측(반사체·CR) → 내측 반사체,
  외측(CR·1·2열) → 외측 반사체, 영구 → 영구 반사체. Q = h A_w,share (T_solid − T_bypass).
- **반환:** 반환값 (Q_to_fluid, Q_from_solid)은 합이 같다. Core가 이를 검사하고 유체·고체에 한 번씩 넘긴다.

### Core 조립 정격 결과 (`tests/components/test_core_role3.py`)

조건: 유체 3유로 + 고체 + 출력 + 결합, 시험 경계, 157.1 kg/s, 259 °C, 350 MW 균일, 우회 D_h 4 mm [시험 가정]. 미지수 288개, 정상해 약 4 s.
- **헬륨 열수지:** 외면 단열이라 350 MW가 모두 헬륨으로 가고, 혼합 출구는 961.0 K다.
- **우회 흡수 약 38 MW:** 이때 평균 유로 출구가 960.8 K다. 보고서 [2] C §5.1의 "우회 39.1 MW 흡수 시 유효 채널 출구 687 °C"와 같은 상황이다.
- **연료 최고온도:** 컴팩트 중심 약 1129 K (856 °C). 균일 출력이라 R3(Ex. 2, 실제 분포) 1261 K보다 낮다.
- **고온 유로 출구:** 약 1018 K. 반경 연결이 없어 열이 빠지지 않는다.
- **과도:** 고온 기둥 단독, 출력 20 % 증가 5 s, 희소 Jacobian. 매 시각 Σ C dT/dt = P − Q_to_fluid (상대 1e-6), 연료 온도 상승.

### 남은 일 (2단계)

1. **출력분포:** 실제 R1 분포(참가자 한정 파일) 또는 팀이 정할 분포. 현재 균일 [가정].
2. **배럴 외면 연결 (Role 4·5):** riser(140) 대류·RPV 링크는 `outer_link`로 받는다. core_benchmark의 RCCS 30 °C 경계까지 가려면 Role 4 Vessel·RCCS-0이 필요하다.
3. **점동특성·붕괴열:** `physics/neutronics.py`, `StatefulCorePower`. 사고 시나리오와 함께 진행한다 (D7·D8, 반응도 계수 미확보).
4. **검증 근거:** web 형상계수(n_pair 6)와 노드 위치(web 가운데)는 [가정]이다. R1/R3 단위격자 결과로 확인해야 한다.
5. **Ex. 2c 물성:** 온도·조사량 의존 열전도도는 조사량 맵이 필요하다.

## 공유 파일 변경 (Role 5 확인 필요)

1. **`physics/heat_transfer/__init__.py` 버그 수정:** 공개 함수 이름(`conduction`, `convection`, `radiation`)이 같은 이름의
   하위 모듈과 겹친다. 하위 모듈이 처음 import되면 패키지 속성이 함수에서 모듈로 바뀌어, 두 번째 호출부터
   `'module' object is not callable` 오류가 났다 (커널 구현 후 처음 드러남). 하위 모듈을 먼저 불러 둔 뒤 base 함수를
   바인딩하도록 두 줄을 추가했다. 공개 이름과 서명은 그대로다.
2. **`tests/contracts/test_contract_v3.py`:** "미구현이면 NotImplementedError" 검사에서 열전달 4개 함수를 빼고,
   유한한 양의 값을 두 번 연속 돌려주는지 검사로 바꿨다. 다른 Component의 미구현 검사는 그대로다.

## 변수표 patch

`docs/registry_patches/role3/variable_patch.csv` 29행 (update 26, add 3).
- **update:** `h_conv_core_wall`, `Nu_fluid`, `Re_fluid`, `Pr_he`, `Ri_fluid`(코드 참조), `Re_min_safe`, `Re_tr_mid`, `Re_tr_width`,
  `dT_eps_nc`, `eps_T_ratio_clamp`, `m_smooth_max`, `include_nc`(값 동일, 구현 기본값), `sigma_SB`(5.6704e-8 → CODATA 정확값)
- **add:** `Nu_lam_floor` (4.364 / 8.235)
- **2단계 update 13행:**
  - 고체 상태·출력: `T_core_solid`, `T_core_fuel_max`(컴팩트 중심 기준), `P_core_solid`, `V_core_solid`
  - 간극: `epsilon_core_gap` 0.85, `Q_dot_core_gap`
  - 벽·대류: `T_core_wall`(Derived로 결정), `A_core_wall`, `Q_dot_core_wall`
  - 전도: `Q_dot_solid_cond`, `R_solid_cond`
  - 우회·용기: `Q_dot_core_bypass_wall`(단열 대신 고체 열결합으로 결정), `Q_dot_core_solid_to_vessel`
- **2단계 add 2행:** `n_core_web_pairs` 6 [가정], `map_core_bypass_heat` [가정]
- **중복 피함:** `k_core_solid`, `c_p_core_solid`, `rho_core_solid`, `f_core_power`, `delta_core_gap`, `mode_core_power`는 role1 patch가 이미 다룬다.
- `mode_nu`는 role1 patch가 이미 `mceligot`로 연결하므로 중복 행을 만들지 않았다.
`equation_patch.csv`는 만들지 않았다. 커널·상관식은 대수 함수다.

## 시험

```bash
python -m pytest -q tests/physics      # 58 passed (Role 2·3 물리 커널)
python -m pytest -q tests/components/test_core_role3.py   # 8 passed (약 9 s)
python -m pytest -q                    # 186 passed (약 52 s)
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
