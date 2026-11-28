# Role 2 1차 헬륨 물성·열수력 인수인계

브랜치는 두 개다. 1단계 `role2/helium-properties`(`95ddb10`)는 `role1/mhtgr350-inputs`(`abd6d48`) 위에,
2·3단계 `role2/flow-pressure-loss`는 `role3/heat-transfer-base`(`ee4b007`) 위에 있다.
시험에 Role 1 입력과 Role 2 물성을 쓰므로 Role 1 → Role 2 물성 → Role 3 → Role 2 유동 순서로 병합한다.

## 진행 상태

| 단계 | 파일 | 상태 |
| --- | --- | --- |
| 1. 헬륨 물성 provider | `physics/properties.py` | 구현 (이 문서) |
| 2. 1D FVM 보존 커널 | `physics/flow.py` | 구현 (아래) |
| 3. 압력손실 | `physics/pressure_loss.py` | 구현 (아래) |
| 4. 노심 유체 | `models/core/fluid.py` | 구현: `ParallelChannelCoreFluid` (아래) |
| 5. 플레넘·배관·순환기 | `models/{plenum,piping,compressor}/` | 구현 (아래). 밸브(MLSV)는 미구현 |

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

## 1D FVM 보존 커널 (`physics/flow.py`)

팀 보고서 [2] 파트 B §3(엇갈림 격자, C∞ 상류차분, DASSL 잔차)을 Registry를 모르는 순수 함수로 구현했다.

- **격자:** 셀량 (N,), 면 유량 (N+1,). 면 유량은 정격 유동 방향이 양수, 면 0이 장치 입구, 면 N이 출구다.
  포트 유량은 장치로 들어오면 양수이므로 `port_inflows_from_faces(m_faces)` = (+ṁ_0, −ṁ_N).
- **잔차:**

| 함수 | 식 | 단위 |
| --- | --- | --- |
| `mass_residual(V, drho_dp, drho_dh, p_sys_dot, h_dot, m_in, m_out)` | V(ρ_p ṗ_sys + ρ_h ḣ) − (ṁ_in − ṁ_out) | kg/s |
| `energy_residual(rho, V, h, h_dot, p_sys_dot, m_in, h_in*, m_out, h_out*, q_wall)` | ρVḣ − Vṗ_sys + ṁ_out(h*_out − h) − ṁ_in(h*_in − h) − Q | W |
| `momentum_residual(p_left, p_right, dp_acc, dp_fric, dp_form, dp_grav, dp_comp, include_inertia=False, l_over_A, m_dot_dot)` | (p_R − p_L) + ΔP_acc + ΔP_fric + ΔP_form + ΔP_grav − ΔP_comp [+ (l/A) dṁ/dt] | Pa |
| `cell_residuals(V=, rho=, drho_dp=, drho_dh=, h=, h_dot=, p_sys_dot=, q_wall=, m_faces=, h_inlet_boundary=, h_outlet_boundary=, m_eps=)` | 한 유로 N셀의 질량·에너지 잔차와 면 엔탈피를 한 번에 | dict |

- **면 엔탈피:** `smooth_upwind`/`face_enthalpies`는 h* = w h_L + (1 − w) h_R, w = ½[1 + tanh(ṁ/ṁ_ε)]이다.
  면 0의 왼쪽 값은 입구 포트 stream 엔탈피, 면 N의 오른쪽 값은 출구 포트(역류 시 유입) stream 엔탈피로 받는다.
- **스무딩 폭:** ṁ_ε = `pressure_loss.mdot_epsilon(ṁ_nom,path)` = 10⁻⁴ ṁ_nom,path (변수표 `eps_mdot_rel`). 경로마다 정격 유량으로 정한다.
- **관성:** 기본 `delta_inert = false`. 면 유량은 `face_flow_state_kind(False)` = `"algebraic"`로 등록한다. 관성을 켜면 `"differential"`이 되고 `face_inertia_length`로 l/A를 만든다.
- **보존:** `R_energy + h R_mass`가 보존형 에너지식과 항등이다 (시험으로 확인). Q_wall은 유체로 들어오면 양수.

## 압력손실 (`physics/pressure_loss.py`)

팀 보고서 [2] 파트 C §2와 파트 B §3.2를 구현했다. 압력 항은 정격 유동 방향 기준으로 하류 압력을 낮추면 양수다.

| 함수 | 식 |
| --- | --- |
| `friction_dp(m_dot, rho, mu, length, D_h, A, roughness=0, phi_s=1)` | (f Re) μ L ṁ / (2 ρ D_h² A), Churchill (1977). ṁ = 0에서 층류 선형식으로 매끄럽게 |
| `face_friction_dp(m_dot, left, right)` | 인접 셀의 반 길이씩 각 셀 물성·기하로 합산 |
| `churchill_f_re(Re, rel_roughness, phi_s)`, `darcy_friction_factor(...)` | f·Re (Re = 0에서 64/Φs), Darcy f (진단) |
| `form_loss_dp(m_dot, rho, A_ref, K_forward, K_reverse=None, m_eps=)` | K_eff ṁ \|ṁ\|_sm / (2 ρ A²), 정·역류 K를 tanh로 전환 |
| `sudden_expansion_K`, `sudden_contraction_K` | (1 − β)², (1/μ_c − 1)², μ_c = 0.62 + 0.38 β³ |
| `gravity_dp(rho, dz_elev)`, `face_gravity_dp(...)` | ρ g Δz_elev, 반셀 가중 ½ g (ρ_i Δz_i + ρ_{i+1} Δz_{i+1}) |
| `acceleration_dp(m_dot, A, rho_L, rho_R)`, `junction_acceleration_dp(...)` | ṁ²/A²(1/ρ_R − 1/ρ_L), ½(ṁ²/(ρ_R A_R²) − ṁ²/(ρ_L A_L²)) |
| `mdot_epsilon(m_nominal_path, eps_rel=1e-4)` | 경로별 스무딩 폭 |

- **형상계수 Φs:** 원관 1 (`PHI_S_CIRCULAR`), 평행평판 틈새 2/3 (`PHI_S_PARALLEL_PLATES`).
- **마찰 모드:** `mode_friction = "churchill"`만 지원한다. 다른 값은 ValueError.
- **f·Re 형식:** 보고서는 Re ≤ 100에서 선형 점성식으로 분기하라고 제안했다. 이 구현은 f·Re = 8[(8/Φs)¹² + Re¹²(A+B)^−1.5]^(1/12)가
  Re → 0에서 유한하다는 점을 이용해 분기 없이 같은 값을 낸다 (0 나눗셈 없음, Re ≤ 100에서 64/Φs와 상대오차 1e-9 이내).
- **입력 정책:** 길이·면적·D_h·K가 음수 또는 비유한이면 ValueError (설정 오류). ρ·μ가 NaN이면 NaN을 전달하고, 유한한 0 이하면 ValueError.
- **보고서 [2] C 검산 재현:**
  - **f:** 0.0203 / 0.0214 / 0.0223 / 0.0225 / 0.0232 (Re 55,490 / 43,800 / 36,710 / 35,330 / 30,760)
  - **거칠기 민감도:** ε 0.01 mm에서 +11.0 %
  - **대형 냉각채널 (D 15.88 mm, 20.58 g/s, 7.93 m, 40분할):**
    - 마찰 14.04 kPa (보고서 14.18, −1.0 %)
    - 열가속 1.50 kPa
    - 중력 −0.33 kPa
  - **급축소 극한 K:** 0.376

## 노심 유체 (`models/core/fluid.py`)

`ParallelChannelCoreFluid`가 `CoreFluid` 계약과 `HeliumInventory`를 구현한다. 상위 `Core`에 `fluid=`로 주입한다.

```python
from models.core.fluid import ParallelChannelCoreFluid, initial_from_operation
fluid = ParallelChannelCoreFluid(geometry, he, initial=initial_from_operation(operation),
                                 mode_bypass="type_I_fixed_fraction",
                                 channel_overrides={"bypass": {"hydraulic_diameter_m": ..., "shape": "parallel_plates"}})
core = Core(fluid=fluid, solid=..., power=..., coupling=...)
```

- **유로:** `geometry.core.channels`의 average·hot·bypass를 병렬로 두고, 세 유로 모두 `core.elevation_profile_m`(18면, 17셀)을 쓴다.
  셀 면적은 `cell_flow_area_m2`(없으면 `flow_area_m2`), 수력직경은 `cell_hydraulic_diameter_m`(없으면 `hydraulic_diameter_m`),
  형상은 `shape`(circular | parallel_plates), 거칠기는 `roughness_m`(기본 0). 값이 null이면 `config/geometry.yaml: missing required input core.channels.<유로>.<키>` 오류.
- **등록 상태** (`{name}` 기본 `core.fluid`):

| 키 | kind | shape |
| --- | --- | --- |
| `{name}.{유로}.h_he` | differential | (17,) |
| `{name}.{유로}.p_fluid_cell` | algebraic | (17,) |
| `{name}.{유로}.m_dot_fluid_face` | algebraic (`delta_inert=True`면 differential) | (18,) |
| `{name}.helium_inlet.{p, m_dot, h_outflow, h_instream}`, `{name}.helium_outlet.{...}` | algebraic | 스칼라 |
| `{name}.bypass.dp_bypass_cal` (type_I_fixed_fraction만) | algebraic | 스칼라 |

  `p_sys`는 읽기만 한다 (지정 하부 플레넘이 등록). 포트는 `helium_inlet`(상부 플레넘 150 쪽), `helium_outlet`(하부 플레넘 230 쪽)이고
  Core 밖에서는 `fluid.helium_inlet/outlet`이다.
- **잔차** (ID `{name}.{유로}.mass[i]`, `.energy[i]`, `.momentum[j]`, `{name}.helium_*.flow_sum`, `.outflow_enthalpy`, `{name}.bypass.type_I_fraction`):
  - 유로마다 질량 17 + 에너지 17 + 운동량 18
  - 포트마다 유량 합(포트 유량 = 유로 경계면 유량 합)과 유출 엔탈피 혼합(나가는 유로 유량으로 매끄럽게 가중) 2개
  - 나머지 포트 식 2개는 접속부(StreamJunction)가 준다.
- **운동량 경계:** 플레넘–유로 면은 베르누이 동압(`plenum_boundary_acceleration_dp`) + 정·역류 K로, 입구 정압강하 1.5 q, 출구 순변화 0
  (보고서 [2] D §4.3, 기본 `inlet_K=(0.5, 1.0)`, `outlet_K=(1.0, 0.5)`). 내부 면은 열가속 + 단면 변화(`face_acceleration_dp`), 반셀 마찰·중력, 선택적 `face_K`.
- **우회:** `mode_bypass="type_I_fixed_fraction"`(S0 벤치마크)은 분율식 ṁ_bypass,0 − f Σṁ_k,0 = 0과 미지수 `dp_bypass_cal`을 더한다
  (f 기본 = Σ `type_I_fraction` = 0.11). `"calibrated_parallel"`은 `bypass_K_cal`을 필수로 받아 우회 출구 면 K에 더하고 압력차로 유량을 푼다.
- **열률:** `residual(..., heat_rate=)`의 heat_rate는 유체로 들어오는 W 배열 (51,)이고 순서는 `cell_order` (average 17 → hot 17 → bypass 17)다.
- **Role 3 결합용:** `coupling_state(t, y, context)`는 셀 순서대로 `T_bulk, rho, mu, k, cp, m_dot(셀 양 면 평균), D_h, A_flow, dz`
  (각 (51,)), `p_sys`, `shape`(유로 형상 문자열)를 돌려준다. `McEligotConvection.coefficient(state, geometry)`에 바로 넣을 수 있다.
- **출력:**
  - `{name}.{유로}.T_he`, `{name}.{유로}.m_dot_inlet`
  - `{name}.T_outlet_mixed`, `{name}.dp_plenum_to_plenum`, `{name}.bypass_fraction`, `{name}.dp_bypass_cal`
- **정격 시험 결과** (실제 geometry, R1 물성, 157.1 kg/s, 259 °C, 350 MW 균일, 우회 D_h 4 mm [시험 가정]):
  - **출구 온도:** 혼합 961.0 K. 평균 유로 출구 1014.2 K로, 보고서 [2] C §5.1의 "단열 우회 시 유효 채널 1014 K"와 같다.
  - **우회:** 분율 0.110, `dp_bypass_cal` 13.8 kPa
  - **플레넘 간 Δp:** 24.5 kPa (R3 19.95 kPa)

### 남은 입력·후속 (노심 유체)

1. **우회 유로 수력직경 (Role 1):** geometry의 `core.channels.bypass.hydraulic_diameter_m`이 null이다 (후보 4 mm 간극, 102 mm 제어봉 채널).
   정해야 실제 실행이 된다. 시험은 4 mm 평행평판으로 가정했다.
2. **하부 블록 수력직경 (Role 1):** 셀 13–16(하부반사체·전이·유량분배·포스트)도 연료 냉각공 D_h 15.7 mm로 계산돼 마찰이 과대평가된다.
   `cell_hydraulic_diameter_m`(전이 68 mm, 유량분배 67 mm, 포스트 116 mm 원 등)을 넣으면 Δp가 R3 값에 가까워질 것으로 본다.
3. **플레넘 연결:** 상부(150)·하부(230) 플레넘 Component와 StreamJunction이 다음 단계(5)다. 그때 `p_sys` 소유와 국부 기준식이 생긴다.
4. **Jacobian 비용 (Role 5):** 노심 유체만 166 미지수다. 밀집 유한차분 Jacobian으로는 3유로 과도 2 s에 약 50 s가 걸렸다
   (시험은 고온 유로 1개로 줄임). 실제 builder는 `jacobian_sparsity`(희소 경로, #10)를 넘겨야 한다.

## 플레넘·배관·순환기 (5단계)

### 배관 `Piping` (`models/piping/piping.py`)

`Piping(nodes["hot_duct"], he, name="hot_duct", initial=dict(m_dot, T, p), inlet=("plenum", 0.5, 1.0), outlet=("pipe", 0, 0), heat_source=None)`
- **config:** `flow_network.nodes.<이름>`의 `cells`, `flow_area_m2`, `hydraulic_diameter_m`, `length_m`, `z_inlet_m`, `z_outlet_m`
  (선택 `roughness_m`, `shape`). null이면 `ModelInputMissing: config/geometry.yaml: missing required input flow_network.nodes.<이름>.<키>`.
- **상태·잔차:** 노심 유체 한 유로와 같은 1D FVM이다 (질량 N, 에너지 N, 운동량 N+1, 포트 4).
- **경계 면:** `("plenum", K_f, K_r)`는 플레넘 베르누이 ±q + K, `("pipe", K_f, K_r)`는 K만 넣는다.
- **벽 열:** 기본 단열. `heat_source(t, y, context, pipe)`가 (N,) W를 주면 에너지식에 넣는다. 루프 시험의 SG-0 모양 열침원이 이 경로다.
- **구현 정리:** 노심 유체와 배관은 같은 면 압력 항 함수 `pressure_loss.channel_face_terms`를 쓴다 (리팩터링 뒤 노심 유체 시험 결과 동일).

- **열 링크용 추가 (Role 4 단계에서):** `Piping.z_faces`(셀 경계 고도), `Piping.coupling_state()`(노심 유체와 같은 키). riser 벽 열은 `heat_source`가 링크 열률을 읽어 받는다.
- **SG He 측용 추가 (Role 4 SG-0 단계):** config `z_faces_m`(비균일 수직 셀, `cells`·`length_m`·`z_*` 대신), `residual_with_heat(t, y, ydot, context, q)`(상위 조립기가 넘긴 벽 열).

### 플레넘 `Plenum` (`models/plenum/plenum.py`)

`Plenum(nodes["lower_plenum"], he, name="lower_plenum", initial=..., port_labels=("from_core", "to_hot_duct"), owns_p_sys=True, steady_pressure_closure=None)`
- **포트:** 실제 접속 수만큼 `helium_ports[i]`. 고도는 `port_elevations`(기본 `[z_inlet_m, z_outlet_m]`)이고, 라벨로 `plenum.port("to_core")`처럼 찾는다.
- **잔차:**
  - 질량 V(ρ_p ṗ_sys + ρ_h ḣ) − Σṁ_i, 에너지 ρVḣ − Vṗ_sys + Σ s⁺(ṁ_i)(h − h_instream,i) − Q
  - 포트마다 압력 p_port,i = p + ρ g (z_ref − z_i)와 유출 엔탈피 h_outflow,i = h
  - z_ref는 `z_bottom_m`·`z_top_m`의 중간이다.
- **지정 하부 플레넘 (`owns_p_sys=True`):**
  - 전역 `p_sys`를 differential·owner `primary_loop`로 한 번 등록하고, 국부 기준식 `{name}.pressure_reference` (p − p_sys = 0)를 더한다.
  - 정상 초기화용 인스턴스에 `steady_pressure_closure=p0`을 주면 자기 질량식을 p_sys − p0 = 0으로 바꾼다 (D2 기준압 폐합).
  - 과도 인스턴스는 원래 질량식을 쓴다.
- **체적 미확보 노드 (M10):** 110, 300, 330, 350은 `volume_m3`가 null이라 `ModelInputMissing`이 난다 (시험으로 확인).

### 순환기 `Compressor` (`models/compressor/`)

체적 없는 두 포트 장치(`helium_inlet` 흡입, `helium_outlet` 토출). `config`는 operation.yaml `compressor` dict.
- **`mode="prescribed_flow"` (S0 전용, 팀 결정):** `m_dot_set`과 `m_rated`가 필수다. 미지수 `{name}.dp_comp`와 유량식을 더하고, 승압은 결과로 나온다.
- **`mode="performance"`:** `performance`(예: `QuadraticPerformance`)와 `speed_schedule(t)`가 필수다. 유량은 루프 운동량이 정한다.
  `QuadraticPerformance`는 ΔP = (ρ/ρ_R) ΔP_R (c0 α² + c1 α ν − c2 ν|ν|)이다. 변수표 부호 규약(c0 + c1 − c2 = 1)을 따르고 계수 기본값은 없다.
- **`mode="coastdown_inertia"`:** 회전자 관성 `I_circ_rotor`가 미확보(M09)라 `ModelInputMissing`을 낸다.
- **잔차:** 질량 ṁ_in + ṁ_out = 0, 운동량 p_out − p_in − ΔP = 0, 양방향 유출 엔탈피 2개(+ 유량식).
- **헬륨 가열 [가정]:** Q = Q_rated (|ṁ|_sm/ṁ_rated)³, Q_rated는 `rated_heat_to_helium_W` 3.1 MW. 정격 승온은 3.8 K로 PSID와 같다.

### Role 2 장치만으로 닫은 1차 루프 시험 (`tests/components/test_primary_loop_role2.py`)

- **구성:**
  - 실제 노드: 하부 플레넘(p_sys 소유) → hot duct → [SG 대용 배관, SG-0 모양 열침원 물 193 °C·UA 1.7 MW/K]
    → [340 대용 배관] → 순환기(유량 지정) → [110 대용 배관] → cold annulus → RV 입구 플레넘 → riser → 상부 플레넘
    → 노심 유체(평균·고온 유로) → 하부 플레넘
  - 연결은 실제 `StreamJunction` 11개, 고도 체인은 geometry 값으로 닫힌다. [대용] 3개는 Role 4 SG와 미확보 노드(M10·M11) 자리다.
- **정상 결과** (260 미지수, 약 3 s):
  - p_sys = 6.3373 MPa (기준압 폐합), 하부 플레넘 국부압 = p_sys
  - 열수지: 노심 350 + 순환기 3.1 = 열침원 353.1 MW
  - 순환기 승압 41 kPa (형상손실 K 미보정이라 PSID 91 kPa보다 작음)
  - 재고량 약 1,650 kg (대용 체적 포함)
- **과도 결과** (출력 10 % 증가, 2 s, 희소 Jacobian 패턴 사용, 약 25 s):
  - 닫힌 루프 질량이 1e-7 이내로 보존되고, 내부에너지 변화 = ∫(노심 + 순환기 − 열침원) dt (2 % 이내, 사다리꼴 적분)
  - p_sys 상승

### 남은 일 (5단계)

1. **MLSV(노드 350, 역류 차단 체크밸브):** 밸브 모델은 미구현이다. 변수표 `K_valve`, `dp_eps_valve`, `alpha_min_valve`가 있다.
2. **K 보정:** PSID 구간 설계손실(91 kPa)에 맞춰 형상손실 K를 S0 system_loop 정상해에서 역산한다 (`pressure_loss_calibration`, 보고서 식 7).
3. **미확보 체적·면적 (M10·M11):** 110·300·330·340·350 체적과 340 단면적.
4. **순환기 감속:** `I_circ_rotor`, 성능곡선(M09)이 확보되면 `coastdown_inertia`(α 미분 상태)를 구현한다. 사고 시나리오와 함께 진행한다.
5. **희소 Jacobian 패턴:** 루프 시험은 정상해의 유한차분 Jacobian에서 0이 아닌 원소를 패턴으로 썼다. builder(Role 5)는 접속 구조로 패턴을 만드는 편이 안전하다.

## 다른 담당이 쓰는 법

- **Role 2 유체 모델(노심 채널, 플레넘, 배관, SG He 측):** EOS는 전역 `p_sys`와 셀 `h_he`로 `density`, `drho_dp`, `drho_dh`를 부른다.
  질량식의 `V (∂ρ/∂p · ṗ_sys + ∂ρ/∂h · ḣ)` 항에 그대로 쓴다.
- **Role 3 대류 상관식:** 벌크 온도에서 `viscosity`, `conductivity`, `cp`로 Re·Pr을 만든다. R1 Nu 식(McEligot)은 Role 3 범위다.
- **Role 5 builder:** `operation.yaml`/`meta.helium_property_set`과 변수표 `mode_he_eos`로 `make_helium_properties(mode)`를 한 번 만들어 모든 Component에 같은 객체를 주입한다.

## 변수표 patch

`docs/registry_patches/role2/variable_patch.csv` 45행 (update 41, add 4).
- **물성 6행:** `T_he`, `T_he_enthalpy_ref`, `rho_he`, `drho_dp_he`, `drho_dh_he`, `k_he`.
  값은 바꾸지 않고 코드 참조와 식을 연결한다. `k_he`의 상태 "KTA 채택 제안"은 R1 식 채택으로 바꾼다.
- **유동·압력손실 15행:**
  - 상태와 잔차: `m_dot_fluid_face`, `h_fluid_upwind`, `R_fluid_mass`, `R_fluid_energy`, `R_fluid_momentum`
  - 압력손실 항: `dp_fluid_fric`, `dp_fluid_form`, `dp_fluid_grav`, `dp_fluid_acc`
  - 마찰 설정: `f_D_fluid`(상태 "상관식 충돌" → Churchill 채택), `mode_friction`, `Phi_fluid_laminar`(1 / 2/3), `K_form_fluid_rev`
  - 수치 설정: `delta_inert`, `eps_mdot_rel`
- **노심 유체 7행:**
  - update 5행: `h_he`, `p_fluid_cell`, `m_dot_fluid_port`, `m_dot_core_bypass`, `f_core_bypass`
  - add 2행: `dp_bypass_cal`, `K_core_channel_boundary`
- **플레넘·배관·순환기 17행:**
  - 배관: `h_pipe`, `p_pipe`, `m_dot_pipe`, `Q_dot_pipe_wall`
  - 플레넘: `h_plenum`, `p_plenum`, `m_dot_plenum_port`, `mode_plenum`
  - 순환기: `dp_circ`, `mode_circ_curve`, `c0/c1/c2_circ_curve`(예비 계수 유지), `mode_coastdown`
  - 공통: `g` = 9.80665
  - add 2행: `m_dot_circ_set`, `mode_circ_heat`
- **다른 patch와 겹치는 행:** `p_sys`, `V_plenum`, `Q_dot_circ`, `I_circ_rotor`는 role1 patch가 값·코드 참조를 이미 다루므로 중복 행을 만들지 않았다.
  `p_sys` 등록 위치는 `models/plenum/plenum.py` (`owns_p_sys=True`)다.
`c_p_he`, `c_v_he`, `R_he`, `mu_he`, `mode_he_eos`는 role1 patch가 이미 R1 값과 물성 파일로 연결하므로 중복 행을 만들지 않았다.
`docs/registry_patches/role2/equation_patch.csv`는 만들지 않았다. 질량·에너지·운동량 잔차는 방정식표에 이미 있는 저마하수 계약식이고,
이번에는 그 식을 계산하는 함수만 구현했다. 실제 Component(노심 유체·플레넘·배관)를 등록할 때 방정식 ID·잔차 수를 patch로 낸다.

## 시험

```bash
python -m pytest -q tests/physics     # 58 passed (물성 16 + 열전달 17 + 압력손실 14 + 유동 11)
python -m pytest -q tests/components  # 17 passed (노심 유체 9 + 플레넘·배관·순환기·루프 8)
python -m pytest -q                   # 178 passed (약 45 s, 과도 시험 2개가 약 35 s)
python main.py --demo                 # 최종 상태 [316.19074506, 306.64335141, 47.73696821]
```

`tests/physics/test_helium_properties.py`가 확인하는 것:
Protocol 만족, 파일 계수와 R = R_u/M·c_v = c_p − R 일관성, Role 1 `helium.csv` 51행 일치, R1 원식(bar 단위) 재현,
엔탈피 왕복과 이상기체 관계, 엔탈피 기준 변경 불변성, 밀도 도함수와 중앙차분 일치, broadcasting·스칼라 출력,
비물리 입력의 NaN 처리(경고 없음), 적용 범위 판정, `operation.yaml` 정격 분포 10개 지점의 범위·밀도 일치, 팩토리 모드 선택, 잘못된 계수 거부.
실제 HTGR 열수력 정확도 시험은 아니다.

`tests/physics/test_pressure_loss.py` (12개)와 `test_flow.py` (11개)가 확인하는 것:
- **압력손실:**
  - 보고서 [2] C의 f 표·거칠기 민감도·대형 냉각채널 압력강하 분해(마찰·가속·중력) 재현
  - 층류 극한 64/96과 저유량 선형성, ṁ = 0에서 기울기 연속
  - Churchill 천이 구간 연속성: 격자를 반으로 줄이면 기울기 변화도 반으로 줄어듦
  - 반셀 가중, 형상손실의 방향 전환·단조성, 급확대·급축소 K, 중력·가속 부호, 입력 정책
- **유동 (보고서 B §5.3 검증 계획):**
  - ① 17셀 가열 채널 정상해 질량·에너지 잔차 < 10⁻¹⁰
  - ② 밀폐 정적 가열 Vṗ_sys = (γ−1)Q = 0.6664 Q
  - ③ 유동 역전 연속성
  - 보존형 항등, 역류 시 하류 경계 엔탈피, 포트 부호, 운동량 잔차·관성 스위치, shape 검사

`tests/components/test_core_fluid.py` (9개)는 실제 geometry·R1 물성과 시험 전용 경계(포트 조건·p_sys 이완식)·지정 열률로 확인한다:
- **구조:** 미지수 166개 정방, 포트 이름, `cell_order`, 정상·과도 반복행렬 full rank
- **정격 정상해:**
  - 열수지(혼합 출구 = 입구 + Q/(ṁ c_p))
  - 우회 0.11
  - 평균 유로 출구 1014 K, 유로 안 유량 일정
  - 고온·평균 유로 출구 차이 < 30 K, 우회 유로 무가열
  - Δp 10–40 kPa 범위
- **재고량·결합 상태:** HeliumInventory 체적이 geometry 체적과 일치, `coupling_state` shape
- **역류:** 하부 플레넘 900 K가 입구 포트로 나감
- **과도:** 출력 10 % 증가 2 s DASSL, 질량 보존 M(t) − M(0) = ∫ 포트 유입
- **조립과 입력 오류:**
  - Core 상위 조립기 주입 (최소 고체·결합)
  - 입력 오류 메시지: 우회 D_h null, `bypass_K_cal` 누락, `initial` 누락, heat_rate shape

## 면 압력 항 공개 (2026-11-27, Role 5 K 보정용)

- `Piping.face_pressure_terms(t, y, context)`, `ParallelChannelCoreFluid.face_pressure_terms(...)`(유로별 dict), `SGHeliumSide.face_pressure_terms(...)`.
  - 잔차와 같은 `channel_face_terms` 계산(공용 `_face_terms`)으로 면별 `acc·fric·form·grav` [Pa]와 `m`(면), `rho`·`A`(셀)를 준다.
  - Type I 우회 유로는 `dp_cal`도 준다.
  - 잔차에는 영향이 없다 (진단·보정 전용).
- 연결부 K는 builder가 `geometry.junctions.<J>.K_form`(증분)을 인접 유로 면의 경계 K에 더해 넣는다.
  - 배관·SG: `inlet`/`outlet` 튜플
  - 노심: 유로별 `face_K`
- S0 보정값은 `config/calibration/s0_pressure_loss_K.yaml`이다. 변수표 patch에 `K_form_fluid_fwd` 행을 추가했다 (46행).
