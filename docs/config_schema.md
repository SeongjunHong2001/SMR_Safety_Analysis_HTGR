# 설정 스키마 v3

기본 파일은 입력 **양식**이다. `null`과 빈 목록은 미정이다. 실제 HTGR
builder는 필수값이 빠진 채 실행하면 파일 경로와 키를 포함한 오류를 내야
한다. 현재 `main.py --demo`는 오직 `config/solver.yaml`을 읽는다.

| 파일/객체 | 키 | 단위·shape | 조건 | 소유자 |
| --- | --- | --- | --- | --- |
| `geometry.yaml` Core | `core.axial_cells` | 양의 정수 | 실제 Core 사용 시 필수 | 1/2/3 |
|  | `core.channel_length_m` | m, 양의 스칼라/채널별 | Core 유동 시 필수 | 1/2 |
|  | `core.flow_area_m2` | m², 양의 스칼라/face 배열 | Core 유동 시 필수 | 1/2 |
|  | `core.elevation_profile_m` | m, cell 경계 `N+1` | 중력항 사용 시 필수 | 1/2 |
| `geometry.yaml` SG | `heat_exchanger.axial_cells` | 양의 정수 | SG-2 사용 시 필수 | 1/4 |
|  | `heat_exchanger.heat_transfer_area_m2` | m², 양의 셀별 배열 | SG 열교환 시 필수 | 1/4 |
| `geometry.yaml` Compressor | `compressor.location` | 연결 인스턴스 ID | 실제 계통 연결 시 필수 | 1/2 |
| `operation.yaml` Reactor | `reactor.initial_thermal_power_W` | W, 비음수 | 지정 총출력 mode에서 필수 | 1/3 |
| `operation.yaml` Primary | `primary.initial_pressure_Pa` | 절대 Pa, 양수 | 기준 위치 확인 후 필수 | 1/2/5 |
|  | `primary.initial_temperature_K` | K, 양수 | 실제 헬륨 초기화 시 필수 | 1/2 |
|  | `primary.initial_mass_flow_kg_s` | kg/s, 부호는 연결 방향 | 정상 운전점 선택 시 필수 | 1/2 |
| `operation.yaml` Compressor | `compressor.initial_rotation_rate_rad_s` | rad/s, 비음수 | 축속도 상태 선택 시 필수 | 1/2 |
| `scenarios/steady_state.yaml` | `name`, `conditions` | 문자열, 이름→입력 매핑 | 실제 정상 초기화에서 필수 | 1/5 |
| `scenarios/loss_of_flow.yaml` | `name`, `events[]` | 문자열, 순서 있는 이벤트 | 실제 LOFC에서 필수 | 1/5 |
| `events[]` | `t_event`, `target`, `action`, `value` | s, 인스턴스 ID, 제어 종류, SI 값 | 각 이벤트에 필수; 실측 이력 필요 | 1/2/5 |
| `solver.yaml` | `method`, `t_start`, `t_end`, `dt`, `rtol`, `atol` | DASSL, s,s,s,1,상태 단위 | 데모에는 현재 값 사용; 실제 모델 별도 검증 | 5 |

인스턴스 ID는 중복 없이 각 Component를 지정하고 Registry 키의 접두어가
된다. 배열은 `axial_cells=N`에 따라 cell `(N,)`, face `(N+1,)`로
검증한다. SG/플레넘/배관의 실제 상세 형상·물성 키는 담당자가 자료를
확인해 추가하며, 기본값을 문헌 수치처럼 채우지 않는다.

모델 선택은 변수표의 기존 `mode_nu`, `include_nc`, `mode_friction`,
`mode_he_eos`, `mode_sg_he_crossflow`, `mode_sg_level` 및 새
`mode_core_power`를 사용한다. 이 값은 파라미터/설정이며 DAE 상태가
아니다. 선택 가능한 mode, 필요한 추가 입력, 적용범위는 provider가
실제 구현될 때 검증한다. `u_sg_heat_sink`는 SG-0의 외부 신호이며
자동 기본값이 없다. `m_dot_epsilon`도 `eps_mdot_rel`과 검증된 정격
유량이 있어야만 결정한다.

기존 YAML 키 이름은 그대로 유지한다. 변수표 코드명과 YAML 이름의
대응은 `docs/variable_registry.csv`의 `code_reference` 및 Excel
`GitHub 이름 대조` 시트를 따른다. 실제 builder가 추가되면 위
스키마의 조건부 항목을 해당 mode별로 검증하고 fixture 자동 사용을
막는 시험을 추가한다.

## 실제 builder가 읽는 키 (2026-11-27, Role 5 `system/builder.py`)

builder는 실행 전에 조건의 계산 영역에 필요한 키를 모두 검사한다(`check_inputs`).
`null`이나 누락이 있으면 파일·키·시나리오·조건을 모아 `MissingInputsError`(종료 코드 2)로 멈춘다.
아래 키는 실제 Role 1–4 모델과 builder가 쓰는 것이다. 값의 근거 표시([원문]·[계산]·[도면해석]·[가정])는 각 YAML 주석에 있다.

### `geometry.yaml`

| 키 | 단위·shape | 쓰는 곳 |
| --- | --- | --- |
| `core.elevation_profile_m` | m, face (N+1,) = 18, 유동 방향(위 → 아래) | 노심 유체·고체 공통 격자 |
| `core.channels.<average\|hot\|bypass>.node_id` | 노드 ID 200/210/220 | builder: 세 노드를 Core 하나로 묶음 |
| `core.channels.*.flow_area_m2` 또는 `cell_flow_area_m2` | m², 스칼라 또는 cell (N,) | 노심 유체 (셀 배열이 우선) |
| `core.channels.*.hydraulic_diameter_m` 또는 `cell_hydraulic_diameter_m` | m, 스칼라 또는 (N,) | 마찰·대류, 결합의 젖은 면적 4A/D_h·Δz |
| `core.channels.*.shape` | `circular` \| `parallel_plates` | 층류 f·Re, Nu |
| `core.channels.bypass.type_I_fraction` | 영역별 분율, 합 0.11 | 우회 고정 분율, 결합의 우회 면적 배분 |
| `core.fuel_block`, `core.radial_m`, `core.layout`, `core.block_porosity`, `core.cell_region` | m, 개수 | 노심 고체(단위격자·링), 배럴·riser·RPV 반경 |
| `flow_network.nodes.<이름>.node_id`, `type` | ID, `pipe`\|`branch`\|`pump`\|`valve` | builder: 장치 종류 |
| pipe: `flow_area_m2`, `hydraulic_diameter_m`, `length_m`, `cells`, `z_inlet_m`, `z_outlet_m` | m², m, m, 정수, m | Piping (또는 `z_faces_m` 비균일 셀) |
| branch·valve: `volume_m3`, `z_inlet_m`, `z_outlet_m` (`z_bottom_m`·`z_top_m` 선택) | m³, m | Plenum (valve는 `mlsv: open`일 때 무손실 플레넘) |
| `junctions.<J>.from`, `to`, `K_form` (선택) | 노드 ID, 무차원 증분 | builder: 연결 조립. K_form은 인접 유로 면의 기본 경계 K에 정·역 대칭으로 더함 |
| `heat_structures.HS-06.outer` | 노드 ID | builder: 배럴·RPV 열링크를 붙일 riser |
| `heat_exchanger.axial_cells`, `heat_transfer_area_m2` | 9, m² (9,) | SG-0 |
| `heat_exchanger.helium_side.{fsh,ees}_node_id`, `elevation_profile_m` | ID, m (10,) | builder·SG He 측 셀 고도 |
| `heat_exchanger.helium_side.mean_flow_area_m2`, `hydraulic_diameter_m` | m², m | SG He 측 1D 면적(공극률 × 전면)·D_h. `min_free_flow_area_m2`는 SG-2 상관식용 |
| `rpv.inner_diameter_m`, `outer_diameter_m` | m | RPV 측벽 |
| `rccs.boundary.radius_m`, `emissivity` | m, 1 | RCCS-0 |
| `pressure_loss_calibration.segments.<구간>.{junctions, nodes, design_Pa}` | ID 목록, Pa | K 보정 목표 (`system/calibration.py`) |

### `operation.yaml`

| 키 | 단위 | 쓰는 곳 |
| --- | --- | --- |
| `reactor.power_distribution.{type, axial_fraction, channel_fraction}` | 1, (17,), 유로별 | 지정 출력 분포. 조건의 `reactor.power_distribution`과 같아야 함 |
| `compressor.*` (`rated_heat_to_helium_W` 등) | W, rad/s, Pa | 순환기 입열 Q = Q_rated (ṁ/ṁ_rated)³ |
| `boundaries.sg_heat_sink.{mode, phase_change, water_temperature_K, overall_U_estimate_W_m2K}` | –, K, W/(m²·K) | SG-0 (조건의 물 온도와 같아야 함) |
| `initial_distribution.stations[].{node, T_K, p_Pa}` | K, Pa | 초기 추정 온도, 지점 압력 비교 |
| `helium_inventory.reference_kg` | kg | 재고량 비교 기준 (R22 2,585 kg, 정체 공간 포함 총재고와 비교) |
| `helium_inventory.stagnant_volume_temperature_K` | K | 정체 공간 He 온도 [가정]. 비교용 총재고 추정에만 쓰며 과도 질량 보존 진단에는 쓰지 않음 |

### `scenarios/steady_state.yaml`

- `id: S0`, `default_condition`, `conditions.<이름>`
- **공통:**
  - `domain.{nodes, heat_structures, closed_junctions}`
  - `reactor.{thermal_power_W, power_distribution}`, `bypass_model`
  - `helium_properties` (`파일#helium`), `solid_properties`
  - `rccs.{model_level: RCCS-0, boundary_temperature_K}`
  - `comparison` (참고), `acceptance` (판정)
- **개방 영역 (core_benchmark):**
  - `inlet_boundary.{junction, mass_flow_kg_s, temperature_K}`
  - `outlet_boundary.{junction, pressure_Pa}`
  - `benchmark_subcases`, `default_subcase` (2b 고정 물성만 구현)
- **폐루프 (system_loop):**
  - `pressure_closure.{type: pressure, location, pressure_Pa, replaced_residual_id}`
  - `circulator.{boundary_type: prescribed_mass_flow, mass_flow_kg_s}`, `mlsv: open`
  - `sg_heat_sink.{model_level: SG-0, calibration_target_helium_outlet_K}`
- **`solution`:**
  - `method: solve_steady_state`, `jacobian: sparse_finite_difference`
  - `newton.{max_iter, rtol, atol}`, `acceptance.scaled_residual_max`
  - `hold_transient.{t_end_s, dt_initial_s, rtol, atol, max_step_s, max_order}`
  - 데모용 `solver.yaml`과 별개다.

### 덮어쓰기 파일

- **`config/assumptions/*.yaml` (`--assumptions`):**
  - `meta`·`geometry`·`operation`만 쓴다.
  - config의 `null`/누락 키만 채운다. 원문 값을 바꾸면 오류다.
  - 메타데이터의 `assumptions`에 기록된다. 형식은 `config/assumptions/README.md`.
- **`config/calibration/*.yaml` (`--calibration`, `python main.py --calibrate`로 생성):**
  - `geometry.junctions.<J>.K_form`만 담는다. `meta.input_fingerprint`(유효 입력값 지문)가 현재 입력과 다르면 적용하지 않고 멈춘다.
  - `meta.kind: calibration`, `meta.basis.condition`, `meta.depends_on_assumptions`, 입력 해시를 담는다.
  - 다른 조건에 넘기거나 의존 가정 파일 없이 넘기면 오류다.
