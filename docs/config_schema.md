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
