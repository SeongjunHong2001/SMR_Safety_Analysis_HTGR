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

## STEP B — 정상상태 builder와 실행기 (2026-11-27)

### 실행

```bash
python main.py --demo                                         # 기존 3변수 데모 (변경 없음)
python main.py --scenario steady_state                        # 기본 조건 core_benchmark. 가정 파일 없으면 미정 입력 목록, 종료 코드 2
python main.py --scenario steady_state --condition system_loop \
    --assumptions config/assumptions/provisional_s0.yaml      # 정상해 → 검증 → 무외란 유지 → output/results/
```

출력은 `output/results/steady_state_<condition>_{states.csv, outputs.csv, meta.json}`이다.
- `states.csv`: `time_s` + Registry 이름/shape 평탄화 열 (`name[i]`). 행 = 무외란 유지 적분의 수용 시각 (첫 행 = 정상해).
- `outputs.csv`: Component `outputs()` 평탄화.
- `meta.json`: 시나리오·조건·모드, 입력 파일과 SHA-256, 적용한 [가정] 키·값(영역 안/밖 구분), builder 모델 가정,
  코드 버전(commit·branch·dirty), 상태 단위(SI), 출력 단위(이름에서 추정), 솔버 통계(Newton 반복·잔차·최악 잔차 ID,
  DASSL 수용/거절 스텝·차수·Newton 실패), 보존 진단, comparison(참고)·acceptance(판정).

### 구성

| 파일 | 역할 |
| --- | --- |
| `system/builder.py` | `load_inputs(scenario, condition, assumptions=)` → `ScenarioInputs`; `check_inputs` (미정 입력 전체 목록); `build_system(inputs, mode="steady"/"transient", U_sg=, power=)` → `BuiltSystem` (system, parts, inventory, flow_order, connections, 가정 기록); `validate_built` |
| `system/boundaries.py` | 개방 영역 경계 `FixedSystemPressure`(ṗ + (p − p_ref)/τ = 0), `PortBoundary`(유량 또는 압력 + 유입 엔탈피) — 원형에서 이전 |
| `system/simulation.py` | `settings_from_scenario`(시나리오 `solution` 절), `solve_steady`, `verify_steady`, `run_hold`, `summary`, `compare`, `flatten_names`, `write_outputs`, `run_scenario` |
| `validation/conservation.py` | `helium_mass`, `helium_energy`, `solid_energy`, `boundary_rates`, `steady_energy_balance`, `transient_balance` |
| `config/assumptions/provisional_s0.yaml` | null 입력의 [가정]값 (우회 D_h 4 mm 평행평판, SG He D_h 0.05 m, M10 체적, 340 유로) |
| `config/scenarios/steady_state.yaml` `solution` | Newton·희소 Jacobian·허용 잔차·무외란 유지 설정. `pressure_closure.replaced_residual_id: lower_plenum.mass` |

### 조립 규칙

- **노드 → 장치:** `core.channels.*.node_id`(200·210·220) → `Core` 하나, SG He 측 노드(310·320) → `HeatExchanger` 하나.
  그 밖의 `flow_network.nodes`는 type으로 고른다: pipe → `Piping`, branch → `Plenum`, pump → `Compressor`(유량 지정),
  valve → `Plenum`(`mlsv: open`일 때만, 무손실 [가정]).
- **유체 연결:** `geometry.junctions` 중 양 끝이 영역 안이고 `closed_junctions`가 아닌 것을 쓴다.
  - 같은 장치 쌍의 병렬 연결은 하나로 묶는다 (J06–J08 → `J06_J08`).
  - 장치 내부 연결(J15)은 생략한다.
  - 영역을 자르면서 경계로 지정되지 않은 연결부가 있으면 오류다.
- **경계 손실:** 배관 끝이 플레넘에 붙으면 급축소 0.5 / 급확대 1.0 + 베르누이(`("plenum", K, K_rev)`)를 쓴다.
  그 밖의 연결부 K는 미보정(0)이다.
- **p_sys:**
  - 폐루프는 `pressure_closure.location` 플레넘이 단일 소유한다. 정상 모드에서는 그 플레넘의 질량식 1개를 p_sys − p0로 바꾼다.
  - 개방 영역은 `FixedSystemPressure`를 쓴다.
  - 조립 후 검사: 등록 1회, differential, owner `primary_loop`, 폐합 잔차 수, 과도 질량식 유지, 모든 유체 포트 연결 1회.
- **열경계:** HS-06 바깥 노드(riser)에 배럴→riser, riser→RPV(He 간극), RPV→RCCS-0를 연결한다.
  영역 안 미구현 열구조체(HS-09·10·11)는 단열로 두고 메타데이터에 남긴다.
- **[가정] 파일:** null인 키만 채운다. 원문 값을 바꾸려 하면 오류다. 적용 항목 중 영역 안 장치에 해당하는 것만 "사용"으로 기록한다.

### 정상 초기화 순서 (steady_state.yaml `solution`)

1. **정상 풀이:** mode="steady" 계를 ydot = 0 Newton으로 푼다 (max 30회, 희소 FD).
   - 패턴은 초기 추정값에서 y와 ẏ를 한 열에 함께 섭동해 한 번 만든다.
2. **검증:** mode="transient" 계를 보정 U 고정으로 만들고, 정상해를 이름별로 옮겨 ydot = 0 잔차를 다시 계산한다 (≤ 1e-6).
3. **유지 적분:** 과도 계를 일치 초기화한 뒤 DASSL로 10 s 무외란 유지 적분한다.
   - 표류(상태 척도)와 전역 질량·에너지 수지를 기록한다.
   - `consistent_initialize`의 Newton 미지수는 [대수 y, 미분 ẏ] 순서라, 패턴 열을 `z_order`로 바꿔 넘긴다.
4. **보조 풀이:** pseudo-transient fallback은 미구현이다. 정상 폐합·U 보정 계는 지수 2라 그대로 적분할 수 없다.

### S0 결과 (2026-11-27, 가정 입력 포함)

| | core_benchmark | system_loop |
| --- | --- | --- |
| 미지수 | 347 | 481 / 과도 480 |
| 정상해 (Newton, 시간) | 4회, 잔차 5.3e-11, 3 s | 4회, 잔차 1.0e-7, 6 s |
| 무외란 10 s | 12 스텝, 표류 9e-14, 질량 2e-15, 에너지 2e-12 | 12 스텝, 표류 3.6e-8, 질량 4e-16, 에너지 7e-10 |
| 노심 출구 / 연료 최고 | 960.5 K / 1131.7 K | 960.3 K / 1131.5 K |
| 노심 Δp / 배럴 평균 | 23.1 kPa / 568 K | 23.3 kPa / 568 K |
| RPV 최고 / RCCS | 463.2 K / 428 kW | 463.0 K / 427 kW |
| SG 제거열 / U / 순환기 승압 | — | 352.67 MW / 474 W/(m²·K) / 36.8 kPa |
| He 재고량 | 808 kg | 1,762 kg |

출력 +10 % 램프(0.5 s) 3 s 과도에서 질량 오차는 2.5e-10이다. 에너지 오차는 램프 꺾임에서 최대 7.8e-4(사다리꼴 적분)이고 끝값은 7.8e-5다.

## 압력손실 K 보정 (2026-11-27)

```bash
python main.py --calibrate --assumptions config/assumptions/provisional_s0.yaml      # → config/calibration/s0_pressure_loss_K.yaml
python main.py --scenario steady_state --condition system_loop     --assumptions config/assumptions/provisional_s0.yaml --calibration config/calibration/s0_pressure_loss_K.yaml
```

- **목표:** `geometry.pressure_loss_calibration.segments.*.design_Pa`(PSID R19 Table 5.3-5)를 구간별 비가역 손실(마찰 + 형상손실)로 본다.
  - 정수압과 가역 동압 변화(플레넘 베르누이, 열가속)는 목표에 넣지 않는다 (geometry 주석의 [가정]).
- **K 위치:** 구간 연결부 목록에서 면을 실을 수 있는 첫 연결부(`builder.junction_face`)에 둔다.
  - 하류 유로 입구 면을 먼저 쓰고, 없으면 상류 유로 출구 면을 쓴다. SG 내부 J15와 순환기 쪽 J19·J01은 면이 없다.
  - 노심 묶음은 연료 유로 J06·J07에 같은 K를 준다. Type I 고정 분율 우회 J08은 분율식이 유량을 정하므로 뺀다.
- **K_form 의미:** 면의 기본 경계 K(플레넘 급축소 0.5·급확대 1.0, 그 밖 0)에 더하는 증분이다.
  - 음수도 허용하지만 면 합은 0 이상이어야 한다. 예: J02 −0.09 → 입구 K 0.41.
  - 정·역 대칭으로 더한다 [가정] (변수표 `K_form_fluid_rev` 주석과 같다). LOFC 역류 전에 다시 본다.
- **갱신 (식 7):** K ← K + (Δp_target − Δp_irr)/q, q = ṁ|ṁ|/(2ρA²).
  - 저마하수라 밀도가 국부압과 무관하고 유량은 지정이라, 비노심 구간은 한 번에 맞는다.
  - 노심은 평균·고온 유로 분배 때문에 한 번 더 반복한다. 실제로 3회 반복에 걸렸다(정상해 약 6.4 s씩).
- **손실 배정:**
  - 유로 면 마찰과 내부 면 형상손실 → 노드 구간 (목록에 없는 120·340은 출구 연결부 구간)
  - 경계 면 형상손실 → 그 연결부 구간
  - 노심은 평균 유로 경로를 쓴다.
- **보정 파일 규칙:**
  - `meta.kind: calibration`이다.
  - 다른 조건에 넘기면 오류다. `basis.condition`이 그 조건과 다르면 막는다. 예: core_benchmark에 노심 K가 들어가면 R1/R3 비교가 깨진다.
  - `depends_on_assumptions` 파일을 함께 넘기지 않아도 오류다.
  - 메타데이터에는 가정과 따로 `calibration` 항목으로 남는다.
- **회귀 시험:** `tests/integration/test_calibration.py`는 현재 입력으로 다시 보정한 K가 저장 파일과 같은지 본다.
  입력이 바뀌어 실패하면 `--calibrate`로 파일을 다시 만든다.

현재 값 (2026-11-27 #31 형상 보완 뒤 재보정. 괄호 안은 #30 가정 입력 때 값):

| 구간 | 목표 Pa | K 연결부 | K_form 증분 | 비고 |
| --- | --- | --- | --- | --- |
| mc_outlet_to_core_inlet | 3,400 | J02 | −0.090 (−0.090) | 기본 플레넘 경계 손실만으로 목표 초과 → 입구 K 0.41 |
| riser | 4,800 | J05 | +1.654 (+1.649) | |
| core | 34,500 | J06, J07 | +14.81 (+12.11) | 하부 블록 셀별 D_h(0.068–0.116 m) 반영으로 마찰이 줄어 K 증가 |
| core_outlet_plenum | 15,800 | J12 | +4.529 (+4.529) | hot duct 입구 면 |
| hot_duct | 4,800 | J13 | +0.477 (+0.477) | |
| sg_inlet_plenum | 2,100 | J14 | +26.05 (+50.47) | SG 평균 유로면적 6.66 m² 동압 기준 |
| sg_bundle | 20,000 | J16 | +409.2 (+789.7) | 번들 마찰(D_h 0.0577) 외 나머지를 출구 면에 몰아 둠 → SG-2 교차류 상관식(A_min 3.85 m²)이 대체 |
| sg_outlet | 3,400 | J17 | +9.383 (+5.865) | 340 단면 2.49 m² (가정 2.0) |
| mc_inlet | 2,000 | J18 | +4.225 (+1.403) | |

**결과 (system_loop, 가정 2개 + 보정 K):**
- **순환기 승압:** 91.93 kPa = 비가역 90.80 + 가역 순합 1.13 kPa. PSID 91.0 kPa보다 1.0 % 높다.
- **지점 압력:** PSID 대비 3.5 kPa 이내다 (시험). 지점 3과 지점 8의 차이는 지점 정의 때문이다.
- **노심:** Δp 34.7 kPa, 연료 최고 1119.2 K, 배럴 평균 559 K, RPV 최고 462.5 K, RCCS 425 kW.
- **He 재고량:** 유로 2,314.9 kg + 정체 공간 243.6 kg = 2,558.5 kg. R22 2,585 kg 대비 −1.0 %다.
- **유지 적분:** 10 s 표류 2.3e-8, 질량 오차 4e-16, 에너지 오차 5e-10.
- **core_benchmark:** 가정 파일 없이 돈다. 노심 Δp 20.7 kPa로 R3 19.95 kPa와 +3.6 % 차이다 (보정 전 23.1 kPa).
  연료 최고는 1124.7 K다.

### 남은 일 (Role 5)

- **K 보정 재실행:** M10·M11·우회 D_h가 원문 값으로 바뀌면 다시 보정한다.
  SG 번들은 SG-2 교차류 상관식이 들어오면 K 대신 분포 마찰로 바꾼다.
- **다포트 분기:** SCS 운전(J20–J23)에는 포트별 고도가 필요하다. 지금은 장치마다 입구 1, 출구 1이다.
- **STEP C:** `system/scenario.py`, LOFC.
- **변수표 patch 병합:** role1–4 patch를 마스터에 넣는다.
- **pseudo-transient fallback:** 아직 없다.
