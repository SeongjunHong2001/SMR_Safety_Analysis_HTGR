# Role 1 형상·운전점·검증 자료 인수인계

기준은 `main` `8b98291`(#10)이고 브랜치는 `role1/mhtgr350-inputs`다. 처음에는 `main` #8(`0c79a97`)
스냅샷에서 작업했고, #9(변수표 460행)·#10(조립·솔버)을 반영해 옮겼다. 근거 문서는 담당 1
「[1] 기준 설계와 계통 구성 조사 보고서」 판 1.0(2026-11-26)과 판 1.1 §3.4(M08),
OECD/NEA 벤치마크 R1(NEA/NSC/R(2017)4), INL 결과 보고서 R3·R4다.

## 변경 파일

| 파일 | 내용 |
| --- | --- |
| [config/geometry.yaml](../../config/geometry.yaml) | 노심 3유로 × 17셀, 수력 노드 100~350·400~420·500·900~930, 연결부 J01~J91, 열구조체 HS-01~14, SG(PSID 350관), RPV, RCCS-0, 플랜트 고도(M08), 재료 |
| [config/operation.yaml](../../config/operation.yaml) | 초기값 5개, 정격 압력·온도 분포(지점 10개), SG·RCCS 경계값, He 재고량, 출처별 참고값 |
| [config/scenarios/steady_state.yaml](../../config/scenarios/steady_state.yaml) | S0 조건 2개(`core_benchmark` 기본, `system_loop`), 풀이 방법, 미정 입력 목록 |
| [data/properties/r1_benchmark_properties.yaml](../../data/properties/r1_benchmark_properties.yaml) | R1 Appendix IV 헬륨·흑연·금속·세라믹·단열재 물성식, Ex. 2a·2b 고정 물성 |
| [data/properties/graphite.csv](../../data/properties/graphite.csv), [helium.csv](../../data/properties/helium.csv) | 위 식으로 계산한 비조사 흑연(H-451, 2020)과 헬륨(0.1·5.0·6.39 MPa) 표 |
| [data/validation/benchmark/manifest.yaml](../../data/validation/benchmark/manifest.yaml) | S0 비교 자료(R1 정의, R3 Table IV, R4 그림)의 위치·경계조건 |
| [docs/registry_patches/role1/variable_patch.csv](../registry_patches/role1/variable_patch.csv) | 변수표 변경안 160행 (수정 154, 추가 6). 기준은 460행 마스터 |
| `tests/inputs/test_*.py` 4개 | 입력 검증 테스트 38개 |

`docs/registry_patches/role1/equation_patch.csv`는 만들지 않았다. 새 방정식이 없다.

## 공통 규약

- 단위는 SI다. 압력은 절대압 Pa, 온도는 K, 길이는 m다. 원문 단위는 YAML 주석에 적었다.
- z = 0은 R1 RPV 하부 외면이고 위쪽이 +다. SG 쪽 노드는 음수 z를 가진다(최저 SG 용기 바닥 −18.2 m).
- 배열은 정격 유동 방향(입구 → 출구) 순서다. 노심과 SG He 측은 고도가 감소한다. cell은 (N,), face는 (N+1,)이다.
- 수력 노드의 `z_inlet_m`·`z_outlet_m`은 연결부 고도다. 루프를 따라 이어져 Σ Δz = 0이다. `z_bottom_m`·`z_top_m`은 체적이 차지하는 높이 범위다.
- `null`은 미정이며 0이 아니다. 주석의 `[원문]`·`[도면]`·`[계산]`·`[가정]`·`[제안]`·`[미확보]`가 값의 상태다.
- PyYAML은 `350.0e6`처럼 지수 부호가 없는 표기를 문자열로 읽는다. 큰 값은 정수로 쓰거나 `3.5e+8`처럼 부호를 붙인다. 테스트가 이 실수를 잡는다.

## 담당자별 입력

### 담당 2 — 헬륨 유동, 순환기, 배관, 플레넘

| 키 | 단위·shape | 내용 |
| --- | --- | --- |
| `geometry: core.elevation_profile_m` | m, (18,) | 노심 face 고도 13.4339 → 1.9356 |
| `geometry: core.channels.{average,hot,bypass}` | | 유로면적, `cell_flow_area_m2` (17,), 수력직경, 체적. 우회는 Type I 분율 7개와 제어봉 채널 반경 |
| `geometry: flow_network.nodes.<인스턴스>` | | `node_id`, `cells`, 면적, 수력직경, 길이, 체적, 입출구 고도 |
| `geometry: junctions`, `pressure_loss_calibration` | Pa | 연결부(from/to 노드)와 PSID 손실 배분 9구간 (합 90.8 kPa). `K_form`은 S0에서 보정 |
| `operation: primary.*` | Pa, K, kg/s | 초기 압력 6337300 Pa(하부 플레넘), 노심 입구 532.15 K, 157.1 kg/s |
| `operation: initial_distribution` | | 지점별 p·T·ρ 초기 추정값. 정수압 ρgΔz는 빠져 있다 (루프 합 약 −0.4 kPa) |
| `operation: compressor.*` | rad/s, Pa, W | 649.26 rad/s, 승압 91 kPa, 입열 3.1 MW. 관성·성능곡선 미확보(M09) |
| `data/properties: helium` | | R1 Table AIV.13 식 (이상기체 R = 2077.46, c_p 5195) |

`p_sys`는 `lower_plenum`(노드 230) 인스턴스가 등록한다. 미확보 체적(110·300·330·340·350)은 `null`이다.

### 담당 3 — 노심 고체, 출력, 열전달

| 키 | 내용 |
| --- | --- |
| `geometry: core.cell_region`, `r1_reporting_layer` | 셀 영역(연료 10셀)과 R1 Fig. II.7 보고 층 번호 |
| `geometry: core.layout`, `fuel_block`, `radial_m`, `block_porosity` | 기둥·블록 수, 연료공·냉각공·콤팩트 치수, 간극 0.1 mm, 등가반경, 공극률 |
| `geometry: heat_structures` HS-01~HS-08 | 재료, 반경·고도 범위, 경계 노드 |
| `data/properties: graphite_H451`, `graphite_2020`, `fixed_properties_ex2ab` | 온도·조사량 의존 식과 Ex. 2b 고정값 |
| `operation: reactor.*` | 350 MW, `power_mode: specified_total`, 출력분포 균일 (`power_distribution`) |

### 담당 4 — SG, RPV, RCCS

| 키 | 내용 |
| --- | --- |
| `geometry: heat_exchanger.*` | PSID Amd.3–6: 350관, 22.225 × 3.302 mm, 164.5 m, 4020 m². He 셀 9개(FSH 1 + EES 8)의 면적·높이·고도, 2차측 대향류 매핑 |
| `geometry: heat_exchanger.vessel`, `plant_elevations` | SG 용기 z −18.2 ~ 8.2, SGS·shroud 고도, 열중심 고도차 15.5 m |
| `geometry: rpv`, `rccs`, HS-06·07·09~14 | RPV 치수·재질, RCCS-0 경계(r 4.633 m, z 18.6813 m, ε 0.80), RCCS-1 참고값 |
| `operation: boundaries.sg_heat_sink` | SG-0 = 일정 온도 액체 물 열침원 (물 466.15 K [제안], 5,000 kg). U는 S0에서 보정, 목표 He 출구 528.15 K·제거열 352.2 MW |
| `operation: boundaries.sg_secondary`, `rccs` | SG-2용 급수 137.4 kg/s·466.15 K·20.7 MPa, 주증기 814.15 K·17.3 MPa. RCCS 30 °C |

### 담당 5 — 조립, 솔버, 시나리오 실행

- `steady_state.yaml`의 `conditions`는 이름 → 조건 세트다. 기본은 `default_condition: core_benchmark`이고, 계통 루프 계산은 `system_loop`를 고른다.
- 각 조건의 `domain`이 계산 노드·열구조체·닫는 연결부를 정한다. 테스트가 경계 누수와 루프 닫힘을 검사한다.
- `system_loop`의 정상 초기화 제약은 기준압이다(`p_sys = 6337300 Pa`). 교체할 질량식 ID는 담당 2와 정한다.
- 새 키 설명을 `docs/config_schema.md`에 옮겨야 한다. 이 문서는 담당 5 소유라 고치지 않았다.

## 결정·제안 사항

| 항목 | 내용 | 상태 |
| --- | --- | --- |
| SG 설계 기준 | PSID Amd.3–6 (350관). 보완본[4]의 439관과 섞지 않음. 변수표 관 수·길이·피치·면적·내면적·벽 체적 모두 PSID (아래 "SG 전열관 기준 변경" 절) | 결정 (2026-11-26, 2026-11-27 재확인) |
| 배열 순서 | 정격 유동 방향(입구 → 출구) | 결정 |
| 노심 출력분포 | 균일 (연료 10층 × 0.1, 기둥당 1/66, 고온채널 첨두계수 1.0) | 결정 (가정) |
| 비연료 블록 유로면적 | 상부반사체 = 연료 냉각공, 하부반사체 R1 0.199/0.177, 전이 0.1941 · 유량분배 0.1885 · 포스트 0.1883 (PSID 유로 치수) | 결정 |
| 압력 기준 (보고서 D2) | PSID 절대압 분포 (하부 플레넘 6.3373 MPa). 벤치마크 비교 조건만 R1 출구 6.39 MPa 경계 | 결정 |
| 6.39 MPa 해석 (M06) | 절대압, 노심 출구 (R1·R4와 같음) | 제안 유지 |
| 헬륨 물성 세트 | R1 값으로 전 계통 통일 (R = 2077.46, c_p = 5195, c_v = 3117.54) | 결정 |
| SG-0 열침원 | 2차측 = 일정량(5.0 t)·일정 온도 액체 물. Q_i = U·A_i·(T_He,i − T_물), 상변화·2차측 유동 없음. U는 S0에서 He 출구 255 °C가 되도록 보정 (추정 421.5 W/(m²·K)) | 결정. 물 온도 193 °C는 제안 |
| 순환기 (S0 루프) | 성능곡선이 없어 유량 157.1 kg/s 지정, 승압은 결과 | 결정 |
| FSH–EES 간격 1.83 m | FSH 셀(310)에 배정 | 결정 |
| SG 셀 면적·관 길이 | 번들 높이 비례 배분 | 결정 |
| S0 조건 구조 | `core_benchmark`(기본, Ex. 2b부터) / `system_loop` | 확인 필요 |
| cross duct 길이, He 재고량 기준, RCCS 경계 온도 | 현재 6.5 m, 2,585 kg, 30 °C | 선택 필요 |

## SG 전열관 기준 변경 (439관 → PSID 350관)

**요약:** 변수표의 SG 전열관 세트를 보완본[4]의 439관(Combustion Engineering 1987 설계서)에서
PSID Amd.3–6의 350관으로 바꾼다. 관 외경·두께는 같고 관 수·관 길이·축방향 피치가 다르다.
관 수만 바꾸면 판본이 섞이므로 같은 판본의 관 세트 전체를 함께 바꾼다.

### 두 판본 비교

| 항목 | 439관 (현재 마스터) | 350관 (변경안) |
| --- | --- | --- |
| 출처 | 보완본[4] p.4 ← CE (1987) DOE-HTGR-86-129 Rev.1 / HFD-42102 Rev.2, Table 2.3-4 | 보고서[1] §3.2 ← R19 PSID (DOE-HTGR-86-024) Amd.3–6, Table 5.3-8 |
| 관 수 | 439 | 350 |
| 관 외경 / 두께 | 22.225 / 3.302 mm | 22.225 / 3.302 mm (같음) |
| 평균 유효 관 길이 | 131.363 m | 164.5 m |
| 횡 / 종 피치 | 38.1 / 35.331 mm | 38.1 / 36.6 mm |
| 전열면적 (외면) | 4,026.54 m² | 4,020 m² |
| 관군 높이 | 표의 구간 합 26.70 ft ≠ 표 총계 27.70 ft, FSH/EES 전이부 관 길이 TBD (보완본[4] p.24) | EES 7.86 m, FSH 0.861 m, 사이 간격 1.83 m |
| 정격점 | He 157.2 kg/s, 685.6 → 254.8 °C, 352.1 MW | He 157.3 kg/s, 686 → 255 °C, 352.2 MW |

### 바꾸는 이유

1. **계통 전체와 같은 판본:** 보고서[1]은 해석 형상을 "벤치마크 기하(설계 개정 1986–88, PSID Amd. 3–6)"로 한정한다(§6.1).
   노심, 압력 분포, 순환기, RCCS, SG 고도(M08)가 모두 이 판본이다. SG 관만 439관으로 두면 관 세트와
   번들 높이·SG 셀 고도·운전점이 서로 다른 문서에서 온다.
2. **축방향 형상:** SG He 측 9개 셀의 높이·고도·셀별 면적은 PSID 번들 높이와 Fig. 5.3-10 고도로 만들었다.
   439관 설계서는 관군 높이 합이 표 총계와 1 ft 어긋나고 전이부 관 길이가 TBD여서, 보완본[4]도
   "축방향 길이 단정 금지"로 적었다. 이 자료로는 SG 셀 고도를 만들 수 없다.
3. **팀 결정:** SG 형상 기준 PSID (2026-11-26), 관 수 포함 재확인 (2026-11-27).

**판본 선후에 대해:** 보완본[4]와 `main` #9 비고는 439관을 350관보다 뒤의 설계로 본다. 이 변경은
"더 늦은 판본"이 아니라 "계통 전체와 같은 판본"을 기준으로 고른 것이다. 나중에 439관 판으로 계산하려면
관 세트뿐 아니라 번들 높이·SG 고도·운전점을 모두 그 판본으로 함께 바꿔야 한다.

### 영향

- **SG-0 (현재 열침원 모델):** 거의 없다. 전열면적이 0.16 % 작아지지만 U를 S0에서 He 출구 255 °C에 맞춰 보정한다.
- **SG-2 (상세 모델):** 관 길이 +25 %, 관내 질량유속 +26 % (2,048 vs 1,631 kg/(m²·s)), 관당 유량 0.393 vs 0.313 kg/s,
  축방향 피치 +3.6 %. 2차측 압력강하·비등 시작 위치와 He 교차류 상관식 입력이 달라진다.
- **439관 판에만 있는 값:** 관군 He 마찰 9.651 kPa, 관내 물 마찰 0.946 MPa, 급수 오리피스(20개 직렬, 6.5 mm)는
  값을 남기되 PSID 형상과 함께 쓰지 않는다. PSID 대응값은 아직 없다.

### 변경 항목

| 변수 | 마스터 (439관) | 변경안 (350관) |
| --- | --- | --- |
| `N_sg_tubes` | 439 | 350 |
| `L_sg_tube_effective` | 131.363 m | 164.5 m |
| `L_sg_tube_fsh` / `L_sg_tube_ees` | 13.308 / 118.055 m | 16.241 / 148.259 m (번들 높이 비례 배분) |
| `S_sg_tube_longitudinal` | 0.035331 m | 0.0366 m |
| `A_sg_tube_outer_total` | 4026.539 m² | 4020.0 m² |
| `A_sg_tube_outer_fsh` / `A_sg_tube_outer_ees` | 407.882 / 3618.658 m² | 396.883 / 3623.117 m² (번들 높이 비례 배분) |
| `A_sg_tube_inner_total` | 2830.06 m² | 2825.48 m² (350 × π × d_i × 164.5) |
| `V_sg_tube_wall_total` | 11.32 m³ | 11.302 m³ (350 × π/4 × (d_o² − d_i²) × 164.5) |
| `D_sg_tube_outer`, `t_sg_tube_wall`, `D_sg_tube_inner`, `S_sg_tube_transverse` | 값 같음 | 출처만 PSID로 |
| `*_psid_ref` 4행 (`N_sg_tubes_psid_ref` 등) | PSID 값을 따로 보관 | 값 유지, 상태를 "기본 행에 채택 — 판본 대조용"으로 |

기본 행 비고의 "PSID 값은 `*_psid_ref`로 분리" 문구는 "439관 값은 보완본[4] p.4에 보존"으로 바꾼다.
`config/geometry.yaml`의 `heat_exchanger.tubes`는 이미 350관 세트이고, 입력 테스트가 관 수 350과
내경·내면적·벽 체적의 계산 일치를 검사한다.

## 미정 입력

- **SG-0 물 온도:** 제안값 193 °C(급수 온도)를 확인해야 한다.
- **블록 유로:** CRE(금속 plenum element)의 유로 치수와 RSC 기둥 하부 블록의 유로는 자료가 없어 각각 연료·표준 기둥과 같게 두었다.
- **조사량 맵:** R1 블록별 조사량 맵(Ex. 2c 가변 물성용)은 벤치마크 참가자 한정이다.
- **체적과 형상:** 미확보 체적(M10), SG He 측 최소 유로면적·환형부 단면적(M11), cross duct 길이(M07), riser 치수 해석(M04)이 남아 있다.
- **비교 자료:** Ex. 2 전체 참가자 비교 결과는 아직 없다. 현재 INL R3·R4만 확보했다.
- **순환기:** 관성과 성능곡선(M09)이 없다.

## 변수표 patch 요약

`variable_patch.csv`의 `old_kind`·`old_value`·단위·shape는 마스터 CSV에서 자동으로 옮겼다.
`reason`이 `[initial]`로 시작하는 행은 `initial` 칸, 나머지는 `reference_value` 칸에 대한 변경이다.

전체 160행 (수정 154, 추가 6)이고 기준은 `main` #10의 460행 마스터다. 병합을 보류할 행은 없다 (D2 결정으로 `p_sys` 보류 해제).

- **SG (수정 36행):** 마스터의 보완본 Rev.1 값을 PSID 값으로 바꾼다. 관 수 439 → 350, 유효 길이, 피치, 외면적과 함께 Rev.1 관 세트로 계산된 내면적(2830.06 → 2825.48 m²)·벽 체적(11.32 → 11.302 m³)도 PSID 관 세트로 다시 계산했다. Rev.1에만 있는 값(관군 마찰, 급수 오리피스 등)은 값을 유지하고 혼용 금지 비고만 단다. `UA_sg`를 SG-0 열침원의 보정 대상으로 연결한다.
- **노심·유체 형상, 운전점 (수정 29행):** 새 config 키로 코드 참조를 잇고 값을 채운다. 셀별 유로면적에 PSID 블록 값을 반영했다.
- **R1 물성·간극·출력분포 (수정 11행):** 물성 파일로 연결하고 헬륨 세트를 R1 값으로 통일한다(`R_he` 2077.46, `c_v_he` 3117.54). 출력분포는 균일 분포다.
- **추가 6행:** FSH–EES 간격, 열중심 고도차, 벤치마크 출구 압력 경계, 블록 공극률, SG-0 물 온도·물 양.

### `main` #9 대조 결과 (수정 78행)

`main` #9(`84caa8f`)는 같은 보고서[1]의 기준설계 데이터표를 변수표에 넣었다 (신규 88행, 기존 40행 보완).
role1 patch를 이 460행 기준으로 다시 만들면서 다음을 처리했다.

- **중복 통합 (4행):** role1이 추가하려던 행 중 #9 신규 행과 같은 것은 해당 행의 update로 바꿨다.
  `f_core_bypass_region` → `f_core_bypass_type_i_region`, `dp_loss_target_segment` → `dp_loop_segment_design_ref`,
  `T_core_he_inlet_bc` → `T_core_he_inlet_benchmark_ref`, `Q_dot_circ_he_rated` → `Q_dot_circ`.
  #9가 값 칸을 비우고 비고에 값을 적은 행(우회 분율, 구간 손실)은 그 방식을 따랐다.
- **SG 판본 (5행):** 팀 결정에 따라 기본 행(`N_sg_tubes` 등)은 PSID 값으로 바꾸고, #9의 `*_psid_ref` 4행은
  값을 두고 상태만 "기본 행에 채택 — 판본 대조용"으로 바꾼다. 기본 행 비고의 "PSID 값은 `*_psid_ref`로 분리" 문구도
  "Rev.1(439관) 값은 보완본[4] p.4에 보존"으로 바꾼다. `H_sg_vessel_psid_ref`의 "SG–노심 고도차 미확보"는 M08로 확보했다.
- **config 키 연결 (69행):** #9 기준값 중 config에 같은 값이 있는 64행은 값은 그대로 두고 코드 참조만 잇는다.
  patch 생성 스크립트가 config 값과 마스터 기준값이 다르면 멈추도록 해서 모두 일치를 확인했다 (상대 오차 2×10⁻⁴ 이내).
  값 칸이 빈 5행(`p_loop_node_design_ref`, `V_he_loop_free`, `I_circ_rotor`, `p_psrv_reseat`, `mu_he`)도 코드 참조만 잇는다.
- **#9 비고와 다른 점:** #9는 `N_core_axial`·`L_core_channel` 비고에 "config는 확정 전 null 유지"라고 적었지만
  이 브랜치의 `geometry.yaml`은 담당 1 제안값(17셀, 11.4983 m)을 채웠다. 격자 확정은 담당 2·3과 협의한다.
- **연결하지 않은 #9 행:** 사고 시나리오 행(`t_scram`, `dt_flow_ramp`, `p_sys_scenario_target` 등)은 `loss_of_flow.yaml` 작업 때,
  판정기준 행은 `T_rpv_limit_normal`만 S0에 있고 나머지는 사고 시나리오 때 연결한다. `H_rccs_chimney_psid_ref`(14 m, R23 29 m와 상충)는
  config에서 `null`로 둔다.

## 시험

```bash
python -m pytest -q                      # 103 passed (main 65 + 입력 검증 38)
python -m pytest -q tests/inputs         # 형상 14, 운전점 9, S0 8, 물성 7
python scripts/compare_variable_registry.py   # CSV=460, Excel=460, 불일치 0
python main.py --demo                    # 최종 상태 [316.19074506, 306.64335141, 47.73696821]
```

입력 테스트는 배열 길이, 양의 기하, 노드·연결부 참조, 루프 고도 닫힘, 파일 간 값 일치,
보고서 검산(열수지·압력 루프·재고량 밀도), 물성 표와 식의 일치를 확인한다.
실제 HTGR 물리 정확도 시험은 아니다.

## 다음 연결

1. 담당 2·3·4는 위 키로 G1 단독 모델을 만든다. S0 `core_benchmark`의 Ex. 2b(고정 물성)부터 시작하면 참가자 한정 조사량 맵 없이 진행할 수 있다.
2. 담당 5는 `variable_patch.csv`를 검토해 마스터 CSV·Excel에 병합한다. 그리고 `config_schema.md`에 새 키를 반영한다.
3. 남은 선택 사항(S0 조건 구조, SG-0 물 온도, cross duct 길이, He 재고량 기준, RCCS 경계 온도)을 정한다.
4. 담당 1은 남은 형상 자료(M10·M11·M07·M04)를 확보한다. 사고 시나리오(`loss_of_flow.yaml`)는 이후 작업한다.
