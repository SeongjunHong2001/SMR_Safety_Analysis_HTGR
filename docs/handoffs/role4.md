# Role 4 SG·압력용기·RCCS 인수인계

브랜치는 `role4/vessel-rccs-s0`이다. RPV·RCCS-0·S0 core_benchmark는 `8530b9b`(Role 3 노심 고체 `8a44cf7` 위)에 있고,
SG-0과 S0 system_loop는 그 위에서 커밋 전이다.
병합 순서는 Role 1 → Role 2 물성 → Role 3 열전달 → Role 2 유동·장치 → Role 3 노심 고체 → Role 4다.

## 진행 상태

| 항목 | 파일 | 상태 |
| --- | --- | --- |
| RPV 측벽 | `models/vessel/vessel.py` `Vessel` | 구현 (아래) |
| 표면–유체 대류 링크 | `models/vessel/vessel.py` `SurfaceFluidConvection`, `helium_gap_conductance`, `axial_overlap_matrix` | 구현 |
| RCCS-0 고정 온도 열침원 | `models/rccs/rccs.py` `RCCS0Sink` | 구현 |
| RCCS-1 공기 순환 | `models/rccs/rccs.py` `RCCS` | 미구현 (패널·굴뚝 입력 상충, chimney_height_m null) |
| SG-0 (HeatExchanger 자식) | `models/heat_exchanger/{helium_side,wall,secondary_side,coupling}.py` | 구현 (아래). SG-1/2(관벽 축열, 2차측 물·증기)는 미구현 |
| 상·하부 헤드 (HS-09 MCSS·하부헤드, HS-10 UPTPS·상부헤드) | — | 미구현 (단열) |

## RPV `Vessel`

`Vessel(geometry["rpv"], R1SolidProperties(), name="rpv", z_faces=..., initial=dict(T=...), links=(("riser_rpv", "to_a"), ("rpv_rccs", "to_a")))`
- **상태:** `{name}.T` (N,) differential. RPV 측벽을 반경 1노드 × 축 N셀로 둔다. 두께 133 mm, k 40이라 Biot가 작다.
- **포트:** `inner_surface[i]`, `outer_surface[i]`는 같은 키를 가리키는 별칭이다 (단일 벽 노드, 저장항 중복 없음).
- **에너지식:** ρ c_p V dT/dt − Σ 링크 열률 − 축방향 전도 = 0. 열은 `context.heat_rates["링크.to_a|to_b"]` (N,)로만 받는다.
  링크 열률이 없으면 `ModelInputMissing`.
- **물성:** `R1SolidProperties`의 `pressure_vessel` (SA533B, R1 Table AIV.12: k 40, ρ 7500, c_p 580, ε 0.85).

## RCCS-0 `RCCS0Sink`

`RCCS0Sink(geometry["rccs"], vessel, R1SolidProperties(), T_sink=operation["boundaries"]["rccs"]["boundary_temperature_K"])`
- **열 경로:** RPV 외면 → 원자로 공동 → RCCS 경계(r 4.633 m, ε 0.80, 30 °C, R1 Table I.10).
  - 복사 ε_eff σ A (T⁴ − T_s⁴), ε_eff = 1/(1/ε_RPV + (r_o/r_s)(1/ε_RCCS − 1)) = 0.735 (변수표 0.74)
  - 정체 공기 전도 2π k_air Δz (T − T_s)/ln(r_s/r_o) (k_air 0.068)
- **링크:** `.link()`가 ThermalLink `rpv_rccs`를 돌려준다 (to_a = −Q는 RPV로, to_b = [ΣQ]는 열침원으로). 독립 상태는 없다.

## 표면–유체 대류 링크 `SurfaceFluidConvection`

- **역할:** 고체 표면 온도 배열(side a)과 `Piping` 유체(side b)를 잇는 ThermalLink 전달 함수다.
- **격자 대응:** 축방향 격자가 달라도 `axial_overlap_matrix`(고도 겹침 비율 W)로 대응시킨다.
  유체 쪽 온도 T_f,i = Σ W_ij T_j, 열 배분 = Wᵀ Q라서 보존된다.
- **대류 계산:**
  - 대류는 `McEligotConvection`이고 벽온도를 두 번 평가한다. 고체 쪽 직렬 컨덕턴스(`series`)를 받을 수 있다.
  - `helium_gap_conductance`는 동심 원통 헬륨 간극(전도 + 선형화 복사)이다.
- **`Piping` 추가 (Role 2 파일):** `z_faces`와 `coupling_state()`를 넣었다. 키는 노심 유체와 같다.

## S0 core_benchmark 첫 계산

steady_state.yaml `conditions.core_benchmark`의 영역과 경계를 실제 Component로 조립했다. 처음에는 시험 파일에서 손으로 조립했고,
지금은 Role 5 builder(`system/builder.py`)가 config와 연결부 표에서 같은 조립을 한다 (원형 `tests/integration/s0_core_benchmark.py`는 삭제, 결과 동일).
- **노드 경로:** RV 입구 플레넘(130) → riser(140, 5셀) → 상부 플레넘(150) → Core(유체 3유로 + 고체 + 출력 + 결합) → 하부 플레넘(230)
- **경계 조건:**
  - 입구 J03: 157.1 kg/s, 259 °C
  - 출구 J12: 6.39 MPa
  - p_sys = 출구 압력 고정 (개방 영역이라 하부 플레넘이 소유하지 않음)
- **열 링크:**
  - `barrel_riser`: 배럴 17셀 → riser 5셀. 배럴 반두께 전도 + 대류
  - `riser_rpv`: riser → RPV 5셀. 대류 + 정체 He 간극 3.201–3.275 m
  - `rpv_rccs`: RPV → 30 °C
- **규모:** 347 미지수, 정상해 약 8 s.

| 결과 | 값 | 참고 |
| --- | --- | --- |
| 노심 출구 혼합 | 960.5 K | R1 Table I.8 ≈ 960.15 K |
| 열수지 | 157.1 × Δh = 350 MW − Q_RCCS (상대 1e-7) | |
| RCCS 열손실 | 약 428 kW | 측벽만. PSID 정상 887 kW는 외기 21 °C·실기 기하 |
| RPV 최고 | 463 K | 판정 기준 643.15 K 미만. PSID 정상 내면 495 K |
| 배럴 평균 | 569 K | R3 595 K |
| 노심 Δp | 23.2 kPa | R3 19.95 kPa (하부 블록 D_h 미반영) |
| 연료 최고 (컴팩트 중심) | 1132 K | 균일 출력. R3 1261 K (실제 분포) |
| riser 출구 | 약 535 K | 배럴에서 2.4 MW를 받고 RPV로 0.43 MW를 줌 |

**[시험 가정]**
- 우회 D_h는 4 mm 평행평판으로 두었다.
- 상·하부 헤드는 단열로 두었다.
- riser 외벽 저항은 무시했다.
- riser는 12개 직사각 유로를 면적으로 합쳤고, 대류 면적은 배럴 외면과 riser 외경 기준이다.

## SG-0 (`models/heat_exchanger/`)

Role 5 `HeatExchanger` 상위 조립기에 네 부품을 주입한다:
`HeatExchanger(helium_side=SGHeliumSide(...), wall=NoStorageWall(...), secondary_side=FixedTemperatureWaterSink(...), coupling=SG0Coupling(...), name="sg")`.

- **헬륨 측 `SGHeliumSide`:** SG He 측 9셀(FSH 1 + EES 8) 1D FVM이다. Role 2 `Piping`을 비균일 셀(`z_faces_m`)로 쓴다.
  - 고도는 `heat_exchanger.helium_side.elevation_profile_m`(FSH–EES 간격은 FSH 셀에 포함)이고, 면적은 shroud 환형 9.23 m²(관 단면 미차감)다.
  - 수력직경은 geometry `hydraulic_diameter_m` 0.05766 m [계산]을 쓴다 (2026-11-27 #31, R19 Table 5.3-8). 인자로 덮어쓸 수 있고, 둘 다 없으면 `ModelInputMissing`.
  - 1D 셀 면적은 `mean_flow_area_m2` 6.6622 m²(공극률 0.7218 × 9.23)다. 최소 자유유로면적 3.8458 m²는 SG-2 교차류 상관식용으로 남긴다.
  - 포트는 `helium_inlet`·`helium_outlet`이고, HeliumInventory와 `coupling_state()`를 제공한다.
- **관벽 `NoStorageWall`:** SG-0은 관벽 축열이 없다. 상태·잔차 없이 순열률이 0인지만 확인한다 (에너지 누락 방지).
- **2차측 `FixedTemperatureWaterSink`:** 일정 온도 액체 물이다 (operation.yaml `boundaries.sg_heat_sink`: 466.15 K [제안], 상변화 없음, 물 포트 없음).
  - `calibrate_outlet_K=528.15`(정상 초기화)는 미지수 `{name}.U_overall`과 잔차 T_He,출구 − 528.15 K를 둔다.
  - `U=값`(과도)은 상태가 없다. 둘 중 하나만 받는다.
- **결합 `SG0Coupling`:** Q_i = U A_i (T_He,i − T_물), (to_helium, to_wall_he, to_wall_sec, to_secondary) = (−Q, Q, −Q, Q).
  A_i는 `heat_transfer_area_m2`(관 외면 4020 m², 높이 비례)다. T_He,i는 셀 값(1차 상류에서 셀 출구 상태)이다.
- **단독 시험:**
  - 고정 U에서 He 온도가 이산 해석해 T_i = (ṁc_pT_{i−1} + UA_iT_w)/(ṁc_p + UA_i)와 일치한다.
  - U 보정 시 출구 528.15 K, 제거열 352.2 MW.

## S0 system_loop 첫 계산

steady_state.yaml `conditions.system_loop`의 폐루프를 실제 Component로 조립했다. 아래 표는 손 조립 원형(삭제) 결과다.
Role 5 builder 판은 110·300·330·340·350을 `config/assumptions/provisional_s0.yaml`의 [가정] 체적·유로로 넣어 481 미지수이고,
순환기 승압 36.8 kPa, He 재고 1,762 kg으로 바뀌었다 (나머지는 같은 수준). 최신 값은 `docs/handoffs/role5.md`.
- **루프:** 하부 플레넘(p_sys 소유) → hot duct → SG-0 → [340 대용] → 순환기(유량 157.1 kg/s 지정, 3.1 MW) → [110 대용] → cold annulus
  → RV 입구 → riser(배럴·RPV 링크) → 상부 플레넘 → Core → 하부 플레넘. RPV·RCCS-0도 포함한다.
- **정상 초기화:** 하부 플레넘 질량식을 p_sys − 6.3373 MPa로 바꾼다(D2). SG-0 U를 He 출구 528.15 K에 맞춘다.
- **과도:** 원래 질량식과 보정 U(고정)를 쓴다. 정상해를 이름별로 복사해 일치 초기화한다.
- **규모:** 456 미지수, 정상해 약 20 s.

| 결과 | 값 | PSID 참고 |
| --- | --- | --- |
| SG 제거열 | 352.67 MW (= 노심 350 + 순환기 3.1 − RCCS 0.43, 상대 1e-8) | 352.2 MW |
| SG-0 U (보정) | 약 474 W/(m²·K) | LMTD 추정 421.5 (셀 값 공식이라 더 큼) |
| SG He 출구 | 528.15 K (보정 목표) | 255 °C |
| 노심 입구 / 출구 | 534.4 K / 960.3 K | 532.15 / 960.15 K |
| 순환기 승압 | 33.7 kPa | 91 kPa (K 미보정) |
| SG He 측 Δp | −0.38 kPa (중력 회복 > 마찰) | 번들 20 kPa (K 미보정, D_h 가정) |
| RCCS / RPV 최고 | 427 kW / 463 K | — / 판정 643 K 미만 |
| He 재고량 | 약 1,694 kg | 2,585 kg (미확보 체적 M10) |

- **과도** (출력 +10 %, 3 s, 희소 Jacobian 패턴): 닫힌 루프 질량 보존(1e-8), p_sys 상승, 연료 온도 상승.

**[시험 가정]**
- 수력직경: 우회 4 mm, SG He 측 0.05 m로 두었다.
- 생략·대용: 300 플레넘은 생략하고 hot duct가 SG He 입구로 바로 이어진다. 330·340·110은 대용 배관으로 두었다.
- 헤드와 밸브: 상·하부 헤드는 단열이고 MLSV 밸브는 넣지 않았다.

## 남은 일

1. **K 보정:** PSID 구간 설계손실(순환기 91 kPa, SG 번들 20 kPa 등)에 맞춰 형상손실 K를 S0 system_loop 정상해에서 역산한다 (Role 2와 함께).
   SG He 측 최소 유로면적·수력직경(M11)이 확보되면 교차류 상관식(SG-2)을 쓴다.
2. **헤드 열구조체 HS-09·HS-10과 RCCS 상부 범위:** R1 RCCS 경계는 z 18.68 m까지인데, 모델은 riser 범위(0.85–13.59 m)만 다룬다.
3. **RCCS-1:** 패널·굴뚝 입력(굴뚝 높이 상충)이 확정되면 공기 순환 상태를 둔다.
4. **RPV 온도 차이:** PSID 내면 495 K와 우리 463 K의 차이를 확인해야 한다. 정체 간극을 전도로만 다뤄서일 수 있다 (자연대류 미포함).

## 변수표 patch

`docs/registry_patches/role4/variable_patch.csv` 18행 (update 14, add 4).
- **update 9행:** `T_rpv_wall`(단일 노드 상태), `T_rpv_surface`(별칭), `C_rpv_wall`, `k/rho/c_p_rpv_wall`(값 동일),
  `A_rccs_vessel_rad`, `Q_dot_rccs_surface_rad`, `Q_dot_rccs_vessel_loss`
- **add 2행:** `k_rccs_cavity_air` 0.068, `G_riser_rpv_gap` [가정]
- **SG-0 update 5행:** `h_sg_he_inlet`, `h_sg_he_outlet`, `dp_sg_he`, `Q_dot_sg_removed`, `T_sg_wall`(SG-0 미등록)
- **SG-0 add 2행:** `U_sg0_overall`(보정 미지수), `D_h_sg_he` 0.05 [시험 가정]
- `epsilon_rccs_effective`, `T_rccs_benchmark_sink`, `epsilon_rpv_surface`는 role1 patch가 다룬다.
`equation_patch.csv`는 만들지 않았다 (고체 에너지식은 방정식표의 RCCS-0·용기 식 범위).

## 시험

```bash
python -m pytest -q tests/components/test_vessel_rccs.py     # 5 passed
python -m pytest -q tests/components/test_sg0.py             # 4 passed
python -m pytest -q tests/integration/test_s0_core_benchmark.py   # 4 passed (Role 5 builder 경로, 약 12 s)
python -m pytest -q tests/integration/test_s0_system_loop.py      # 4 passed (약 50 s)
python -m pytest -q                                           # 212 passed (약 107 s)
```
