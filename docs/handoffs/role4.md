# Role 4 SG·압력용기·RCCS 인수인계

브랜치는 `role3/core-solid-power-coupling`이고, Role 2 유동·장치 커밋 `d7130ec` 위에 Role 3 노심 고체와 함께 미커밋으로 있다.
커밋할 때 Role 3 파일과 Role 4 파일을 나눈다.
병합 순서는 Role 1 → Role 2 물성 → Role 3 열전달 → Role 2 유동·장치 → Role 3 노심 고체 → Role 4다.

## 진행 상태

| 항목 | 파일 | 상태 |
| --- | --- | --- |
| RPV 측벽 | `models/vessel/vessel.py` `Vessel` | 구현 (아래) |
| 표면–유체 대류 링크 | `models/vessel/vessel.py` `SurfaceFluidConvection`, `helium_gap_conductance`, `axial_overlap_matrix` | 구현 |
| RCCS-0 고정 온도 열침원 | `models/rccs/rccs.py` `RCCS0Sink` | 구현 |
| RCCS-1 공기 순환 | `models/rccs/rccs.py` `RCCS` | 미구현 (패널·굴뚝 입력 상충, chimney_height_m null) |
| SG (HeatExchanger 자식 helium_side·wall·secondary_side, coupling) | `models/heat_exchanger/` | 미구현 — 다음 단계 (system_loop에 필요) |
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

## S0 core_benchmark 조립 원형과 첫 계산 (`tests/integration/s0_core_benchmark.py`)

steady_state.yaml `conditions.core_benchmark`의 영역과 경계를 실제 Component로 손 조립했다. Role 5 builder의 원형이다.
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

## 남은 일

1. **SG (system_loop용):** `helium_side`·`wall`·`secondary_side`·`coupling` 구현. SG-0은 일정 온도 물(193 °C 제안)과 U 보정을 쓴다.
   He 측 최소 유로면적(M11)이 없어 SG-2 교차류 상관식은 아직 쓸 수 없다.
2. **헤드 열구조체 HS-09·HS-10과 RCCS 상부 범위:** R1 RCCS 경계는 z 18.68 m까지인데, 모델은 riser 범위(0.85–13.59 m)만 다룬다.
3. **RCCS-1:** 패널·굴뚝 입력(굴뚝 높이 상충)이 확정되면 공기 순환 상태를 둔다.
4. **RPV 온도 차이:** PSID 내면 495 K와 우리 463 K의 차이를 확인해야 한다. 정체 간극을 전도로만 다뤄서일 수 있다 (자연대류 미포함).

## 변수표 patch

`docs/registry_patches/role4/variable_patch.csv` 11행 (update 9, add 2).
- **update 9행:** `T_rpv_wall`(단일 노드 상태), `T_rpv_surface`(별칭), `C_rpv_wall`, `k/rho/c_p_rpv_wall`(값 동일),
  `A_rccs_vessel_rad`, `Q_dot_rccs_surface_rad`, `Q_dot_rccs_vessel_loss`
- **add 2행:** `k_rccs_cavity_air` 0.068, `G_riser_rpv_gap` [가정]
- `epsilon_rccs_effective`, `T_rccs_benchmark_sink`, `epsilon_rpv_surface`는 role1 patch가 다룬다.
`equation_patch.csv`는 만들지 않았다 (고체 에너지식은 방정식표의 RCCS-0·용기 식 범위).

## 시험

```bash
python -m pytest -q tests/components/test_vessel_rccs.py     # 5 passed
python -m pytest -q tests/integration/test_s0_core_benchmark.py   # 3 passed (약 8 s)
python -m pytest -q                                           # 194 passed (약 60 s)
```
