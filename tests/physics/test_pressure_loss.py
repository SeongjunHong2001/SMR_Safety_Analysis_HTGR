"""Role 2 압력손실 항 (physics/pressure_loss.py) 단독 시험.

팀 보고서 [2] 파트 C §2·§5의 검산값(Churchill f, 거칠기 민감도, 표준 대형 냉각채널 압력강하 분해) 재현,
층류 극한·부호·연속성·입력 처리를 확인한다. 실제 HTGR 수력 정확도 시험은 아니다.
"""
import math

import numpy as np
import pytest

from physics import pressure_loss as pl
from physics.properties import make_helium_properties

HE = make_helium_properties()


@pytest.mark.parametrize("Re,f_report", [(55490, 0.0203), (43800, 0.0214), (36710, 0.0223),
                                         (35330, 0.0225), (30760, 0.0232)])
def test_churchill_matches_team_report_table(Re, f_report):
    # 보고서 [2] C 표 5.2 (매끈한 관, ε = 0)
    assert pl.darcy_friction_factor(Re) == pytest.approx(f_report, abs=6e-5)


def test_roughness_sensitivity_matches_team_report():
    # 보고서 [2] C §2.3: D_h 15.9 mm, Re 5×10⁴에서 매끈관 f = 0.0208, ε 0.01 mm면 +11.0 %
    smooth = pl.darcy_friction_factor(5.0e4)
    assert smooth == pytest.approx(0.0208, abs=1e-4)
    assert pl.darcy_friction_factor(5.0e4, 0.01e-3 / 15.9e-3) / smooth - 1 == pytest.approx(0.110, abs=0.002)


def test_laminar_limit_without_division():
    assert pl.churchill_f_re(0.0) == pytest.approx(64.0, rel=1e-12)
    assert pl.churchill_f_re(0.0, phi_s=pl.PHI_S_PARALLEL_PLATES) == pytest.approx(96.0, rel=1e-12)
    for Re in (1e-6, 1.0, 100.0, 1000.0):
        assert pl.churchill_f_re(Re) == pytest.approx(64.0, rel=1e-9 if Re <= 100 else 1e-3)
    with pytest.raises(ValueError):
        pl.churchill_f_re(-1.0)


def test_friction_is_linear_viscous_at_low_flow_and_odd_in_flow():
    rho, mu, L, D = 3.4, 4.3e-5, 0.793, 0.01588
    A = math.pi * D ** 2 / 4
    for m in (1e-7, 1e-6, 1e-5):                    # Re ≪ 100
        expected = 32.0 * mu * L * m / (rho * D ** 2 * A)
        assert pl.friction_dp(m, rho, mu, L, D, A) == pytest.approx(expected, rel=1e-9)
        assert pl.friction_dp(-m, rho, mu, L, D, A) == pytest.approx(-expected, rel=1e-12)
    assert pl.friction_dp(0.0, rho, mu, L, D, A) == 0.0
    # ṁ = 0을 지나며 도함수가 연속 (층류 선형 기울기)
    h = 1e-9
    slope0 = (pl.friction_dp(h, rho, mu, L, D, A) - pl.friction_dp(-h, rho, mu, L, D, A)) / (2 * h)
    assert slope0 == pytest.approx(32.0 * mu * L / (rho * D ** 2 * A), rel=1e-6)


def test_friction_is_smooth_through_transition():
    rho, mu, L, D = 3.4, 4.3e-5, 1.0, 0.01588
    A = math.pi * D ** 2 / 4
    # Re ≈ 190 … 1.9×10⁴. 천이 구간(Re 2000–4000)은 가파르지만 연속이다: 격자를 반으로 줄이면
    # 기울기의 최대 상대 변화도 반으로 줄어야 한다 (계단이면 줄지 않음).
    def max_slope_change(n):
        m = np.linspace(1e-4, 0.01, n)
        dp = pl.friction_dp(m, rho, mu, L, D, A)
        assert np.all(np.diff(dp) > 0)
        slope = np.diff(dp) / np.diff(m)
        return np.max(np.abs(np.diff(slope)) / slope[:-1])
    coarse, fine = max_slope_change(2001), max_slope_change(4001)
    assert fine < 0.6 * coarse and fine < 0.05


def test_standard_core_channel_pressure_drop_breakdown_matches_team_report():
    # 보고서 [2] C §5.3: 대형 냉각채널 D 15.88 mm, ṁ 20.58 g/s, L 7.93 m, 6.39 MPa, 259 → 687 °C, 40분할
    N, L, D, m, p = 40, 7.93, 0.01588, 0.02058, 6.39e6
    A = math.pi * D ** 2 / 4
    T_face = np.linspace(532.2, 960.2, N + 1)
    T_cell = 0.5 * (T_face[1:] + T_face[:-1])
    rho = HE.density(p, HE.enthalpy(p, T_cell))
    mu = HE.viscosity(p, T_cell)
    dp_fric = np.sum(pl.friction_dp(m, rho, mu, L / N, D, A))
    rho_face = HE.density(p, HE.enthalpy(p, T_face))
    dp_acc = pl.acceleration_dp(m, A, rho_face[0], rho_face[-1])
    dp_grav = np.sum(pl.gravity_dp(rho, -L / N))        # 하향 유동
    assert dp_fric == pytest.approx(14.18e3, rel=0.02)   # 보고서 14.18 kPa (구현 14.04)
    assert dp_acc == pytest.approx(1.50e3, rel=0.01)
    assert dp_grav == pytest.approx(-0.33e3, rel=0.03)


def test_face_friction_uses_half_cells():
    left = dict(rho=4.0, mu=4.0e-5, dz=0.8, D_h=0.0159, A=2e-4)
    right = dict(rho=3.5, mu=4.4e-5, dz=0.4, D_h=0.0159, A=2e-4)
    m = 0.02
    expected = (pl.friction_dp(m, 4.0, 4.0e-5, 0.4, 0.0159, 2e-4) + pl.friction_dp(m, 3.5, 4.4e-5, 0.2, 0.0159, 2e-4))
    assert pl.face_friction_dp(m, left, right) == pytest.approx(expected, rel=1e-12)
    assert pl.face_friction_dp(m, None, right) == pytest.approx(pl.friction_dp(m, 3.5, 4.4e-5, 0.2, 0.0159, 2e-4))


def test_form_loss_direction_and_smoothing():
    rho, A, m_eps = 3.4, 1.0, 1e-3
    q = lambda m: m * abs(m) / (2 * rho * A ** 2)  # noqa: E731
    assert pl.form_loss_dp(10.0, rho, A, 1.0, 0.5, m_eps=m_eps) == pytest.approx(1.0 * q(10.0), rel=1e-6)
    assert pl.form_loss_dp(-10.0, rho, A, 1.0, 0.5, m_eps=m_eps) == pytest.approx(0.5 * q(-10.0), rel=1e-6)
    assert pl.form_loss_dp(0.0, rho, A, 1.0, 0.5, m_eps=m_eps) == 0.0
    m = np.linspace(-5 * m_eps, 5 * m_eps, 1001)
    dp = pl.form_loss_dp(m, rho, A, 1.0, 0.5, m_eps=m_eps)
    assert np.all(np.diff(dp) > 0)                  # 단조, 계단 없음


def test_expansion_contraction_coefficients():
    assert pl.sudden_contraction_K(1.0, 1e-12) == pytest.approx(0.376, abs=1e-3)   # 보고서 [2] C: RELAP5 극한 0.376
    assert pl.sudden_expansion_K(1e-12, 1.0) == pytest.approx(1.0)
    assert pl.sudden_expansion_K(1.0, 1.0) == 0.0 and pl.sudden_contraction_K(1.0, 1.0) == pytest.approx(0.0)
    assert pl.sudden_expansion_K(0.5, 1.0) == pytest.approx(0.25)
    with pytest.raises(ValueError):
        pl.sudden_expansion_K(2.0, 1.0)


def test_gravity_and_acceleration_signs():
    assert pl.gravity_dp(5.0, 2.0) == pytest.approx(5.0 * pl.GRAVITY * 2.0)       # 상향: 하류 압력 감소
    assert pl.face_gravity_dp(5.0, -1.0, 3.0, -1.0) == pytest.approx(-0.5 * pl.GRAVITY * 8.0)
    assert pl.acceleration_dp(10.0, 1.0, 5.0, 3.0) > 0                              # 가열(밀도 감소)이면 양수
    assert pl.acceleration_dp(-10.0, 1.0, 5.0, 3.0) == pl.acceleration_dp(10.0, 1.0, 5.0, 3.0)
    assert pl.junction_acceleration_dp(2.0, 1.0, 0.5, 4.0, 4.0) == pytest.approx(0.5 * (4.0 / (4.0 * 0.25) - 4.0 / 4.0))


def test_input_policy():
    with pytest.raises(ValueError):
        pl.friction_dp(0.01, 3.4, 4e-5, 1.0, 0.0, 1e-4)        # D_h = 0
    with pytest.raises(ValueError):
        pl.friction_dp(0.01, -3.4, 4e-5, 1.0, 0.01, 1e-4)      # ρ < 0
    with pytest.raises(ValueError, match="mode_friction"):
        pl.friction_dp(0.01, 3.4, 4e-5, 1.0, 0.01, 1e-4, mode_friction="blasius")
    assert np.isnan(pl.friction_dp(0.01, np.nan, 4e-5, 1.0, 0.01, 1e-4))   # 비물리 상태 → NaN 전달
    assert pl.mdot_epsilon(157.1) == pytest.approx(0.01571)
    with pytest.raises(ValueError):
        pl.mdot_epsilon(0.0)


def test_face_acceleration_reduces_to_both_limits():
    m, A, rl, rr = 2.0, 1.5e-4, 5.0, 4.0
    assert pl.face_acceleration_dp(m, A, A, rl, rr) == pytest.approx(pl.acceleration_dp(m, A, rl, rr), rel=1e-12)
    assert pl.face_acceleration_dp(m, A, 2 * A, 4.5, 4.5) == pytest.approx(
        pl.junction_acceleration_dp(m, A, 2 * A, 4.5, 4.5), rel=1e-12)
    assert pl.face_acceleration_dp(-m, A, 2 * A, rl, rr) == pl.face_acceleration_dp(m, A, 2 * A, rl, rr)


def test_plenum_boundary_bernoulli_convention():
    # 보고서 [2] D §4.3: 입구 정압강하 (1 + K_in) q, 출구 정압변화 (−1 + K_out) q = 0 (K_in 0.5, K_out 1.0)
    m, A, rho, eps = 0.02058, 1.98e-4, 5.7, 1e-6
    q = m * m / (2 * rho * A * A)
    inlet = pl.plenum_boundary_acceleration_dp(m, A, rho, side="inlet") + pl.form_loss_dp(m, rho, A, 0.5, 1.0, m_eps=eps)
    outlet = pl.plenum_boundary_acceleration_dp(m, A, rho, side="outlet") + pl.form_loss_dp(m, rho, A, 1.0, 0.5, m_eps=eps)
    assert inlet == pytest.approx(1.5 * q, rel=1e-6)          # p_cell = p_plenum − 1.5 q
    assert outlet == pytest.approx(0.0, abs=1e-6 * q)          # p_plenum = p_cell
    rev_in = pl.plenum_boundary_acceleration_dp(-m, A, rho, side="inlet") + pl.form_loss_dp(-m, rho, A, 0.5, 1.0, m_eps=eps)
    assert rev_in == pytest.approx(0.0, abs=1e-6 * q)          # 역류로 입구 Face가 출구가 됨: 순변화 0
    with pytest.raises(ValueError):
        pl.plenum_boundary_acceleration_dp(m, A, rho, side="middle")
