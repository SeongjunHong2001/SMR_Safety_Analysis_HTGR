"""Role 2 1D FVM 보존 커널 (physics/flow.py) 단독 시험.

팀 보고서 [2] 파트 B §5.3의 단위 검증 계획을 따른다:
① 단일 1D 가열 채널 정상상태 질량·에너지 보존, ② 밀폐 정적 가열의 −V ṗ_sys 압력 상승률 (γ − 1) Q,
③ 유동 역전(ṁ = 0 통과) 연속성과 헬륨 물성 연계. 실제 HTGR 열수력 정확도 시험은 아니다.
"""
import numpy as np
import pytest

from physics import flow
from physics.properties import make_helium_properties

HE = make_helium_properties()
P_SYS = 6.3373e6


def _channel_state(T_cells):
    h = HE.enthalpy(P_SYS, T_cells)
    return dict(h=h, rho=HE.density(P_SYS, h), drho_dp=HE.drho_dp(P_SYS, h), drho_dh=HE.drho_dh(P_SYS, h))


def test_smooth_upwind_limits():
    m_eps = 1e-3
    assert flow.smooth_upwind(1.0, 10.0, 20.0, m_eps) == pytest.approx(10.0)
    assert flow.smooth_upwind(-1.0, 10.0, 20.0, m_eps) == pytest.approx(20.0)
    assert flow.smooth_upwind(0.0, 10.0, 20.0, m_eps) == pytest.approx(15.0)
    m = np.linspace(-5e-3, 5e-3, 101)
    assert np.all(np.diff(flow.smooth_upwind(m, 10.0, 20.0, m_eps)) <= 0)


def test_heated_channel_steady_state_conserves_mass_and_energy():
    # ① N셀 가열 채널, 정상 (ẏ = 0): 셀 엔탈피가 열수지대로면 질량·에너지 잔차가 0
    N, m, Q_total = 17, 137.7, 3.5e8 * 0.98
    q = np.full(N, Q_total / N)
    h_in = HE.enthalpy(P_SYS, 532.15)
    h_cells = h_in + np.cumsum(q) / m                  # 1차 상류: 셀 i의 상태 = 출구 상태
    st = _channel_state(HE.temperature(P_SYS, h_cells))
    m_faces = np.full(N + 1, m)
    out = flow.cell_residuals(V=np.full(N, 15.4 / N), rho=st["rho"], drho_dp=st["drho_dp"], drho_dh=st["drho_dh"],
                              h=h_cells, h_dot=0.0, p_sys_dot=0.0, q_wall=q, m_faces=m_faces,
                              h_inlet_boundary=h_in, h_outlet_boundary=h_cells[-1], m_eps=flow_eps(m))
    assert out["R_mass"].shape == (N,) and out["h_face"].shape == (N + 1,)
    assert np.max(np.abs(out["R_mass"])) < 1e-10 * m
    assert np.max(np.abs(out["R_energy"])) < 1e-10 * Q_total
    # 출구 온도: 532.15 K + Q/(ṁ c_p)
    assert HE.temperature(P_SYS, h_cells[-1]) == pytest.approx(532.15 + Q_total / (m * HE.cp_const), rel=1e-12)
    assert np.sum(flow.energy_flux_balance(m_faces, out["h_face"])) == pytest.approx(Q_total, rel=1e-12)


def flow_eps(m_nom):
    from physics.pressure_loss import mdot_epsilon
    return mdot_epsilon(m_nom)


def test_closed_rigid_volume_pressure_rise_is_gamma_minus_one():
    # ② 밀폐·정적 가열 (면 유량 0): V ṗ_sys = (γ − 1) Q = R/c_v Q (보고서 B §2.4, §5.3)
    V, Q, T = 20.0, 1.0e6, 900.0
    st = _channel_state(np.array([T]))
    rho, rp, rh, h = st["rho"][0], st["drho_dp"][0], st["drho_dh"][0], st["h"][0]
    # 미지수 (ḣ, ṗ): 에너지 ρVḣ − Vṗ = Q, 질량 V(ρ_p ṗ + ρ_h ḣ) = 0 → 잔차 함수로 선형계 구성
    def res(x):
        h_dot, p_dot = x
        r_m = flow.mass_residual(V, rp, rh, p_dot, h_dot, 0.0, 0.0)
        r_e = flow.energy_residual(rho, V, h, h_dot, p_dot, 0.0, h, 0.0, h, Q)
        return np.array([r_m, r_e])
    J = np.column_stack([res([1.0, 0.0]) - res([0.0, 0.0]), res([0.0, 1.0]) - res([0.0, 0.0])])
    h_dot, p_dot = np.linalg.solve(J, -res([0.0, 0.0]))
    gamma_minus_1 = HE.R / (HE.cp_const - HE.R)
    assert V * p_dot == pytest.approx(gamma_minus_1 * Q, rel=1e-10)
    assert gamma_minus_1 == pytest.approx(0.6664, abs=1e-4)
    assert np.allclose(res([h_dot, p_dot]), 0.0, atol=1e-6)


def test_enthalpy_form_equals_conservative_form():
    # R_energy + h R_mass = d(ρVh)/dt − Vṗ + Σ면 엔탈피 플럭스 − Q (T_ref 무관, 보고서 B §3.2)
    rng = np.random.default_rng(0)
    N = 6
    T = np.linspace(600.0, 900.0, N)
    st = _channel_state(T)
    V = np.full(N, 0.3)
    h_dot = rng.normal(0, 50.0, N)
    p_dot = 1.2e3
    m_faces = np.array([2.0, 1.9, 1.7, 1.8, 2.1, 2.0, 1.95])
    q = rng.normal(1e5, 1e4, N)
    out = flow.cell_residuals(V=V, rho=st["rho"], drho_dp=st["drho_dp"], drho_dh=st["drho_dh"], h=st["h"],
                              h_dot=h_dot, p_sys_dot=p_dot, q_wall=q, m_faces=m_faces,
                              h_inlet_boundary=HE.enthalpy(P_SYS, 550.0), h_outlet_boundary=st["h"][-1],
                              m_eps=1e-4)
    rho_dot = st["drho_dp"] * p_dot + st["drho_dh"] * h_dot
    d_rhoVh = V * (rho_dot * st["h"] + st["rho"] * h_dot)
    conservative = d_rhoVh - V * p_dot + flow.energy_flux_balance(m_faces, out["h_face"]) - q
    assert np.allclose(out["R_energy"] + st["h"] * out["R_mass"], conservative, rtol=1e-12, atol=1e-6)


def test_reverse_flow_takes_downstream_boundary_enthalpy():
    N = 3
    h = np.array([3.0e6, 4.0e6, 5.0e6])
    h_in, h_out = 2.5e6, 2.7e6
    fwd = flow.face_enthalpies(np.full(N + 1, 1.0), h, h_in, h_out, 1e-6)
    rev = flow.face_enthalpies(np.full(N + 1, -1.0), h, h_in, h_out, 1e-6)
    assert np.allclose(fwd, [h_in, *h])
    assert np.allclose(rev, [*h, h_out])
    m_in_port, m_out_port = flow.port_inflows_from_faces(np.full(N + 1, -1.0))
    assert m_in_port == -1.0 and m_out_port == 1.0        # 역류: 출구 포트로 유입


def test_flow_reversal_is_continuous():
    # ③ ṁ이 0을 지나도 에너지 잔차가 연속 (C∞ 상류차분)
    N = 4
    st = _channel_state(np.linspace(700.0, 900.0, N))
    m_eps = 1e-4
    vals = []
    for m in np.linspace(-1e-3, 1e-3, 401):
        out = flow.cell_residuals(V=0.5, rho=st["rho"], drho_dp=st["drho_dp"], drho_dh=st["drho_dh"], h=st["h"],
                                  h_dot=0.0, p_sys_dot=0.0, q_wall=0.0, m_faces=np.full(N + 1, m),
                                  h_inlet_boundary=HE.enthalpy(P_SYS, 650.0),
                                  h_outlet_boundary=HE.enthalpy(P_SYS, 950.0), m_eps=m_eps)
        vals.append(out["R_energy"])
    vals = np.array(vals)
    step = np.max(np.abs(np.diff(vals, axis=0)), axis=0)
    span = np.max(np.abs(vals), axis=0)
    assert np.all(step < 0.05 * span)


def test_momentum_residual_and_inertia_switch():
    r = flow.momentum_residual(6.37e6, 6.35e6, dp_acc=1.5e3, dp_fric=14e3, dp_form=1.4e3, dp_grav=-0.33e3)
    assert r == pytest.approx(-2.0e4 + 1.5e3 + 14e3 + 1.4e3 - 0.33e3)
    assert flow.momentum_residual(6.29e6, 6.38e6, dp_comp=9.0e4) == pytest.approx(0.0)   # 순환기 승압
    lA = flow.face_inertia_length(0.8, 1.35, 0.8, 1.35)
    assert lA == pytest.approx(0.8 / 1.35)
    base = flow.momentum_residual(1.0, 1.0)
    assert flow.momentum_residual(1.0, 1.0, include_inertia=False, l_over_A=lA, m_dot_dot=5.0) == base
    assert flow.momentum_residual(1.0, 1.0, include_inertia=True, l_over_A=lA, m_dot_dot=5.0) == pytest.approx(5.0 * lA)
    assert flow.face_flow_state_kind(False) == "algebraic"
    assert flow.face_flow_state_kind(True) == "differential"
    with pytest.raises(TypeError):
        flow.face_flow_state_kind(0)


def test_shape_checks():
    with pytest.raises(ValueError, match="N\\+1"):
        flow.face_enthalpies(np.zeros(3), np.zeros(3), 0.0, 0.0, 1e-4)
    with pytest.raises(ValueError):
        flow.cell_residuals(V=1.0, rho=1.0, drho_dp=1.0, drho_dh=1.0, h=np.zeros((2, 2)), h_dot=0.0, p_sys_dot=0.0,
                            q_wall=0.0, m_faces=np.zeros(3), h_inlet_boundary=0.0, h_outlet_boundary=0.0, m_eps=1e-4)
