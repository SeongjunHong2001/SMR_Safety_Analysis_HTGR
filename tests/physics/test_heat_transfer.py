"""Role 3 공통 열전달 커널과 헬륨 대류 상관식 단독 시험.

physics/heat_transfer/{conduction,convection,radiation,correlations}.py의 식·부호·shape·정의역,
팀 보고서 [2] 파트 C의 검산값 재현, 천이·냉각·역류에서의 연속성을 확인한다.
실제 HTGR 열전달 정확도 시험은 아니다.
"""
import math
import warnings

import numpy as np
import pytest

from interfaces.contracts import HeatTransferCorrelation
from physics import heat_transfer as ht
from physics.heat_transfer.correlations import (NU_LAM_CIRCULAR, NU_LAM_PARALLEL_PLATES,
                                                McEligotConvection, make_convection_correlation)
from physics.heat_transfer.radiation import STEFAN_BOLTZMANN, effective_emissivity_gray_enclosure
from physics.properties import make_helium_properties

HE = make_helium_properties()


def _he_state(p, T_b, T_w, m_dot, *, rho=True):
    state = dict(m_dot=m_dot, T_bulk=T_b, T_wall=T_w, mu=HE.viscosity(p, T_b),
                 k=HE.conductivity(p, T_b), cp=HE.cp(p, T_b))
    if rho:
        state["rho"] = HE.density(p, HE.enthalpy(p, T_b))
    return state


# ------------------------------------------------------------------ 기본 커널
def test_kernel_formulas_and_sign_convention():
    assert ht.conduction(400.0, 300.0, 2.0, 3.0, 0.5) == pytest.approx(2.0 * 3.0 * 100.0 / 0.5)
    assert ht.conduction(300.0, 400.0, 2.0, 3.0, 0.5) == pytest.approx(-1200.0)
    assert ht.convection(400.0, 300.0, 10.0, 2.0) == pytest.approx(2000.0)
    assert ht.radiation(1000.0, 300.0, 0.8, 1.0, STEFAN_BOLTZMANN) == pytest.approx(
        0.8 * STEFAN_BOLTZMANN * (1000.0 ** 4 - 300.0 ** 4))
    assert ht.heat_transfer_coefficient(90.0, 0.3, 0.015) == pytest.approx(1800.0)
    assert ht.convection(350.0, 350.0, 10.0, 2.0) == 0.0


def test_public_functions_stay_callable_after_submodule_import():
    import physics.heat_transfer.conduction  # noqa: F401  (예전에는 이 import가 공개 함수를 모듈로 덮었다)
    import physics.heat_transfer.radiation  # noqa: F401
    for _ in range(3):
        assert callable(ht.conduction) and callable(ht.convection) and callable(ht.radiation)
        assert np.isfinite(ht.conduction(400, 300, 1, 1, 1))


def test_stefan_boltzmann_matches_registry_value():
    assert math.isclose(STEFAN_BOLTZMANN, 5.6704e-8, rel_tol=1e-5)   # 변수표 sigma_SB


def test_kernels_broadcast_and_return_scalars():
    T_a = np.array([[400.0], [500.0]])
    T_b = np.array([300.0, 350.0, 380.0])
    assert ht.conduction(T_a, T_b, 1.0, 1.0, 1.0).shape == (2, 3)
    assert ht.convection(T_a, T_b, np.array([1.0, 2.0, 3.0]), 1.0).shape == (2, 3)
    assert ht.radiation(T_a, T_b, 0.5, np.array([1.0, 2.0, 3.0]), STEFAN_BOLTZMANN).shape == (2, 3)
    assert np.ndim(ht.conduction(400, 300, 1, 1, 1)) == 0


def test_invalid_geometry_and_coefficients_raise():
    with pytest.raises(ValueError):
        ht.conduction(400, 300, 1, 1, 0.0)          # L = 0
    with pytest.raises(ValueError):
        ht.conduction(400, 300, -1.0, 1, 1)         # k < 0
    with pytest.raises(ValueError):
        ht.convection(400, 300, 1, -2.0)            # A < 0
    with pytest.raises(ValueError):
        ht.radiation(400, 300, 1.2, 1, STEFAN_BOLTZMANN)   # epsilon > 1
    with pytest.raises(ValueError):
        ht.heat_transfer_coefficient(10, 0.3, 0.0)  # D_h = 0
    assert ht.convection(400, 300, 5.0, 0.0) == 0.0  # 면적 0은 허용 (사용하지 않는 셀)


def test_nonphysical_temperatures_and_nan_coefficients_propagate_nan():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        q = ht.convection(np.array([400.0, -1.0, np.nan, 400.0]), 300.0, np.array([1.0, 1.0, 1.0, np.nan]), 1.0)
        r = ht.radiation(np.array([0.0, 500.0]), 300.0, 0.5, 1.0, STEFAN_BOLTZMANN)
    assert q[0] == 100.0 and np.all(np.isnan(q[1:]))
    assert np.isnan(r[0]) and r[1] > 0


def test_effective_emissivity_rpv_to_rccs_matches_registry():
    # 변수표 epsilon_rccs_effective 0.74 (RPV ε 0.85, RCCS ε 0.80, 동심원통 r 3.408 / 4.633 m)
    eps = effective_emissivity_gray_enclosure(0.85, 0.80, 3.408 / 4.633)
    assert eps == pytest.approx(0.7351, abs=1e-4)
    assert effective_emissivity_gray_enclosure(1.0, 1.0, 0.5) == pytest.approx(1.0)
    assert effective_emissivity_gray_enclosure(0.8, 0.8, 1.0) == pytest.approx(1 / (2 / 0.8 - 1))
    with pytest.raises(ValueError):
        effective_emissivity_gray_enclosure(0.8, 0.8, 1.5)


# ------------------------------------------------------------------ 대류 상관식
def test_correlation_satisfies_protocol_and_factory():
    assert isinstance(McEligotConvection(), HeatTransferCorrelation)
    assert make_convection_correlation("mceligot", include_nc=False).include_nc is False
    with pytest.raises(ValueError, match="unsupported or unimplemented mode_nu"):
        make_convection_correlation("dittus_boelter")


@pytest.mark.parametrize("r_T,expected", [(1.15, 95.2), (1.30, 89.6)])
def test_rated_mceligot_matches_team_report(r_T, expected):
    # 보고서 [2] C 표: Re = 5×10⁴, T_w/T_b = 1.15–1.30에서 Nu = 95.2–89.6 (Pr ≈ 0.66)
    mu, k, D, A = 4.0e-5, 0.3, 0.0159, 1.0
    state = dict(m_dot=5.0e4 * A * mu / D, T_bulk=800.0, T_wall=800.0 * r_T, mu=mu, k=k, cp=0.66 * k / mu)
    out = McEligotConvection(include_nc=False).evaluate(state, dict(D_h=D, A_flow=A))
    assert out["Nu"] == pytest.approx(expected, abs=0.1)
    assert out["mceligot_nominal"]


def test_transition_weights_match_team_report():
    mu, k, D, A = 4.0e-5, 0.3, 0.0159, 1.0
    corr = McEligotConvection(include_nc=False)
    w = []
    for Re in (2300.0, 3150.0, 4000.0):
        state = dict(m_dot=Re * A * mu / D, T_bulk=800.0, T_wall=800.0, mu=mu, k=k, cp=0.66 * k / mu)
        w.append(float(corr.evaluate(state, dict(D_h=D, A_flow=A))["w_tr"]))
    # 보고서 값은 tanh 안에 Re_b를 쓴 결과다. 구현은 역류 대칭·매끄러움을 위해 Re_safe = √(Re² + 100²)를 쓰므로
    # Re = 3150에서 0.502로 0.002 차이가 난다 (의도한 차이).
    assert w == pytest.approx([0.0077, 0.500, 0.9923], abs=3e-3)
    assert w[1] == pytest.approx(0.5 * (1 + math.tanh((math.hypot(3150.0, 100.0) - 3150.0) / 350.0)), rel=1e-9)


def test_rated_core_channel_with_helium_properties():
    # R1 연료 냉각채널 (d = 15.88 mm) 정격 근처: 노심 입구 6.37 MPa, 벌크 750 K, 벽 850 K
    p, D = 6.3718e6, 2 * 0.00794
    A = math.pi * D ** 2 / 4
    m_hole = 137.7 / (65 * 102)                  # [가정] 평균채널 유량을 기둥 65 × 큰 냉각공 102로 나눔
    out = McEligotConvection().evaluate(_he_state(p, 750.0, 850.0, m_hole), dict(D_h=D, A_flow=A))
    assert 3.0e4 < out["Re"] < 1.0e5
    assert 0.65 < out["Pr"] < 0.68
    assert out["Ri"] < 0.1                         # 강제대류 지배
    assert out["Nu"] == pytest.approx(out["Nu_turb"], rel=1e-6)
    assert 1500 < out["h_conv"] < 3000


def test_zero_and_reverse_flow_are_finite_and_symmetric():
    p, D = 6.39e6, 0.0159
    geo = dict(D_h=D, A_flow=math.pi * D ** 2 / 4)
    corr = McEligotConvection()
    h0 = corr.coefficient(_he_state(p, 900.0, 950.0, 0.0), geo)
    assert np.isfinite(h0) and h0 > 0
    for m in (1e-3, 5e-3, 0.05):
        fwd = corr.coefficient(_he_state(p, 900.0, 950.0, m), geo)
        rev = corr.coefficient(_he_state(p, 900.0, 950.0, -m), geo)
        assert fwd == pytest.approx(rev, rel=1e-12)


def test_laminar_floor_and_shape():
    p = 6.39e6
    st = _he_state(p, 900.0, 900.0, 0.0, rho=False)
    circ = McEligotConvection(include_nc=False).evaluate(st, dict(D_h=0.0159, A_flow=2e-4))
    plates = McEligotConvection(include_nc=False).evaluate(st, dict(D_h=0.004, A_flow=1e-3, shape="parallel_plates"))
    assert circ["Nu"] == pytest.approx(NU_LAM_CIRCULAR, rel=1e-6)
    assert plates["Nu"] == pytest.approx(NU_LAM_PARALLEL_PLATES, rel=1e-6)
    with pytest.raises(ValueError, match="shape"):
        McEligotConvection().evaluate(st | {"rho": 3.4}, dict(D_h=0.0159, A_flow=2e-4, shape="square"))


def test_plofc_natural_convection_matches_team_report():
    # 보고서 [2] C §3.2: 6.39 MPa, T_b = 900 K, ΔT = 20–150 K → Ra_Dh 3.6×10³–2.7×10⁴, Nu 4.36 → 약 6.7
    p, D = 6.39e6, 0.0159
    geo = dict(D_h=D, A_flow=math.pi * D ** 2 / 4)
    low = McEligotConvection().evaluate(_he_state(p, 900.0, 920.0, 0.0), geo)
    high = McEligotConvection().evaluate(_he_state(p, 900.0, 1050.0, 0.0), geo)
    off = McEligotConvection(include_nc=False).evaluate(_he_state(p, 900.0, 1050.0, 0.0), geo)
    assert low["Ra"] == pytest.approx(3.6e3, rel=0.05)
    assert high["Ra"] == pytest.approx(2.7e4, rel=0.05)
    assert 6.6 < high["Nu_nc"] < 6.9 and high["Nu"] > high["Nu_nc"]
    assert off["Nu"] == pytest.approx(NU_LAM_CIRCULAR, rel=1e-6)     # include_nc=False 보수 민감도
    assert high["Ri"] > 1.0          # 정체 유량에서 Ri = Gr / Re_safe²는 Re_min_safe(100)로 상한이 정해짐 (약 4)


def test_coefficient_is_continuous_across_transition_and_cooling():
    p, D = 6.39e6, 0.0159
    A = math.pi * D ** 2 / 4
    corr = McEligotConvection()
    mu = HE.viscosity(p, 800.0)
    Re = np.linspace(1500.0, 6000.0, 901)
    h = np.array([corr.coefficient(_he_state(p, 800.0, 850.0, r * A * mu / D), dict(D_h=D, A_flow=A)) for r in Re])
    jumps = np.abs(np.diff(h)) / h[:-1]
    assert jumps.max() < 0.01                                          # 계단 없이 매끄러운 천이
    Tw = np.linspace(780.0, 820.0, 401)                                  # 냉각 ↔ 가열 경계 (T_w = T_b = 800 K)
    m_rated = 5.0e4 * A * mu / D
    h2 = np.array([corr.coefficient(_he_state(p, 800.0, t, m_rated), dict(D_h=D, A_flow=A)) for t in Tw])
    assert (np.abs(np.diff(h2)) / h2[:-1]).max() < 1e-3


def test_correlation_input_errors():
    corr = McEligotConvection()
    st = _he_state(6.39e6, 900.0, 950.0, 0.01)
    with pytest.raises(KeyError, match="D_h"):
        corr.coefficient(st, dict(A_flow=1e-4))
    with pytest.raises(KeyError, match="rho"):
        corr.coefficient({k: v for k, v in st.items() if k != "rho"}, dict(D_h=0.0159, A_flow=2e-4))
    with pytest.raises(ValueError):
        corr.coefficient(st | {"mu": -1.0}, dict(D_h=0.0159, A_flow=2e-4))
    with pytest.raises(ValueError):
        McEligotConvection(Re_tr_width=0.0)
    with pytest.raises(TypeError):
        McEligotConvection(include_nc=1)
    nan_out = corr.coefficient(st | {"T_wall": -5.0}, dict(D_h=0.0159, A_flow=2e-4))
    assert np.isnan(nan_out)
