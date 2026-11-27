"""Role 2 헬륨 물성 provider (physics/properties.py) 단독 시험.

R1 Table AIV.13 식, Role 1 물성표(data/properties/helium.csv)와 정격 분포(config/operation.yaml)와의
일치, 왕복·도함수·broadcasting·정의역 처리를 확인한다. 실제 HTGR 열수력 정확도 시험은 아니다.
"""
import csv
import math
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

from interfaces.contracts import HeliumProperties
from physics.properties import (
    DEFAULT_R1_PROPERTY_FILE,
    R1HeliumProperties,
    make_helium_properties,
)

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def he():
    return R1HeliumProperties.from_yaml()


@pytest.fixture(scope="module")
def r1_yaml():
    with DEFAULT_R1_PROPERTY_FILE.open(encoding="utf-8") as f:
        return yaml.safe_load(f)["helium"]


def test_satisfies_helium_properties_protocol(he):
    assert isinstance(he, HeliumProperties)


def test_constants_come_from_r1_file_and_are_consistent(he, r1_yaml):
    assert he.R == r1_yaml["specific_gas_constant_J_kgK"]
    assert math.isclose(he.R, r1_yaml["universal_gas_constant_J_molK"] / r1_yaml["molar_mass_kg_mol"], rel_tol=1e-6)
    assert he.cp_const == r1_yaml["cp_J_kgK"]
    # 팀 결정(2026-11-26): R1 세트 R 2077.46, c_p 5195, c_v = c_p − R = 3117.54
    assert math.isclose(he.R, 2077.46, rel_tol=1e-9)
    assert math.isclose(float(he.cv(6.39e6, 900.0)), 3117.54, abs_tol=1e-9)
    assert math.isclose(float(he.cp(6.39e6, 900.0) - he.cv(6.39e6, 900.0)), he.R, rel_tol=1e-12)


def test_matches_role1_helium_table(he):
    with (ROOT / "data/properties/helium.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 50
    for r in rows:
        T, p = float(r["temperature_K"]), float(r["pressure_Pa"])
        h = he.enthalpy(p, T)
        assert math.isclose(he.density(p, h), float(r["density_kg_m3"]), rel_tol=2e-5), r
        assert math.isclose(he.cp(p, T), float(r["cp_J_kgK"]), rel_tol=1e-12), r
        assert math.isclose(he.viscosity(p, T), float(r["viscosity_Pa_s"]), rel_tol=2e-5), r
        assert math.isclose(he.conductivity(p, T), float(r["conductivity_W_mK"]), rel_tol=2e-5), r


def test_r1_conductivity_formula_in_bar_units(he):
    # R1 원식: k = 2.682e-3 (1 + 1.123e-3 P) T^(0.71 (1 − 2e-4 P)),  P [bar]
    for T, p_bar in [(532.15, 63.9), (960.15, 63.373), (1500.0, 1.0)]:
        ref = 2.682e-3 * (1 + 1.123e-3 * p_bar) * T ** (0.71 * (1 - 2e-4 * p_bar))
        assert math.isclose(he.conductivity(p_bar * 1e5, T), ref, rel_tol=1e-12)


def test_ideal_gas_and_enthalpy_round_trip(he):
    T = np.linspace(300.0, 1900.0, 41)
    p = np.full_like(T, 6.3373e6)
    h = he.enthalpy(p, T)
    assert np.allclose(he.temperature(p, h), T, rtol=0, atol=1e-10)
    assert np.allclose(p, he.density(p, h) * he.R * T, rtol=1e-13)
    # 기본 T_ref = 0 K에서는 h = c_p T. T = 0 K 자체는 비물리 입력이라 NaN이다.
    assert math.isclose(he.enthalpy(6.0e6, 1000.0), he.cp_const * 1000.0, rel_tol=1e-15)
    assert np.isnan(he.enthalpy(6.0e6, he.T_ref))
    assert R1HeliumProperties.from_yaml(T_ref=298.15).enthalpy(6.0e6, 298.15) == 0.0


def test_enthalpy_reference_shifts_h_but_not_state(he):
    shifted = R1HeliumProperties.from_yaml(T_ref=298.15)
    p, T = 6.3373e6, 960.15
    assert math.isclose(shifted.enthalpy(p, T), he.enthalpy(p, T) - he.cp_const * 298.15, rel_tol=1e-12)
    assert math.isclose(shifted.density(p, shifted.enthalpy(p, T)), he.density(p, he.enthalpy(p, T)), rel_tol=1e-13)
    assert math.isclose(shifted.drho_dh(p, shifted.enthalpy(p, T)), he.drho_dh(p, he.enthalpy(p, T)), rel_tol=1e-13)


@pytest.mark.parametrize("p,T", [(1.0e5, 300.0), (6.3373e6, 532.15), (6.3373e6, 960.15), (7.18e6, 1900.0)])
def test_density_derivatives_match_finite_differences(he, p, T):
    h = he.enthalpy(p, T)
    dp, dh = p * 1e-6, abs(h) * 1e-6
    fd_p = (he.density(p + dp, h) - he.density(p - dp, h)) / (2 * dp)
    fd_h = (he.density(p, h + dh) - he.density(p, h - dh)) / (2 * dh)
    assert math.isclose(he.drho_dp(p, h), fd_p, rel_tol=1e-7)
    assert math.isclose(he.drho_dh(p, h), fd_h, rel_tol=1e-7)
    assert he.drho_dp(p, h) > 0 > he.drho_dh(p, h)


def test_broadcasting_and_scalar_outputs(he):
    T = np.array([500.0, 700.0, 900.0])
    p = np.array([[5.0e6], [6.0e6]])
    for fn in (he.viscosity, he.conductivity, he.cp, he.cv, he.enthalpy):
        assert np.shape(fn(p, T)) == (2, 3), fn.__name__
    h = he.enthalpy(6.0e6, T)
    for fn in (he.temperature, he.density, he.drho_dp, he.drho_dh):
        assert np.shape(fn(p, h)) == (2, 3), fn.__name__
        assert np.shape(fn(6.0e6, h)) == (3,), fn.__name__
    scalar = he.density(6.0e6, float(h[0]))
    assert np.ndim(scalar) == 0 and isinstance(float(scalar), float)


def test_nonphysical_inputs_give_nan_without_warnings(he):
    p = np.array([6.0e6, 0.0, -1.0, 6.0e6, np.nan, 6.0e6])
    T = np.array([600.0, 600.0, 600.0, 0.0, 600.0, -50.0])
    h = he.enthalpy(6.0e6, T)          # T ≤ 0 → h도 NaN
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        values = {fn.__name__: fn(p, T) for fn in (he.viscosity, he.conductivity, he.cp, he.cv)}
        values["density"] = he.density(p, he.enthalpy(6.0e6, np.full_like(T, 600.0)))
        values["temperature"] = he.temperature(6.0e6, np.array([-1.0e6, 1.0e6]))
    for name, v in values.items():
        if name == "temperature":
            assert np.isnan(v[0]) and np.isfinite(v[1])
            continue
        assert np.isfinite(v[0]), name
        assert np.all(np.isnan(v[1:5]) if name != "density" else np.isnan(v[1:3])), name
    assert np.isnan(h[3]) and np.isnan(h[5])


def test_valid_range_check(he):
    assert he.in_valid_range(6.3373e6, 960.15)
    assert not he.in_valid_range(6.3373e6, 2000.0)
    assert not he.in_valid_range(1.0e4, 600.0)
    assert np.array_equal(he.in_valid_range(6.0e6, np.array([299.0, 300.0, 1900.0, 1901.0])),
                          [False, True, True, False])
    he.check_valid_range(np.array([1.0e5, 7.18e6]), np.array([300.0, 1900.0]))
    with pytest.raises(ValueError, match="outside R1 helium range"):
        he.check_valid_range(6.0e6, np.array([500.0, 2500.0]), where="test cells")


def test_rated_distribution_states_are_in_range_and_consistent(he):
    with (ROOT / "config/operation.yaml").open(encoding="utf-8") as f:
        stations = yaml.safe_load(f)["initial_distribution"]["stations"]
    p = np.array([s["p_Pa"] for s in stations], dtype=float)
    T = np.array([s["T_K"] for s in stations], dtype=float)
    rho = np.array([s["rho_kg_m3"] for s in stations], dtype=float)
    he.check_valid_range(p, T, where="operation.yaml initial_distribution")
    assert np.allclose(he.density(p, he.enthalpy(p, T)), rho, rtol=2e-3)


def test_factory_selects_mode_explicitly():
    props = make_helium_properties("ideal_gas")
    assert isinstance(props, R1HeliumProperties)
    assert make_helium_properties("ideal_gas", T_ref=298.15).T_ref == 298.15
    with pytest.raises(ValueError, match="unsupported mode_he_eos"):
        make_helium_properties("kta_real_gas")


def test_constructor_rejects_bad_coefficients(tmp_path):
    good = dict(R=2077.46, cp_const=5195.0, mu_a=3.674e-7, mu_b=0.7, k_a=2.682e-3, k_b=1.123e-8, k_c=0.71, k_d=2e-9)
    R1HeliumProperties(**good)
    for key, bad in [("R", 0.0), ("cp_const", 1000.0), ("mu_a", -1.0), ("k_c", float("nan"))]:
        with pytest.raises(ValueError):
            R1HeliumProperties(**{**good, key: bad})
    with pytest.raises(ValueError):
        R1HeliumProperties(**good, T_ref=-1.0)
    broken = tmp_path / "bad.yaml"
    broken.write_text("helium:\n  density_model: ideal_gas\n  specific_gas_constant_J_kgK: 2077.46\n", encoding="utf-8")
    with pytest.raises(KeyError, match="missing helium property key"):
        R1HeliumProperties.from_yaml(broken)
    text_number = tmp_path / "text.yaml"
    text_number.write_text(
        "helium:\n  density_model: ideal_gas\n  specific_gas_constant_J_kgK: 2077.46\n  cp_J_kgK: 5195\n"
        "  viscosity: {a: 3.674e-7, b: 0.7}\n  conductivity: {a: 2.682e-3, b: 1.123e-8, c: 0.71, d: 2e-9}\n",
        encoding="utf-8")
    with pytest.raises(TypeError, match="must be a number"):   # PyYAML은 2e-9를 문자열로 읽는다
        R1HeliumProperties.from_yaml(text_number)
