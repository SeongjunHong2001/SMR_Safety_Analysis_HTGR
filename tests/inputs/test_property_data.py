"""data/properties 입력 검증 (Role 1).

r1_benchmark_properties.yaml의 R1 물성식이 숫자로 읽히는지, graphite.csv·helium.csv의
표 값이 그 식과 일치하는지, 물리적으로 그럴듯한 경향(비열 증가, 조사량에 따른 열전도도
감소, 구간 경계 연속성)을 보이는지 확인한다. 물성 provider 구현 시험이 아니다.
"""
import csv
import math
from functools import lru_cache
from numbers import Number
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
PROPS = ROOT / "data" / "properties"


@lru_cache(maxsize=1)
def _props():
    return yaml.safe_load((PROPS / "r1_benchmark_properties.yaml").read_text(encoding="utf-8"))


def _rows(name):
    return list(csv.DictReader((PROPS / name).open(encoding="utf-8")))


# ---- R1 식 계산 (CSV 생성과 시험이 같은 함수를 쓴다) ----
def helium(T, P, props=None):
    he = (props or _props())["helium"]
    rho = P / (he["specific_gas_constant_J_kgK"] * T)
    mu = he["viscosity"]["a"] * T ** he["viscosity"]["b"]
    c = he["conductivity"]
    k = c["a"] * (1 + c["b"] * P) * T ** (c["c"] * (1 - c["d"] * P))
    return {"density": rho, "cp": he["cp_J_kgK"], "viscosity": mu, "conductivity": k}


def _clamp(x, lo, hi):
    return min(max(x, lo), hi)


def graphite_cp(T, props=None):
    """R1은 비열 식에 온도 범위 제한을 두지 않았다 (300–1900 K 제한은 열전도도에만 해당)."""
    cp = (props or _props())["graphite_H451"]["cp"]
    c = cp["c"]
    return cp["multiplier"] * (c[0] + c[1] * T + c[2] / T + c[3] / T ** 2 + c[4] / T ** 3 + c[5] / T ** 4)


def _interp_fluence(table, fluence, fn):
    fl = [row["fluence"] for row in table]
    fluence = _clamp(fluence, fl[0], fl[-1])
    for lo, hi in zip(table, table[1:]):
        if lo["fluence"] <= fluence <= hi["fluence"]:
            w = (fluence - lo["fluence"]) / (hi["fluence"] - lo["fluence"])
            return (1 - w) * fn(lo) + w * fn(hi)
    return fn(table[0])


def h451_conductivity(T, fluence=0.0, props=None):
    k = (props or _props())["graphite_H451"]["conductivity"]
    T = _clamp(T, *k["temperature_clamp_K"])
    return _interp_fluence(k["table"], fluence, lambda row: sum(a * T ** n for n, a in enumerate(row["a"])))


def _segment_value(segments, T):
    for seg in segments:
        lo, hi = seg["T_K"]
        if lo <= T <= hi:
            if "poly" in seg:
                return sum(a * T ** n for n, a in enumerate(seg["poly"]))
            A, B = seg["exp"]
            return A * math.exp(B * T)
    raise ValueError(T)


def g2020_conductivity(T, fluence=0.0, props=None):
    k = (props or _props())["graphite_2020"]["conductivity"]
    T = _clamp(T, *k["temperature_clamp_K"])
    return _interp_fluence(k["table"], fluence, lambda row: _segment_value(row["segments"], T))


# ---- 시험 ----
def _walk(node, path=()):
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, path + (str(key),))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _walk(value, path + (str(i),))
    else:
        yield path, node


def test_yaml_numbers_are_parsed_as_numbers():
    """1e5처럼 부호 없는 지수는 PyYAML이 문자열로 읽는다. 숫자처럼 보이는 문자열이 없어야 한다."""
    for path, value in _walk(_props()):
        if isinstance(value, str):
            try:
                float(value)
            except ValueError:
                continue
            raise AssertionError(f"숫자가 문자열로 읽힘: {path} = {value!r}")
    he = _props()["helium"]
    assert math.isclose(he["specific_gas_constant_J_kgK"],
                        he["universal_gas_constant_J_molK"] / he["molar_mass_kg_mol"], rel_tol=1e-5)


def test_graphite_csv_matches_r1_unirradiated_correlations():
    rows = _rows("graphite.csv")
    assert {r["material_grade"] for r in rows} == {"H-451", "2020"}
    for r in rows:
        T = float(r["temperature_K"])
        grade = r["material_grade"]
        k = h451_conductivity(T) if grade == "H-451" else g2020_conductivity(T)
        rho = _props()["graphite_H451" if grade == "H-451" else "graphite_2020"]["density_kg_m3"]
        assert math.isclose(float(r["conductivity_W_mK"]), k, rel_tol=1e-4), (grade, T)
        assert math.isclose(float(r["cp_J_kgK"]), graphite_cp(T), rel_tol=1e-4), (grade, T)
        assert float(r["density_kg_m3"]) == rho
        assert r["source"].startswith("R1")


def test_helium_csv_matches_r1_table_aiv13():
    rows = _rows("helium.csv")
    assert len({r["pressure_Pa"] for r in rows}) >= 3
    for r in rows:
        ref = helium(float(r["temperature_K"]), float(r["pressure_Pa"]))
        for col, key in (("density_kg_m3", "density"), ("cp_J_kgK", "cp"),
                         ("viscosity_Pa_s", "viscosity"), ("conductivity_W_mK", "conductivity")):
            assert math.isclose(float(r[col]), ref[key], rel_tol=1e-4), (r, col)


def test_helium_density_matches_operation_distribution():
    """operation.yaml 정격 분포의 노심 입구 밀도(R = 2077.3)와 R1 이상기체(R = 2077.46)가 0.1 % 안에서 같다."""
    op = yaml.safe_load((ROOT / "config" / "operation.yaml").read_text(encoding="utf-8"))
    core_inlet = next(s for s in op["initial_distribution"]["stations"] if s["point"] == "core_inlet")
    rho = helium(core_inlet["T_K"], core_inlet["p_Pa"])["density"]
    assert math.isclose(rho, core_inlet["rho_kg_m3"], rel_tol=1e-3)


def test_graphite_trends_are_physical():
    temps = range(300, 1901, 100)
    cps = [graphite_cp(T) for T in temps]
    assert all(a < b for a, b in zip(cps, cps[1:])), "흑연 비열은 온도에 따라 증가"
    assert 600 < cps[0] < 800 and 1800 < cps[-1] < 2200

    for T in (600, 1000, 1400):
        h451 = [h451_conductivity(T, f) for f in (0.0, 0.2, 0.5, 1.0, 3.0, 8.0)]
        assert all(a > b for a, b in zip(h451, h451[1:])), ("H-451 k는 조사량에 따라 감소", T, h451)
    for T in (500, 900, 1500):
        g2020 = [g2020_conductivity(T, f) for f in (0.0, 1.0, 4.0, 20.0)]
        assert g2020[0] > g2020[-1], ("2020 k는 조사량에 따라 감소", T, g2020)


def test_2020_piecewise_fits_are_nearly_continuous():
    """R1 구간식은 경계에서 완전히 이어지지 않는다. 현재 최대 불연속을 기록해 오타를 잡는다."""
    for row in _props()["graphite_2020"]["conductivity"]["table"]:
        segs = row["segments"]
        for left, right in zip(segs, segs[1:]):
            T = left["T_K"][1]
            assert T == right["T_K"][0]
            a, b = _segment_value([left], T), _segment_value([right], T)
            assert abs(a - b) / max(a, b) < 0.08, (row["fluence"], T, a, b)


def test_fixed_and_structural_properties_positive():
    p = _props()
    assert all(v > 0 for v in p["fixed_properties_ex2ab"]["conductivity_W_mK"].values())
    for name in ("alloy_800H", "pressure_vessel_SA533B", "ceramic_tile", "thermal_insulation"):
        m = p[name]
        assert m["density_kg_m3"] > 0 and m["cp_J_kgK"] > 0 and 0 < m["emissivity"] <= 1, name
    for name in ("ceramic_tile", "thermal_insulation"):
        poly = p[name]["conductivity_poly"]
        for T in (300, 800, 1500):
            assert sum(a * T ** n for n, a in enumerate(poly)) > 0, (name, T)
