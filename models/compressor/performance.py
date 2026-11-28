"""Role 2 compressor performance provider contract and a coefficient-explicit map.

검증된 MHTGR-350 주순환기 성능곡선은 아직 없다 (보고서[1] M09). 그래서 기본 계수를 두지 않고,
``QuadraticPerformance``는 모든 계수를 생성자에서 필수로 받는다.
"""
from dataclasses import dataclass
from typing import Protocol


class CompressorPerformance(Protocol):
    def evaluate(self, *, speed, mass_flow, inlet_state, geometry):
        """Return pressure rise and shaft work using a validated map."""
        ...


@dataclass(frozen=True)
class QuadraticPerformance:
    """상사 2차식 성능 맵 (계수는 모두 필수, 기본값 없음).

    ΔP = (ρ/ρ_rated) ΔP_rated [c0 α² + c1 α ν − c2 ν|ν|],  α = ω/ω_rated,  ν = V̇/V̇_rated
    변수표 c0_circ_curve, c1_circ_curve, c2_circ_curve와 같은 부호 규약 (mode_circ_curve = affinity_quad).
    정격점(α = ν = 1, ρ = ρ_rated)에서 ΔP = ΔP_rated가 되려면 c0 + c1 − c2 = 1 이어야 한다
    (변수표 예비값 1.3, 0, 0.3은 제조사 곡선 확보 전 대리값).
    축일: ΔP V̇ / efficiency (efficiency 필수).
    """

    rated_speed: float          # rad/s
    rated_mass_flow: float      # kg/s
    rated_density: float        # kg/m³
    rated_pressure_rise: float  # Pa
    c0: float
    c1: float
    c2: float
    efficiency: float

    def __post_init__(self):
        for name in ("rated_speed", "rated_mass_flow", "rated_density", "rated_pressure_rise"):
            if not getattr(self, name) > 0:
                raise ValueError(f"QuadraticPerformance.{name} must be > 0")
        if not 0 < self.efficiency <= 1:
            raise ValueError("efficiency must be in (0, 1]")
        if abs(self.c0 + self.c1 - self.c2 - 1.0) > 1e-9:
            raise ValueError("c0 + c1 - c2 must equal 1 so the map passes through the rated point")

    def evaluate(self, *, speed, mass_flow, inlet_state, geometry):
        rho = float(inlet_state["rho"])
        alpha = float(speed) / self.rated_speed
        vdot = float(mass_flow) / rho
        nu = vdot / (self.rated_mass_flow / self.rated_density)
        dp = rho / self.rated_density * self.rated_pressure_rise * (
            self.c0 * alpha * alpha + self.c1 * alpha * nu - self.c2 * nu * abs(nu))
        return {"pressure_rise": dp, "shaft_work": dp * vdot / self.efficiency}
