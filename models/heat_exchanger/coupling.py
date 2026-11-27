"""Role 4 SG 열 결합.

``SG0Coupling``  SG-0 일정 온도 물 열침원: Q_i = U A_i (T_He,i − T_물) [W, He에서 빠지는 열 양수].
    HeatExchanger 계약대로 (to_helium, to_wall_from_helium, to_wall_from_secondary, to_secondary) =
    (−Q, +Q, −Q, +Q)를 돌려준다. 벽(NoStorageWall)의 순열률은 0이다.
    T_He,i는 FVM 셀 값(1차 상류차분에서 셀 출구 상태)이고, U는 S0에서 He 출구 528.15 K가 되도록 보정한
    등가값이라 셀 평균 대신 셀 값을 쓴 차이를 흡수한다.
    A_i = geometry heat_exchanger.heat_transfer_area_m2 (관 외면, 높이 비례 배분, 팀 결정).
"""
import numpy as np


class SG0Coupling:
    def __init__(self, helium_side, secondary):
        self.helium_side = helium_side
        self.secondary = secondary

    def heat_removed(self, t, y, context):
        fs = self.helium_side.coupling_state(t, y, context)
        U = self.secondary.U(y, context)
        return U * self.helium_side.area_heat * (fs["T_bulk"] - self.secondary.T_water)

    def evaluate(self, t, y, context):
        Q = self.heat_removed(t, y, context)
        return -Q, Q.copy(), -Q, Q.copy()
