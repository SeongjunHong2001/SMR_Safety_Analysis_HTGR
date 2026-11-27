"""Role 2 장치 모델의 입력 확인 도구 (비공개).

설정 없이 만든 장치(기존 ``Piping()`` 같은 호환 생성)는 등록할 때 ``ModelInputMissing``을 낸다.
``NotImplementedError``의 하위 클래스라 "실제 입력 없이 생산 경로가 실행되지 않는다"는 공통 계약 시험과 맞고,
메시지에 파일 경로와 누락 항목을 적는다.
"""
import math

GEOMETRY = "config/geometry.yaml"


class ModelInputMissing(NotImplementedError):
    """필수 입력 누락 (파일 경로·키를 메시지에 포함)."""


def require(mapping, key, where, *, positive=True):
    """mapping[key]가 유한한 숫자(positive면 > 0)인지 확인하고 float로 돌려준다."""
    value = mapping.get(key) if hasattr(mapping, "get") else None
    if value is None:
        raise ModelInputMissing(f"{GEOMETRY}: missing required input {where}.{key}")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{GEOMETRY}: {where}.{key} must be a finite number, got {value!r}")
    if positive and value <= 0:
        raise ValueError(f"{GEOMETRY}: {where}.{key} must be > 0, got {value!r}")
    return float(value)
