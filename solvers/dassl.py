"""선두 계수를 고정한 가변 간격 BDF로 F(t, y, ydot) = 0을 적분한다."""
from __future__ import annotations

from dataclasses import dataclass, field
from math import factorial

import numpy as np

from .bdf import bdf_predictor, divided_differences
from .error_control import error_weights, rms_norm, step_factor
from .initialization import consistent_initialize
from .newton import NewtonFailure, NewtonOptions, newton_solve


class IntegrationFailure(RuntimeError):
    """일치 초기조건을 얻지 못하거나 시간 적분을 계속할 수 없을 때 발생한다."""


@dataclass
class IntegrationResult:
    """t, y, ydot는 초기점을 포함한다. 반복 횟수, 차수, 오차는 수용한 스텝마다 하나씩 기록한다."""

    t: np.ndarray
    y: np.ndarray
    ydot: np.ndarray
    newton_iterations: np.ndarray
    method: str = "DASSL-style adaptive BDF1-5 (Python)"
    orders: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    error_norms: np.ndarray = field(default_factory=lambda: np.empty(0))
    error_test_failures: int = 0
    newton_failures: int = 0

    @property
    def rejected_steps(self):
        """국소 오차 초과 또는 뉴턴 보정 실패로 거절된 스텝 시도 횟수."""
        return self.error_test_failures + self.newton_failures


class DASSL:
    """정칙 지수 1 미분대수방정식을 1~5차 암시적 BDF로 적분한다."""

    def __init__(self, *, dt: float = 0.1, max_order: int = 5,
                 rtol=1e-5, atol=1e-8, min_step=0.0, max_step=np.inf,
                 max_steps=100_000, max_retries=12,
                 newton_options: NewtonOptions | None = None, residual_scale=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be positive and finite")
        if not isinstance(max_order, (int, np.integer)) or not 1 <= max_order <= 5:
            raise ValueError("max_order must be an integer from 1 to 5")
        if not np.isfinite(min_step) or not 0 <= min_step <= max_step or max_step <= 0:
            raise ValueError("step limits must satisfy 0 <= min_step <= max_step and max_step > 0")
        if not isinstance(max_steps, (int, np.integer)) or max_steps < 1:
            raise ValueError("max_steps must be a positive integer")
        if not isinstance(max_retries, (int, np.integer)) or max_retries < 0:
            raise ValueError("max_retries must be a nonnegative integer")
        if any(np.any(~np.isfinite(tol)) or np.any(tol < 0)
               for tol in (np.asarray(rtol, dtype=float), np.asarray(atol, dtype=float))):
            raise ValueError("rtol and atol must be finite and nonnegative")
        self.dt = float(dt)
        self.max_order = max_order
        self.rtol = np.array(rtol, dtype=float)
        self.atol = np.array(atol, dtype=float)
        self.min_step = float(min_step)
        self.max_step = float(max_step)
        self.max_steps = max_steps
        self.max_retries = max_retries
        self.newton_options = newton_options or NewtonOptions()
        self.residual_scale = residual_scale

    def solve(self, residual, *, t_span, y0, ydot0=None, differential_mask=None,
              jacobian=None) -> IntegrationResult:
        """시작 시각에서 종료 시각까지 순방향 또는 역방향으로 적분하고 수용된 해를 반환한다.

        differential_mask가 있으면 미분 상태를 고정하고 대수 상태와 미분 상태의
        초기 변화율을 구한다. 이때 ydot0는 초기 추정값이다. 마스크 없이 y0와
        ydot0를 함께 주면 두 값이 잔차식을 만족하는지 검사한다. 마스크와 ydot0를
        모두 생략하면 모든 초기 상태를 고정하고 변화율만 구한다.

        jacobian(t, y, ydot, cj)는 보정 방정식의 야코비안 ∂F/∂y + cj·∂F/∂ydot를 반환해야 한다.
        """
        t0, tf = (float(x) for x in t_span)
        if tf == t0 or not np.isfinite(tf - t0):
            raise ValueError("t_span must have distinct finite endpoints and a finite duration")
        y = np.asarray(y0, dtype=float).copy()
        if y.ndim != 1 or not y.size or not np.all(np.isfinite(y)):
            raise ValueError("y0 must be a nonempty finite 1-D state vector")
        for name, tolerance in (("atol", self.atol), ("rtol", self.rtol)):
            if tolerance.shape not in ((), y.shape):
                raise ValueError(f"{name} must be a scalar or match y0")
        yp = None if ydot0 is None else np.asarray(ydot0, dtype=float).copy()
        if yp is not None and (yp.shape != y.shape or not np.all(np.isfinite(yp))):
            raise ValueError("ydot0 must be finite and match y0")
        scale = (np.ones_like(y) if self.residual_scale is None
                 else np.asarray(self.residual_scale, dtype=float))
        if scale.shape != y.shape or np.any(scale <= 0) or not np.all(np.isfinite(scale)):
            raise ValueError("residual_scale must be positive, finite and match y0")
        initial_weights = error_weights(y, atol=self.atol, rtol=self.rtol)

        try:
            if differential_mask is not None or yp is None:
                mask = (np.ones(y.shape, dtype=bool) if differential_mask is None
                        else np.asarray(differential_mask, dtype=bool))
                if mask.shape != y.shape:
                    raise ValueError("differential_mask shape must match y0")
                # 변화율의 보정 척도는 상태 허용오차를 첫 스텝의 후보 간격으로 나누어 구한다.
                derivative_weights = initial_weights[mask] / self.dt
                init_scale = np.concatenate((initial_weights[~mask], derivative_weights))
                y, yp = consistent_initialize(
                    residual, t0, y, mask, ydot_guess=yp, options=self.newton_options,
                    residual_scale=scale, correction_scale=init_scale,
                )
            initial_residual = np.asarray(residual(t0, y, yp), dtype=float).ravel()
        except (NewtonFailure, FloatingPointError) as exc:
            raise IntegrationFailure(f"Consistent initialization failed at t={t0:.16g}: {exc}") from exc
        if initial_residual.shape != y.shape:
            raise ValueError("residual size must match y0")
        initial_norm = np.max(np.abs(initial_residual / scale))
        initial_target = self.newton_options.atol + self.newton_options.rtol
        if not np.isfinite(initial_norm) or initial_norm > initial_target:
            raise IntegrationFailure(
                f"Inconsistent initial conditions at t={t0:.16g}; "
                f"scaled residual={initial_norm:.3e}. Supply a differential_mask "
                "or consistent y0 and ydot0."
            )
        error_weights(y, atol=self.atol, rtol=self.rtol)

        times, states, derivatives = [t0], [y], [yp]
        history_times, history_states = [t0], [y]
        iterations, orders, errors = [], [], []
        error_failures = newton_failures = retries = same_order_steps = 0
        order = 1
        direction = 1.0 if tf > t0 else -1.0
        step = min(max(self.dt, self.min_step), self.max_step)

        while times[-1] != tf:
            t = times[-1]
            if len(orders) == self.max_steps:
                raise IntegrationFailure(f"max_steps={self.max_steps} reached at t={t:.16g}")
            remaining = direction * (tf - t)
            t_new = tf if step >= remaining else t + direction * step
            h = t_new - t
            if h == 0.0:
                raise IntegrationFailure(f"Step cannot advance time at t={t:.16g}; step={step:.3e}")

            try:
                y_new, yp_new, error, count = self._attempt_step(
                    residual, jacobian, t_new, history_times, history_states,
                    derivatives[-1], order,
                )
            except (NewtonFailure, FloatingPointError) as exc:
                newton_failures += 1
                failure = exc
                ratio = 0.25
            else:
                if error <= 1.0:
                    times.append(t_new)
                    states.append(y_new)
                    derivatives.append(yp_new)
                    iterations.append(count)
                    orders.append(order)
                    errors.append(error)
                    # 예측과 차수 비교에 쓰는 이력에는 수용된 스텝만 최신순으로 넣는다.
                    history_times.insert(0, t_new)
                    history_states.insert(0, y_new)
                    del history_times[self.max_order + 2:]
                    del history_states[self.max_order + 2:]
                    if t_new == tf:
                        break
                    same_order_steps += 1
                    next_order, ratio = self._next_step(
                        h, order, error, history_times, history_states, same_order_steps,
                    )
                    if next_order != order:
                        same_order_steps = 0
                    order = next_order
                    step = min(self.max_step, max(self.min_step, abs(h) * ratio))
                    retries = 0
                    continue
                error_failures += 1
                failure = None
                ratio = np.clip(step_factor(error, order), 0.1, 0.8)

            # 같은 시각에서 더 짧은 간격으로 다시 시도한다. 연속 두 번 이상 실패하면 차수도 낮춘다.
            retries += 1
            smaller_step = max(self.min_step, abs(h) * ratio)
            if (retries > self.max_retries or smaller_step >= abs(h) or
                    t + direction * smaller_step == t_new):
                description = str(failure) if failure is not None else f"weighted LTE={error:.3e}"
                raise IntegrationFailure(
                    f"Cannot accept a step from t={t:.16g} after {retries} failed "
                    f"attempt(s); step={abs(h):.3e}, order={order}; {description}"
                ) from failure
            step = smaller_step
            same_order_steps = 0
            if retries >= 2:
                order = max(1, order - 1)

        return IntegrationResult(
            t=np.asarray(times), y=np.vstack(states), ydot=np.vstack(derivatives),
            newton_iterations=np.asarray(iterations, dtype=int),
            orders=np.asarray(orders, dtype=int), error_norms=np.asarray(errors),
            error_test_failures=error_failures, newton_failures=newton_failures,
        )

    def _correct(self, residual, jacobian, t, prediction, slope, cj, previous):
        if not np.all(np.isfinite(prediction)) or not np.all(np.isfinite(slope)) or not np.isfinite(cj):
            raise FloatingPointError("Non-finite BDF predictor or coefficient")

        # 미지수는 예측 상태의 보정량이다. 같은 보정량으로 상태와 시간미분을 함께 갱신한다.
        def equations(delta):
            return residual(t, prediction + delta, slope + cj * delta)

        def iteration_matrix(delta):
            return jacobian(t, prediction + delta, slope + cj * delta, cj)

        weights = error_weights(prediction, previous=previous, atol=self.atol, rtol=self.rtol)
        # 짧은 스텝에서 cj가 커지므로 수치 미분 시 시간미분의 변화를 제한하되,
        # 상태에 가한 변화가 부동소수점 반올림으로 사라지지 않도록 하한을 둔다.
        state_magnitude = np.maximum(1.0, np.abs(prediction))
        derivative_magnitude = np.maximum(1.0, np.abs(slope))
        fd_scales = np.maximum(
            8 * np.sqrt(np.finfo(float).eps) * np.abs(prediction),
            np.minimum(state_magnitude, derivative_magnitude / abs(cj)),
        )
        solution = newton_solve(
            equations, np.zeros_like(prediction), options=self.newton_options,
            residual_scale=self.residual_scale, jacobian=iteration_matrix if jacobian is not None else None,
            correction_scale=weights, fd_scales=fd_scales,
        )
        delta = solution.x
        return prediction + delta, slope + cj * delta, delta, solution.iterations

    def _attempt_step(self, residual, jacobian, t, past_times, past_states, rate, order):
        previous = past_states[0]
        h = t - past_times[0]
        if len(past_times) > 1:
            prediction, slope, cj, error_factor = bdf_predictor(t, past_times, past_states, order)
            y, yp, correction, count = self._correct(
                residual, jacobian, t, prediction, slope, cj, previous,
            )
            error_vector = error_factor * correction
        else:
            # 첫 스텝에는 보간 이력이 부족하다. 1차 BDF 한 스텝과 두 반 스텝을 비교해
            # 대수 상태의 초기 변화율 없이 오차를 추정하고, 수용 시에는 한 스텝의 해를 쓴다.
            y, yp, correction, count = self._correct(
                residual, jacobian, t, previous, np.zeros_like(rate), 1.0 / h, previous,
            )
            midpoint = past_times[0] + 0.5 * h
            if midpoint == past_times[0] or midpoint == t:
                # 중점을 별도 시각으로 표현할 수 없으면 초기 기울기에 따른 변화량과 비교한다.
                error_vector = 0.5 * (correction - h * rate)
            else:
                half, _, first_delta, first_count = self._correct(
                    residual, jacobian, midpoint, previous, np.zeros_like(rate),
                    1.0 / (midpoint - past_times[0]), previous,
                )
                _, _, second_delta, second_count = self._correct(
                    residual, jacobian, t, half, np.zeros_like(rate), 1.0 / (t - midpoint), half,
                )
                # 대수 상태는 두 계산의 끝값이 같을 수 있다. 두 반 스텝의 변화량 차이로
                # 중점 부근의 곡률까지 반영해 이런 경우에도 오차를 감지한다.
                doubling_error = np.abs(first_delta + second_delta - correction)
                midpoint_curvature = np.abs(second_delta - first_delta)
                error_vector = 2.0 * np.maximum(doubling_error, midpoint_curvature)
                count += first_count + second_count

        if not np.all(np.isfinite(y)) or not np.all(np.isfinite(yp)):
            raise FloatingPointError("Non-finite BDF corrector")
        weights = error_weights(y, previous=previous, atol=self.atol, rtol=self.rtol)
        return y, yp, rms_norm(error_vector / weights), count

    def _next_step(self, h, order, error, times, states, same_order_steps):
        next_order = order
        factor = step_factor(error, order)
        if same_order_steps >= order + 1:
            weights = error_weights(states[0], previous=states[1], atol=self.atol, rtol=self.rtol)
            for candidate in (order - 1, order + 1):
                if not 1 <= candidate <= self.max_order or len(times) < candidate + 2:
                    continue
                nodes = (np.asarray(times[:candidate + 2]) - times[0]) / h
                leading_term = divided_differences(nodes, states[:candidate + 2])[-1]
                # 최고차 분할차분으로 간격 h인 등간격 격자에서 후보 차수의 국소 절단 오차를 추정한다.
                estimate = rms_norm(factorial(candidate) * leading_term / weights)
                proposed = step_factor(estimate, candidate)
                if proposed > 1.2 * factor:
                    factor, next_order = proposed, candidate
        max_growth = 1.5 if next_order == order else 1.25
        return next_order, np.clip(factor, 0.5, max_growth)


# 초기 버전의 호출 코드가 쓰는 클래스 이름을 별칭으로 유지한다.
DASSLPrototype = DASSL
