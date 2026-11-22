# SMR Safety Analysis · HTGR

고온가스로(HTGR)의 정상상태와 냉각재 유량 상실(LOFC) 과도현상을 해석하기 위한 Python 프로젝트다. 각 구성요소가 상태변수와 방정식을 정의하면, 계통 전체의 미분대수방정식(DAE) `F(t, y, ydot) = 0`을 조립해 시간에 따른 변화를 계산한다.

현재 실행할 수 있는 예제는 두 열용량 사이의 열전달이다. 이 예제로 상태 등록부터 초기조건 계산, 시간 적분, 결과 저장까지 확인할 수 있다. 실제 원자로를 대상으로 한 계산에는 노형별 구성요소 모델과 입력 데이터가 필요하다. [개발 구현 현황](docs/implementation-status.md)

## 설치와 첫 실행

Python 3.9 이상이 필요하다. 저장소를 내려받은 뒤 프로젝트 최상위 폴더에서 실행한다.

```bash
git clone https://github.com/SeongjunHong2001/SMR_Safety_Analysis_HTGR.git
cd SMR_Safety_Analysis_HTGR
python -m pip install -r requirements.txt
python main.py --demo
```

데모는 고온부와 저온부의 온도를 미분변수로, 두 지점 사이의 열유량을 대수변수로 둔다. 기본 설정은 0초부터 5초까지 적분한다. 터미널에는 수용·거절된 스텝 수와 사용한 최고 차수가 표시된다. 초기 시점과 이후 수용된 각 시점의 결과는 `output/results/demo_dae.csv`에 저장된다.

| CSV 열 | 의미 |
| --- | --- |
| `time_s` | 시뮬레이션 시간(s) |
| `demo.hot_T_K`, `demo.cold_T_K` | 고온부와 저온부 온도(K) |
| `demo.q_dot_W` | 고온부에서 저온부로 전달되는 열유량(W) |

기본 적분 구간, 첫 스텝 후보 간격(`dt`), 오차 허용값(`rtol`, `atol`)은 [솔버 설정](config/solver.yaml)에서 확인할 수 있다. 다른 설정 파일을 사용하려면 `python main.py --demo --solver-config 경로/설정.yaml`을 실행한다. 전체 테스트는 `python -m pytest -q`로 실행한다.

## Python에서 솔버 사용

`DASSL`은 지수 1 DAE를 적분하는 Python 솔버다. 1~5차 후진차분공식(BDF)과 뉴턴법을 사용하며, 가중 RMS 오차에 따라 스텝 간격과 차수를 조정한다.

```python
from examples.demo_system import make_demo_system
from solvers import DASSL, solve_steady_state

system = make_demo_system()
registry = system.registry
result = DASSL(dt=0.1, rtol=1e-6, atol=1e-8).solve(
    system.residual,
    t_span=(0.0, 5.0),
    y0=registry.initial_state(),
    differential_mask=registry.differential_mask(),
)

y_equilibrium, ydot_equilibrium = solve_steady_state(
    system.residual, 0.0, registry.initial_state()
)
```

`differential_mask`를 주면 미분변수의 초깃값을 유지하면서 잔차식을 만족하는 대수변수와 미분변수의 초기 변화율을 구한다. 잔차식이 대수변수의 시간미분에 의존하지 않을 때 사용할 수 있다. 일치하는 초기조건을 알고 있다면 마스크를 생략하고 `y0`와 `ydot0`를 함께 전달한다. `solve_steady_state`는 모든 변화율을 0으로 두고 평형 상태를 구한다.

`result.t`, `result.y`, `result.ydot`에는 초기점을 포함한 시각, 상태, 변화율이 담긴다. `result.orders[i]`, `result.error_norms[i]`, `result.newton_iterations[i]`는 `t[i]`에서 `t[i+1]`까지 수용된 스텝의 정보다. 해석적 야코비안을 사용할 때는 `∂F/∂y + cj·∂F/∂ydot`를 반환하는 `jacobian(t, y, ydot, cj)`를 `solve`에 전달한다. 기본값은 유한차분이다.

## 관련 문서

- [개발 구현 현황](docs/implementation-status.md)에서 기능별 구현 범위와 테스트 결과, 실제 계통 해석에 필요한 작업을 확인할 수 있다.
- [계통 조립과 시간 적분 구조](docs/architecture.md)는 구성요소의 잔차가 전역 솔버로 전달되는 과정을 설명한다.
- 새 구성요소를 추가할 때는 [인터페이스 규약](docs/interface_spec.md)의 상태 등록 방법, 단위, 포트 부호 규약을 따른다.
- [미확정 물리 가정](docs/assumptions.md)에는 대상 노형, 공간 차원, 출력 모델 등 결정이 필요한 항목이 정리돼 있다.
