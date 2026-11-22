import numpy as np
from examples.demo_system import make_demo_system
from main import run_demo
from solvers.dassl import DASSLPrototype
from solvers.initialization import consistent_initialize

def test_demo_system_and_conservation():
    sys=make_demo_system()
    reg=sys.registry
    assert reg.size == 3
    assert reg.differential_mask().tolist()==[True,True,False]
    y0,yp0=consistent_initialize(sys.residual,0.,reg.initial_state(),reg.differential_mask())
    assert np.linalg.norm(sys.residual(0,y0,yp0),ord=np.inf)<1e-6
    sol=DASSLPrototype(dt=0.05,max_order=2,residual_scale=[50.,50.,50.]).solve(
        sys.residual,t_span=(0.,4.),y0=y0,ydot0=yp0)
    # 대수변수: q=5*(Th-Tc)
    assert np.allclose(sol.y[:,2],5*(sol.y[:,0]-sol.y[:,1]),atol=1e-6)
    # 각 수용단계의 동일한 시간에서 전역 DAE가 만족되는지
    for t,y,yp in zip(sol.t[1:],sol.y[1:],sol.ydot[1:]):
        assert np.linalg.norm(sys.residual(t,y,yp),ord=np.inf)<1e-5


def test_default_demo_matches_analytic_temperature_solution():
    result = run_demo(write_output=False)
    equilibrium = np.array([300 + 50 / 3 + 10, 300 + 50 / 3])
    matrix = np.array([[-5 / 20, 5 / 20], [5 / 40, -(5 + 3) / 40]])
    rates, modes = np.linalg.eig(matrix)
    amplitudes = np.linalg.solve(modes, np.array([320.0, 300.0]) - equilibrium)
    exact = equilibrium + (np.exp(np.outer(result.t, rates)) * amplitudes) @ modes.T

    assert result.t[-1] == 5.0
    assert result.orders.max() == 5
    assert np.max(np.abs(result.y[:, :2] - exact)) < 1e-3
    exact_heat_flow = 5 * (exact[:, 0] - exact[:, 1])
    assert np.max(np.abs(result.y[:, 2] - exact_heat_flow)) < 1e-3
