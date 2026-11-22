"""기본 실행은 실제 원자로가 아닌 3변수 Demo DAE를 계산한다."""
import argparse
import csv
from pathlib import Path
import yaml
from examples.demo_system import make_demo_system
from solvers import DASSL
from solvers.newton import NewtonOptions


def run_demo(config_path='config/solver.yaml', write_output=True):
    root = Path(__file__).resolve().parent
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = root / config_path
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))['solver']
    if config.get('method', 'DASSL') not in {'DASSL', 'BDF_DAE_Prototype'}:
        raise ValueError(f"Unsupported solver method: {config['method']}")
    system = make_demo_system()
    options = NewtonOptions(max_iter=config['newton_max_iter'],
                            atol=config['newton_atol'], rtol=config['newton_rtol'])
    adaptive_options = {
        key: config[key] for key in ('rtol', 'atol', 'min_step', 'max_step', 'max_steps', 'max_retries')
        if key in config
    }
    solver = DASSL(dt=config['dt'], max_order=config['max_order'],
                   newton_options=options, residual_scale=[50.0, 50.0, 50.0],
                   **adaptive_options)
    result = solver.solve(system.residual,
                          t_span=(config['t_start'], config['t_end']),
                          y0=system.registry.initial_state(),
                          differential_mask=system.registry.differential_mask())
    if write_output:
        out = root / 'output' / 'results'
        out.mkdir(parents=True,exist_ok=True)
        filepath = out / 'demo_dae.csv'
        with filepath.open('w',encoding='utf-8',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['time_s','demo.hot_T_K','demo.cold_T_K','demo.q_dot_W'])
            for t, state in zip(result.t,result.y):
                writer.writerow([t,*state])
        print(f"[DEMO only] wrote: {filepath}")
    print(f"[DEMO only] Method: {result.method}")
    print(f"Steps: {len(result.orders)} accepted, {result.rejected_steps} rejected; "
          f"max order used: {max(result.orders)}")
    print(f"Final t: {result.t[-1]:.3f}s; final y: {result.y[-1]}")
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HTGR Role 5 starter (DEMO DAE only)')
    parser.add_argument('--demo',action='store_true', help='Run dummy heater/cooler DAE (default)')
    parser.add_argument('--solver-config',default='config/solver.yaml')
    args = parser.parse_args()
    run_demo(args.solver_config)
