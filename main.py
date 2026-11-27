"""실행 진입점.

기본(인자 없음 또는 --demo)은 실제 원자로가 아닌 3변수 Demo DAE를 계산한다 (output/results/demo_dae.csv).
--scenario steady_state 는 config에서 실제 Role 1–4 모델로 S0 정상상태를 조립·풀이한다 (system/simulation.py).
미정 입력이 있으면 파일·키·시나리오를 적고 종료 코드 2로 실패한다.
"""
import argparse
import csv
import sys
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
                   newton_options=options, residual_scale=system.residual_scales(),
                   **adaptive_options)
    result = solver.solve(system.residual,
                          t_span=(config['t_start'], config['t_end']),
                          y0=system.registry.initial_state(),
                          differential_mask=system.registry.differential_mask(),
                          preflight=system.validate_for_solve)
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

def run_scenario_cli(args):
    """--scenario 실행. 성공 0, 미정 입력·미구현 2."""
    from models._inputs import ModelInputMissing
    from system.simulation import run_scenario
    try:
        report = run_scenario(args.scenario, args.condition, assumptions=args.assumptions or (),
                              config_dir=args.config_dir, hold=not args.no_hold, output_dir=args.output_dir)
    except (ModelInputMissing, NotImplementedError, FileNotFoundError, KeyError, ValueError) as exc:
        print(f"[{args.scenario}] not run — missing input or implementation:", file=sys.stderr)
        print(exc, file=sys.stderr)
        return 2
    print("  summary (SI):")
    for key, value in report.summary.items():
        print(f"    {key}: {value:.6g}" if isinstance(value, float) else f"    {key}: {value}")
    for key, row in report.checks.get("acceptance", {}).items():
        print(f"  acceptance {key}: {row['computed']:.6g} <= {row['limit']:.6g} -> {'PASS' if row['pass'] else 'FAIL'}")
    meta = report.steady.built.describe()
    if meta["assumptions"] or meta["calibration"]:
        print(f"  [가정] {len(meta['assumptions'])} provisional input(s), {len(meta['calibration'])} calibrated K value(s) "
              "used — see *_meta.json")
    return 0


def run_calibration_cli(args):
    """--calibrate: S0 system_loop 정상해에서 연결부 K 보정 후 파일로 쓴다. 성공 0, 실패 2."""
    from models._inputs import ModelInputMissing
    from system.builder import load_inputs
    from system.calibration import CalibrationError, calibrate_pressure_losses, write_calibration
    try:
        inputs = load_inputs(args.scenario or 'steady_state', args.condition or 'system_loop',
                             config_dir=args.config_dir, assumptions=args.assumptions or ())
        print(f"[K calibration] {inputs.scenario_name}/{inputs.condition_name}")
        cal = calibrate_pressure_losses(inputs)
    except (ModelInputMissing, NotImplementedError, FileNotFoundError, KeyError, ValueError, CalibrationError) as exc:
        print("[K calibration] failed:", file=sys.stderr)
        print(exc, file=sys.stderr)
        return 2
    print(f"  {'segment':26s} {'target Pa':>10s} {'irrev. Pa':>10s}  K_form increment")
    for seg in cal.segments:
        ks = ', '.join(f"{j}={cal.K[j]:+.4g}" for j in seg.k_junctions)
        print(f"  {seg.name:26s} {seg.design_Pa:10.0f} {cal.losses[seg.name]['irreversible']:10.1f}  {ks}")
    sm = cal.summary
    print(f"  circulator rise {sm['circulator_pressure_rise_Pa']:.0f} Pa = irreversible {sm['sum_irreversible_Pa']:.0f}"
          f" + reversible (gravity, acceleration) {sm['reversible_net_Pa']:.0f}  (PSID 91000)")
    write_calibration(cal, args.calibration_out)
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HTGR primary-system code (Role 5 runner)')
    parser.add_argument('--demo', action='store_true', help='Run dummy heater/cooler DAE (default)')
    parser.add_argument('--solver-config', default='config/solver.yaml', help='demo solver settings')
    parser.add_argument('--scenario', help='scenario file name under <config-dir>/scenarios (e.g. steady_state)')
    parser.add_argument('--condition', help='conditions.<name> (default: scenario default_condition)')
    parser.add_argument('--config-dir', default=None, help='directory with geometry.yaml, operation.yaml, scenarios/')
    parser.add_argument('--assumptions', action='append',
                        help='[가정] file filling null inputs (repeatable), e.g. config/assumptions/provisional_s0.yaml')
    parser.add_argument('--calibration', action='append',
                        help='calibration file (e.g. config/calibration/s0_pressure_loss_K.yaml), repeatable')
    parser.add_argument('--calibrate', action='store_true',
                        help='calibrate junction K to the PSID loss allocation (S0 system_loop) and write a file')
    parser.add_argument('--calibration-out', default='config/calibration/s0_pressure_loss_K.yaml')
    parser.add_argument('--output-dir', default=None, help='default output/results')
    parser.add_argument('--no-hold', action='store_true', help='skip the hold (no-disturbance) DASSL check')
    args = parser.parse_args()
    if args.scenario and args.demo:
        parser.error('--demo and --scenario are exclusive')
    args.assumptions = (args.assumptions or []) + (args.calibration or [])
    if args.calibrate:
        sys.exit(run_calibration_cli(args))
    if args.scenario:
        sys.exit(run_scenario_cli(args))
    run_demo(args.solver_config)
