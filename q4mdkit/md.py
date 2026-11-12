"""Common input and restart preparation for explicit MM and QM/MM scripts."""

import argparse
import math
from pathlib import Path

from q4mdkit.qmmm.xml import read_restart_box, remove_montecarlo_xml


def prepare_md_run(theory_config, md_config, ensemble, statefile=None, argv=(),
                   description=None):
    """Validate a run and return the statefile for the MDConfig runner.

    The caller sets its working directory and creates the Fragment/theory
    explicitly after this function. No ASH imports occur here. Pass
    ``argv=None`` in a command-line script to enable ``--validate-only`` and
    ``--simulation-time``; the default consumes no process arguments.

    Restart cell vectors retain State XML precision. Constant-volume runs
    use a local ``initial_state.xml`` with barostat parameters removed.
    The source State is preserved. Optional ``output.protect_existing``
    prevents concatenating independent trajectories into existing outputs.
    """
    if ensemble not in ('nvt', 'npt', 'nve'):
        raise ValueError("ensemble must be 'nvt', 'npt', or 'nve'")
    title = description or f'{ensemble.upper()} MD'
    parser = argparse.ArgumentParser(description=title)
    parser.add_argument('--validate-only', action='store_true',
                        help='validate inputs without importing ASH or running MD')
    parser.add_argument('--simulation-time', type=float, metavar='PS',
                        help='override the simulation duration in ps')
    args = parser.parse_args(argv)
    if args.simulation_time is not None:
        if not math.isfinite(args.simulation_time) or args.simulation_time < md_config.timestep:
            parser.error('--simulation-time must be finite and at least one timestep')
        setattr(md_config, f'{ensemble}_time', args.simulation_time)

    atom_count = theory_config.validate_inputs()
    if statefile is not None:
        statefile = str(statefile)
        vectors = read_restart_box(statefile, expected_atoms=atom_count)
        if theory_config.periodic:
            theory_config.pbc_vectors = vectors
    if args.validate_only:
        print(f'{title} inputs are consistent.')
        if statefile is not None:
            print(f'Initial state: {statefile}')
        raise SystemExit(0)

    output = Path(md_config.output_dir)
    if md_config.protect_existing:
        existing = [path for path in (output / f'{ensemble}.dcd', Path('OpenMM_MD_final_state.xml'))
                    if path.exists()]
        if existing:
            parser.error('Move existing results before starting a new run: '
                         + ', '.join(str(path) for path in existing))
    output.mkdir(parents=True, exist_ok=True)
    if statefile is not None and ensemble in ('nvt', 'nve'):
        destination = 'initial_state.xml'
        if Path(statefile).resolve() == Path(destination).resolve():
            raise ValueError('Restart source must differ from the local initial_state.xml copy')
        remove_montecarlo_xml(statefile, destination)
        return destination
    return statefile
