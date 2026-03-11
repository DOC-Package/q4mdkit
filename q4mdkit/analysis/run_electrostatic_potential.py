#!/usr/bin/env python3
"""
Compute displacement-induced electrostatic potential fluctuations.

This script calculates:
    φ_a^disp(t) = φ_a(t) - φ_a^env(t)

where:
    φ_a(t) = Σ_{b∈MM} Q_b / ||r_a(t) - R_b(t)||
    φ_a^env(t) = Σ_{b∈MM} Q_b / ||r_a^ref(t) - R_b(t)||

Usage:
    python run_electrostatic_potential.py config.yaml
    python run_electrostatic_potential.py --dcd traj.dcd --top system.pdb --ref ref.gen --qm-index qmatoms --mm-charges PCcharges.dat
"""

import argparse
from pathlib import Path
from q4mdkit.analysis.electrostatic_potential import (
    run_potential_analysis,
    analyze_electrostatic_potential,
)

# Default config file names to search for
DEFAULT_CONFIG_NAMES = [
    'config.yaml',
    'config.yml',
    'expot_config.yaml',
    'electrostatic_potential.yaml',
    'settings.yaml',
]


def find_config_file() -> str | None:
    """Search for a config file in the current directory."""
    for name in DEFAULT_CONFIG_NAMES:
        if Path(name).exists():
            return name
    return None


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description='Compute displacement-induced electrostatic potentials for QM/MM systems',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example with config file:
    python run_electrostatic_potential.py config.yaml

Example with command-line arguments:
    python run_electrostatic_potential.py \\
        --dcd trajectory.dcd \\
        --top system.pdb \\
        --ref reference.gen \\
        --qm-index qm_atoms.ndx \\
        --mm-charges mm_charges.dat \\
        -o output/potential

Output files:
    {prefix}_phi_disp.dat       - Displacement potential φ_a^disp(t) for each QM atom
    {prefix}_phi_disp_summary.dat - Sum and mean of φ_disp over QM atoms
    {prefix}_phi.dat            - Total potential φ_a(t) (optional)
    {prefix}_phi_env.dat        - Envelope potential φ_a^env(t) (optional)
    {prefix}_r_cm.dat           - Center of mass trajectory
    {prefix}_rotation.dat       - Rotation angles from reference
"""
    )
    
    # Config file as positional argument (optional)
    parser.add_argument('config', nargs='?', help='YAML configuration file')
    
    # Command-line options (used if config not specified)
    parser.add_argument('--dcd', '-d', help='DCD trajectory file')
    parser.add_argument('--top', '-t', help='Topology file (PDB)')
    parser.add_argument('--ref', '-r', help='Reference structure (.gen)')
    parser.add_argument('--qm-index', '-q', help='QM atom index file')
    parser.add_argument('--mm-charges', '-c', help='MM charges file')
    parser.add_argument('--output', '-o', default='potential', help='Output prefix')
    parser.add_argument('--start', type=int, default=0, help='Start frame (0-indexed)')
    parser.add_argument('--end', type=int, default=None, help='End frame (exclusive)')
    parser.add_argument('--stride', type=int, default=1, help='Frame stride')
    parser.add_argument('--charge-column', type=int, default=-1,
                        help='Column index for charge in mm_charges file (-1 = last column)')
    parser.add_argument('--no-individual', action='store_true',
                        help='Skip saving individual atom potentials')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output')
    
    args = parser.parse_args()
    
    # Auto-detect config file if not specified
    config_file = args.config
    if config_file is None:
        config_file = find_config_file()
        if config_file:
            print(f"Auto-detected config file: {config_file}")
    
    if config_file:
        # Run from config file
        run_potential_analysis(config_file)
    else:
        # Validate required arguments
        required = ['dcd', 'top', 'ref', 'qm_index', 'mm_charges']
        missing = [arg for arg in required if getattr(args, arg.replace('-', '_')) is None]
        if missing:
            parser.error(f"Missing required arguments: {', '.join('--' + m.replace('_', '-') for m in missing)}. "
                        f"Either provide a config file or all required arguments.")
        
        analyze_electrostatic_potential(
            dcd_file=args.dcd,
            top_file=args.top,
            ref_gen_file=args.ref,
            qm_index_file=args.qm_index,
            mm_charges_file=args.mm_charges,
            output_prefix=args.output,
            save_individual=not args.no_individual,
            verbose=not args.quiet,
            start_frame=args.start,
            end_frame=args.end,
            stride=args.stride,
            charge_column=args.charge_column
        )


if __name__ == "__main__":
    main()
