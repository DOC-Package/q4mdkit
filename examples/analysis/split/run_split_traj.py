#!/usr/bin/env python3
"""
Split trajectory and energy file script.

Usage:
    python run_split_traj.py [config.yaml]

If no config file is specified, defaults to 'split_traj_settings.yaml' in the current directory.
"""

import argparse
import sys
from pathlib import Path
import yaml

from qm4d4crystal.analysis.split_traj import split_trajectory
from qm4d4crystal.analysis.split_energy import split_energy_file


def main():
    parser = argparse.ArgumentParser(
        description="Split DCD trajectory into multiple parts"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="split_traj_settings.yaml",
        help="Path to YAML configuration file (default: split_traj_settings.yaml)"
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)
    
    # Load configuration
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    base_dir = config_path.parent
    
    # Parse paths (relative to config file)
    traj_path = base_dir / config.get('trajectory', 'prod.dcd')
    top_path = base_dir / config.get('topology', 'pentacene.pdb')
    output_dir = traj_path.parent  # Output to the same directory as trajectory
    n_splits = config.get('n_splits', 2)
    
    # Frame range settings
    start_frame = config.get('start_frame', 0)
    n_frames = config.get('n_frames', None)
    
    # Parse energy file settings (optional)
    energy_cfg = config.get('energy', {})
    energy_enabled = energy_cfg.get('enabled', False)
    
    print(f"Running trajectory split with config: {config_path}")
    print(f"  Trajectory: {traj_path}")
    print(f"  Topology: {top_path}")
    print(f"  Output directory: {output_dir}")
    print(f"  Number of splits: {n_splits}")
    print(f"  Start frame: {start_frame}")
    print(f"  N frames: {n_frames if n_frames else 'all'}")
    print()
    
    # Split trajectory
    split_trajectory(
        trajectory_path=str(traj_path),
        topology_path=str(top_path),
        n_splits=n_splits,
        output_dir=str(output_dir),
        start_frame=start_frame,
        n_frames=n_frames,
    )
    
    # Split energy file if enabled
    if energy_enabled:
        energy_file = energy_cfg.get('file', 'qm_energy.dat')
        energy_path = str(base_dir / energy_file)
        energy_interval = energy_cfg.get('energy_interval', 1)
        sample_interval = energy_cfg.get('sample_interval', 1)
        
        print(f"\nSplitting energy file: {energy_path}")
        print(f"  Energy interval: {energy_interval} fs")
        print(f"  Sample interval: {sample_interval} fs")
        print()
        
        split_energy_file(
            energy_path=energy_path,
            n_splits=n_splits,
            output_dir=str(output_dir),
            energy_interval=energy_interval,
            sample_interval=sample_interval,
            start_frame=start_frame,
            n_frames=n_frames,
        )


if __name__ == "__main__":
    main()
