#!/usr/bin/env python3
"""
Split trajectory script.

Usage:
    python run_split_traj.py [config.yaml]

If no config file is specified, defaults to 'split_traj_settings.yaml' in the current directory.
"""

import argparse
import sys
from pathlib import Path
import yaml

from qm4d4crystal.analysis.split_traj import split_trajectory


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
    
    print(f"Running trajectory split with config: {config_path}")
    print(f"  Trajectory: {traj_path}")
    print(f"  Topology: {top_path}")
    print(f"  Output directory: {output_dir}")
    print(f"  Number of splits: {n_splits}")
    print()
    
    split_trajectory(
        str(traj_path),
        str(top_path),
        n_splits,
        str(output_dir)
    )


if __name__ == "__main__":
    main()
