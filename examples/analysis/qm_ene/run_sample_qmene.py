#!/usr/bin/env python3
"""
Sample QM energy file script.

Usage:
    python run_sample_qmene.py [config.yaml]

If no config file is specified, defaults to 'sample_qmene_settings.yaml' in the current directory.
"""

import argparse
import sys
from pathlib import Path
import yaml

from qm4d4crystal.analysis.sample_energy import sample_energy_file


def main():
    parser = argparse.ArgumentParser(
        description="Sample QM energy file to match trajectory frames"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="sample_qmene_settings.yaml",
        help="Path to YAML configuration file (default: sample_qmene_settings.yaml)"
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
    
    # Parse paths and settings
    energy_file = config.get('energy_file', 'qm_energy.dat')
    energy_path = base_dir / energy_file
    
    output_file = config.get('output_file', None)
    if output_file:
        output_path = str(base_dir / output_file)
    else:
        output_path = None
    
    energy_interval = config.get('energy_interval', 1)
    sample_interval = config.get('sample_interval', 1)
    start_frame = config.get('start_frame', 0)
    n_frames = config.get('n_frames', None)
    
    print(f"Running energy file sampling with config: {config_path}")
    print(f"  Energy file: {energy_path}")
    print(f"  Energy interval: {energy_interval} fs")
    print(f"  Sample interval: {sample_interval} fs")
    print(f"  Start frame: {start_frame}")
    print(f"  N frames: {n_frames if n_frames else 'all'}")
    print()
    
    sample_energy_file(
        energy_path=str(energy_path),
        output_path=output_path,
        energy_interval=energy_interval,
        sample_interval=sample_interval,
        start_frame=start_frame,
        n_frames=n_frames,
    )


if __name__ == "__main__":
    main()
