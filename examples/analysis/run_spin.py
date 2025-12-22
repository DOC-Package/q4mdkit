#!/usr/bin/env python3
"""
Compute spin populations for fragments from CDFTB calculations.

Usage:
    python run_spin_populations.py [config.yaml]

If no config file is specified, defaults to 'cdftb_settings.yaml' in the current directory.
"""

import argparse
from pathlib import Path

from qm4d4crystal.analysis.spin import compute_spin_populations_from_cdftb


def main():
    parser = argparse.ArgumentParser(
        description="Compute spin populations for fragments from CDFTB calculations"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftb_settings.yaml",
        help="Path to YAML configuration file (default: cdftb_settings.yaml)"
    )
    args = parser.parse_args()
    
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        exit(1)
    
    print(f"Computing spin populations with config: {config_path}")
    compute_spin_populations_from_cdftb(config_path)


if __name__ == "__main__":
    main()
