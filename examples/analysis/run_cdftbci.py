#!/usr/bin/env python3
"""
CDFTB-CI analysis script.

Compute the CDFTB-CI Hamiltonian and transfer integrals from CDFTB output.

Usage:
    python run_cdftbci.py [config.yaml]

If no config file is specified, defaults to 'cdftb_settings.yaml' in the current directory.
"""

import argparse
from pathlib import Path

from qm4d4crystal.analysis.cdftbci import run_cdftbci_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Run CDFTB-CI analysis on MD trajectory"
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
    
    print(f"Running CDFTB-CI analysis with config: {config_path}")
    run_cdftbci_analysis(config_path)


if __name__ == "__main__":
    main()
