#!/usr/bin/env python3
"""
Constrained DFT-B analysis script.

Usage:
    python run_cdftb.py [config.yaml]

If no config file is specified, defaults to 'cdftb_settings.yaml' in the current directory.
"""

import argparse
from pathlib import Path
from qm4d4crystal.analysis.cdftb import run_cdftb_analysis

def main():
    
    parser = argparse.ArgumentParser(
        description="Run constrained DFT-B analysis on MD trajectory"
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
        sys.exit(1)
        
    print(f"Running CDFTB calculations with config: {config_path}")
    run_cdftb_analysis(config_path)


if __name__ == "__main__":
    main()
