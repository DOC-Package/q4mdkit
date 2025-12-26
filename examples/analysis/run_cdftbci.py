#!/usr/bin/env python3
"""
CDFTB-CI analysis script.

Perform CDFTB calculations with online CDFTB-CI calculation for each frame.

This script performs:
1. CDFTB calculations for each fragment (constrained DFT-B)
2. CDFTB-CI calculation immediately after (on-the-fly)
3. Spin population analysis

Features:
- Uses a single working directory (reused for each frame)
- No H/S matrix storage needed  
- Significantly reduced disk usage
- Optional charge continuation between frames

Usage:
    python run_cdftbci.py [config.yaml]

Default config file: cdftb_settings.yaml
"""

import argparse
from pathlib import Path

from qm4d4crystal.analysis.cdftbci import run_cdftbci_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Run CDFTB with online CDFTB-CI calculation"
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
