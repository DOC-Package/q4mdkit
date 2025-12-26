#!/usr/bin/env python
"""
Run standard DFTB trajectory analysis.

Usage:
    python run_dftb.py                    # Use default config file
    python run_dftb.py config.yaml        # Use specified config file
"""

import sys
from pathlib import Path
from qm4d4crystal.analysis.dftb import run_dftb_analysis

def main():
    # Default config file
    default_config = Path(__file__).parent / "dftb_settings.yaml"
    
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = str(default_config)
    
    print(f"Running DFTB analysis with config: {config_path}")
    run_dftb_analysis(config_path)


if __name__ == "__main__":
    main()
