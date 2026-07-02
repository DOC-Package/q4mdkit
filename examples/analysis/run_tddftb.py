#!/usr/bin/env python
"""
Run TD-DFTB trajectory analysis for excitation energies and oscillator strengths.

Usage:
    python run_tddftb.py                    # Use default config file
    python run_tddftb.py config.yaml        # Use specified config file
"""

import sys
from pathlib import Path
from q4mdkit.analysis.tddftb import run_tddftb_analysis

def main():
    # Default config file
    default_config = Path(__file__).parent / "tddftb_settings.yaml"
    
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = str(default_config)
    
    print(f"Running TD-DFTB analysis with config: {config_path}")
    run_tddftb_analysis(config_path)


if __name__ == "__main__":
    main()
