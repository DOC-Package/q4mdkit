#!/usr/bin/env python
"""
Run ORCA TD-DFT analysis on Quinacridone NVE trajectory.

This script computes absorption spectra from TD-DFT calculations
for each frame in the trajectory and averages them.

Usage:
    conda activate ash
    python run_tdorca.py
"""

import sys
from pathlib import Path

# Add q4mdkit to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from q4mdkit.analysis.tdorca import run_tdorca_analysis

if __name__ == "__main__":
    config_file = "tdorca_settings.yaml"
    print(f"Running ORCA TD-DFT analysis with config: {config_file}")
    print("="*60)
    
    run_tdorca_analysis(config_file)
    
    print("\nAnalysis complete!")
    print("Output files:")
    print("  - tdorca_output/absorption_spectrum.dat")
    print("  - tdorca_output/excitations.dat")
