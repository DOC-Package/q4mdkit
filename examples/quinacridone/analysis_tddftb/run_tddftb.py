#!/usr/bin/env python3
"""
Run TD-DFTB trajectory analysis for quinacridone in DMSO.

Computes excitation energies and oscillator strengths using TD-DFTB
(Linear Response / Casida) for each frame of the NVE trajectory.

Usage:
    python run_tddftb.py

Output:
    - output_tddftb/excitations.dat: Raw excitation data for each frame
"""

import sys
from pathlib import Path

# Add parent directories to path for q4mdkit imports
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))

from q4mdkit.analysis.tddftb import run_tddftb_analysis


def main():
    config_file = Path(__file__).parent / "tddftb_settings.yaml"
    
    print("=" * 60)
    print("  Quinacridone TD-DFTB Trajectory Analysis")
    print("=" * 60)
    print(f"Config: {config_file}")
    print()
    
    run_tddftb_analysis(str(config_file))
    
    print()
    print("Analysis complete!")


if __name__ == "__main__":
    main()
