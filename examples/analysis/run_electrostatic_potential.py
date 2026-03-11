#!/usr/bin/env python
"""
Run electrostatic potential analysis.

Computes displacement-induced electrostatic potential fluctuations
at QM atoms due to MM point charges.

Usage:
    python run_electrostatic_potential.py                         # Use default config file
    python run_electrostatic_potential.py config.yaml             # Use specified config file
"""

import sys
from pathlib import Path
from q4mdkit.analysis.electrostatic_potential import run_potential_analysis

def main():
    # Default config file
    default_config = Path(__file__).parent / "potential_settings.yaml"
    
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = str(default_config)
    
    print(f"Running electrostatic potential analysis with config: {config_path}")
    results = run_potential_analysis(config_path)
    
    # Print summary statistics
    import numpy as np
    phi_disp = results['phi_disp']
    print(f"\nSummary:")
    print(f"  Frames analyzed: {phi_disp.shape[0]}")
    print(f"  QM atoms: {phi_disp.shape[1]}")
    print(f"  φ_disp mean: {np.mean(phi_disp):.4f} V")
    print(f"  φ_disp std:  {np.std(phi_disp):.4f} V")
    print(f"  φ_disp range: [{np.min(phi_disp):.4f}, {np.max(phi_disp):.4f}] V")


if __name__ == "__main__":
    main()
