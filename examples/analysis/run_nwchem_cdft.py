#!/usr/bin/env python
"""
Run NWChem constrained-DFT (CDFT) analysis along a trajectory.

Generates an NWChem CDFT input for each selected trajectory frame, executes
NWChem, and writes total DFT energies and constrained populations.

Usage:
    python run_nwchem_cdft.py                            # default config
    python run_nwchem_cdft.py nwchem_cdft_settings.yaml  # specified config
"""

import sys
from pathlib import Path

from q4mdkit.analysis.nwchem_cdft import run_nwchem_cdft_analysis


def main() -> None:
    default_config = Path(__file__).parent / "nwchem_cdft_settings.yaml"
    config_path = sys.argv[1] if len(sys.argv) > 1 else str(default_config)

    print(f"Running NWChem CDFT analysis with config: {config_path}")
    run_nwchem_cdft_analysis(config_path)


if __name__ == "__main__":
    main()
