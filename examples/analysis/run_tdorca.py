#!/usr/bin/env python3
"""
Run ORCA TD-DFT analysis along a trajectory.

Usage:
    python run_tdorca.py [config.yaml]
"""

import sys
from pathlib import Path

# Add parent directory to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from q4mdkit.analysis.tdorca import run_tdorca_analysis


def main():
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = Path(__file__).parent / "tdorca_settings.yaml"
    
    run_tdorca_analysis(str(config_path))


if __name__ == "__main__":
    main()
