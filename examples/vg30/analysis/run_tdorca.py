#!/usr/bin/env python3
"""Run ORCA TDDFT on trajectory snapshots, following run_tddftb.py."""

import argparse
from pathlib import Path
import sys

# Allow execution from any directory in a source checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "config", nargs="?", type=Path,
        default=Path(__file__).with_name("tdorca_settings.yaml"),
        help="YAML settings (default: tdorca_settings.yaml beside this script)",
    )
    args = parser.parse_args()

    from q4mdkit.analysis.tdorca import run_tdorca_analysis

    if not run_tdorca_analysis(str(args.config)):
        parser.exit(1, "No successful ORCA TDDFT snapshots; check the inputs and work/orca.out.\n")


if __name__ == "__main__":
    main()
