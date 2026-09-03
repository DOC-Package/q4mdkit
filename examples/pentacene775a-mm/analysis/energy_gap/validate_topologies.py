#!/usr/bin/env python3
"""Validate that two Amber topologies differ only in target charges."""

import argparse
from pathlib import Path

import numpy as np

from prepare_charged_topology import (
    load_parmed,
    print_validation_report,
    validate_topologies,
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Validate a neutral/charged Amber topology pair"
    )
    parser.add_argument("--neutral-prmtop", required=True, type=Path)
    parser.add_argument("--charged-prmtop", required=True, type=Path)
    parser.add_argument(
        "--resid", required=True, type=int, help="Target residue number (1-based)"
    )
    parser.add_argument("--expected-delta-charge", type=float)
    parser.add_argument("--charge-tolerance", type=float, default=1.0e-4)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    neutral = load_parmed(args.neutral_prmtop)
    charged = load_parmed(args.charged_prmtop)
    report = validate_topologies(neutral, charged, args.resid)
    print_validation_report(report)
    if args.expected_delta_charge is not None and not np.isclose(
        report.delta_target_charge,
        args.expected_delta_charge,
        atol=args.charge_tolerance,
        rtol=0.0,
    ):
        raise ValueError(
            f"Target charge change {report.delta_target_charge:+.8f} does not match "
            f"expected {args.expected_delta_charge:+.8f} within "
            f"{args.charge_tolerance:g}"
        )


if __name__ == "__main__":
    main()
