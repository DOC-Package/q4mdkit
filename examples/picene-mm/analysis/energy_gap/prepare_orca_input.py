#!/usr/bin/env python3
"""Convert Antechamber's ORCA scaffold into an open-shell AM1 job."""

from __future__ import annotations

import argparse
from pathlib import Path
import re


METHOD_LINE = "! AM1 ROHF SlowConv TightSCF Mulliken"
SCF_BLOCK = ["%scf", "  MaxIter 500", "end"]
XYZ_HEADER = re.compile(r"^(\s*\*\s+xyz)\s+[-+]?\d+\s+\d+\s*$", re.IGNORECASE)


def convert_input(text: str, charge: int, multiplicity: int) -> str:
    if multiplicity < 1:
        raise ValueError("multiplicity must be at least 1")

    lines = text.splitlines()
    method_indexes = [i for i, line in enumerate(lines) if line.lstrip().startswith("!")]
    if len(method_indexes) != 1:
        raise ValueError(
            f"expected exactly one ORCA method line, found {len(method_indexes)}"
        )

    xyz_indexes = [i for i, line in enumerate(lines) if XYZ_HEADER.match(line)]
    if len(xyz_indexes) != 1:
        raise ValueError(
            f"expected exactly one '* xyz CHARGE MULTIPLICITY' line, found {len(xyz_indexes)}"
        )

    lines[method_indexes[0]] = METHOD_LINE
    match = XYZ_HEADER.match(lines[xyz_indexes[0]])
    assert match is not None
    lines[xyz_indexes[0]] = f"{match.group(1)} {charge} {multiplicity}"
    lines[method_indexes[0] + 1 : method_indexes[0] + 1] = SCF_BLOCK
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--charge", required=True, type=int)
    parser.add_argument("--multiplicity", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError("input and output paths must differ")
    converted = convert_input(
        args.input.read_text(),
        charge=args.charge,
        multiplicity=args.multiplicity,
    )
    args.output.write_text(converted)
    print(f"Prepared ORCA AM1 input: {args.output}")


if __name__ == "__main__":
    main()
