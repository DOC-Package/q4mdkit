#!/usr/bin/env python3
"""Validate that the VG30 NVT inputs describe the same periodic system."""

import argparse
import math
from pathlib import Path


EXPECTED_ATOMS = 77 + 920 * 6
EXPECTED_VG30_ATOMS = 77


def load_box_vectors(box_path):
    vectors = []
    for line in Path(box_path).read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            values = [float(value) for value in stripped.split()]
            if len(values) != 3:
                raise ValueError(f"invalid box-vector row in {box_path}: {line}")
            vectors.append(values)
    if len(vectors) != 3:
        raise ValueError(f"expected three box vectors in {box_path}")
    return vectors


def read_prmtop_atom_count(prmtop_path):
    lines = Path(prmtop_path).read_text().splitlines()
    try:
        pointer_flag = next(
            index for index, line in enumerate(lines) if line.strip() == "%FLAG POINTERS"
        )
        return int(lines[pointer_flag + 2].split()[0])
    except (StopIteration, IndexError, ValueError) as error:
        raise ValueError(f"cannot read atom count from {prmtop_path}") from error


def validate_vg30_system(pdb_path, prmtop_path, box_path):
    # If no explicit box file is provided, skip box-vector comparison.
    if box_path is None:
        box_vectors = None
        box_lengths = None
    else:
        box_vectors = load_box_vectors(box_path)
        box_lengths = [math.sqrt(sum(value * value for value in vector)) for vector in box_vectors]

    pdb_lines = Path(pdb_path).read_text().splitlines()
    atom_lines = [
        line for line in pdb_lines if line.startswith(("ATOM  ", "HETATM"))
    ]
    if len(atom_lines) != EXPECTED_ATOMS:
        raise ValueError(
            f"PDB has {len(atom_lines)} atoms; expected {EXPECTED_ATOMS} for VG30 + 920 ACE"
        )

    topology_atoms = read_prmtop_atom_count(prmtop_path)
    if topology_atoms != len(atom_lines):
        raise ValueError(
            f"topology has {topology_atoms} atoms but PDB has {len(atom_lines)}"
        )

    solute_resnames = {line[17:20].strip() for line in atom_lines[:EXPECTED_VG30_ATOMS]}
    solvent_resnames = {line[17:20].strip() for line in atom_lines[EXPECTED_VG30_ATOMS:]}
    if solute_resnames != {"VG3"} or solvent_resnames != {"ACE"}:
        raise ValueError("PDB residue ordering is not 77 VG3 atoms followed by ACE")

    cryst1 = next((line for line in pdb_lines if line.startswith("CRYST1")), None)
    if cryst1 is None:
        raise ValueError(f"PDB has no CRYST1 record: {pdb_path}")
    # If a box file was provided, compare its lengths to the PDB CRYST1 values.
    if box_lengths is not None:
        pdb_lengths = [float(cryst1[start:end]) for start, end in ((6, 15), (15, 24), (24, 33))]
        if any(abs(pdb_value - box_value) > 1e-3 for pdb_value, box_value in zip(pdb_lengths, box_lengths)):
            raise ValueError(
                f"PDB CRYST1 box {pdb_lengths} does not match box file {box_lengths}; "
                "rerun opt-mm with the current box before NVT"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--box", required=True)
    args = parser.parse_args()
    try:
        validate_vg30_system(args.pdb, args.prmtop, args.box)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print("VG30 NVT inputs are consistent.")


if __name__ == "__main__":
    main()
