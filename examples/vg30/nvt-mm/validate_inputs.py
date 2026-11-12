#!/usr/bin/env python3
"""Compatibility CLI for VG30 inputs; validation is implemented in q4mdkit."""

import argparse

from q4mdkit.validation import load_box_vectors, read_prmtop_atom_count, validate_amber_system

EXPECTED_ATOMS = 77 + 920 * 6
EXPECTED_VG30_ATOMS = 77


def validate_vg30_system(pdb_path, prmtop_path, box_path):
    return validate_amber_system(
        pdb_path, prmtop_path, box_path, expected_atoms=EXPECTED_ATOMS,
        residue_blocks=[{'name': 'VG3', 'count': 77}, {'name': 'ACE', 'count': 5520}],
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
