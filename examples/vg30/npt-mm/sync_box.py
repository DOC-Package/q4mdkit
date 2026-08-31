#!/usr/bin/env python3
"""Write QM/MM box vectors from an OpenMM final-state XML file."""

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET


def sync_box(state_path, box_path):
    root = ET.parse(state_path).getroot()
    periodic_box = root.find("PeriodicBoxVectors")
    if periodic_box is None:
        raise ValueError(f"PeriodicBoxVectors not found in {state_path}")

    vectors = []
    for name in ("A", "B", "C"):
        vector = periodic_box.find(name)
        if vector is None:
            raise ValueError(f"box vector {name} not found in {state_path}")
        vectors.append([float(vector.attrib[axis]) * 10.0 for axis in ("x", "y", "z")])

    output = Path(box_path)
    with output.open("w") as handle:
        handle.write("# PBC box vectors from OpenMM final-state XML (Angstrom)\n")
        for vector in vectors:
            handle.write(" ".join(f"{value:.10f}" for value in vector) + "\n")
    return vectors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state_xml")
    parser.add_argument("box")
    args = parser.parse_args()
    sync_box(args.state_xml, args.box)


if __name__ == "__main__":
    main()
