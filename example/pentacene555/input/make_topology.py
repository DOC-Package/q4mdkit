#!/usr/bin/env python
"""
Generate GROMACS topology file for pentacene 4x4x4 supercell
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from prep.build_top import build_topology_file

input_dir = Path(__file__).parent
build_topology_file(
    gro_file=str(input_dir / "pentacene.gro"),
    itp_file="pentacene.itp",
    system_name="Pentacene Crystal",
    atoms_per_molecule=36,  # C22H14 = 36 atoms
    molecule_name="PEN",
    output_file=str(input_dir / "pentacene_v1.top")
)
