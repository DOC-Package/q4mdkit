#!/usr/bin/env python
"""
Generate 4x4x4 supercell GRO file from pentacene CIF
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from prep.cif2gro import cif_to_gro

input_dir = Path(__file__).parent

cif_to_gro(
    cif_file=str(input_dir / "pentacene.cif"),
    output_gro=str(input_dir / "pentacene.gro"),
    supercell=(4, 4, 4),
    resname="PEN"
)
