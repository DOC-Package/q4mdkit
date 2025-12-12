#!/usr/bin/env python
"""
Generate VMD visualization script for pentacene 3x3x3 QM/MM system
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from prep.visualize_qmmm import make_vmd_qmmm

input_dir = Path(__file__).parent

result = make_vmd_qmmm(
    structure_file=str(input_dir / "spiro-ometad.pdb"),
    qmatoms_file=str(input_dir / "qmatoms"),
    output_pdb=str(input_dir / "spiro-ometad_qmmm.pdb"),
    vmd_script=str(input_dir / "visualize_vmd.tcl")
)
