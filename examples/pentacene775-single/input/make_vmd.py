#!/usr/bin/env python
"""
Generate VMD visualization script for pentacene 3x3x3 QM/MM system
"""
import sys
from pathlib import Path
from qm4d4crystal.prep.visualize_qmmm import make_vmd_qmmm

input_dir = Path(__file__).parent

result = make_vmd_qmmm(
    structure_file=str(input_dir / "pentacene.pdb"),
    qmatoms_file=str(input_dir / "qmatoms"),
    output_pdb=str(input_dir / "pentacene_qmmm.pdb"),
    vmd_script=str(input_dir / "visualize_vmd.tcl")
)
    
print(f"\nOutput files:")
print(f"  PDB: {result['pdb']}")
print(f"  VMD script: {result['vmd_script']}")
