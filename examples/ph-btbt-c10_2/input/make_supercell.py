#!/usr/bin/env python
"""
Generate 3x3x3 supercell PDB and mol2 template from ph-btbt-c10 CIF for GAFF2/AMBER
"""
import sys
from pathlib import Path
from q4mdkit.prep.cif2pdb import cif_to_pdb

input_dir = Path(__file__).parent
result = cif_to_pdb(
    cif_file=str(input_dir / "ph-btbt-c10.cif"),
    output_pdb=str(input_dir / "ph-btbt-c10.pdb"),
    output_mol2=str(input_dir / "ph-btbt-c10_template.mol2"),
    output_box=str(input_dir / "ph-btbt-c10.box"),
    supercell=(6, 6, 2),
    resname="BTB",
    mol_name="BTB"
)
    
print(f"\nCreated {result['n_molecules']} molecules ({result['n_atoms']} atoms)")
print(f"PDB: {result['pdb']}")
print(f"MOL2: {result['mol2']}")
print(f"BOX: {result['box']}")