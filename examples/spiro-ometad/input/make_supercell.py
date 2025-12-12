#!/usr/bin/env python
"""
Generate 3x3x3 supercell PDB and mol2 template from spiro-ometad CIF for GAFF2/AMBER

Available functions:
- cif_to_pdb:          Slow but exhaustive (O(n!) graph isomorphism)
- cif_to_pdb_wlhash:   Fast AND accurate (O(n log n) WL-hash matching) ★推奨
- cif_to_pdb_coordsort: Fastest but less accurate (coordinate sorting only)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from prep.cif2pdb import cif_to_pdb_wlhash  # WLハッシュ版（高速・正確）

input_dir = Path(__file__).parent
result = cif_to_pdb_wlhash(
    cif_file=str(input_dir / "spiro-ometad.cif"),
    output_pdb=str(input_dir / "spiro-ometad.pdb"),
    output_mol2=str(input_dir / "spiro-ometad_template.mol2"),
    output_box=str(input_dir / "spiro-ometad.box"),
    supercell=(3, 3, 3),
    resname="SPR",
    mol_name="SPR"
)