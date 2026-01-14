#!/usr/bin/env python
"""
Select QM atoms from pentacene 3x3x3 supercell PDB file
"""
import sys
from pathlib import Path
from qm4d4crystal.prep.select_qmatoms import select_qmatoms

input_dir = Path(__file__).parent
result = select_qmatoms(
    structure_file=str(input_dir / "pentacene.pdb"),
    n_molecules=2,           # QM region: 2 molecules
    n_active_molecules=10,   # Active region: 10 molecules (for optimization)
    output=str(input_dir / "qmatoms"),
    active_output=str(input_dir / "active_atoms"),
    verbose=True,
    distance_threshold=0.3,  # Show candidates within 0.3 nm of best option
    # choice=1,              # Uncomment to auto-select (1-based index)
)
    
print(f"\nReturned {len(result['qm_atoms'])} QM atoms")
print(f"Returned {len(result['active_atoms'])} active atoms")