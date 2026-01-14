"""
Select the central molecule from pentacene supercell PDB file for QM region.
This script selects only one molecule - the one closest to the center.
"""
import sys
from pathlib import Path
from qm4d4crystal.prep.select_qmatoms import select_central_molecule

input_dir = Path(__file__).parent
result = select_central_molecule(
    structure_file=str(input_dir / "pentacene.pdb"),
    n_active_molecules=10,   # Active region: 10 molecules (for optimization)
    output=str(input_dir / "qmatoms"),
    active_output=str(input_dir / "active_atoms"),
    verbose=True,
)
    
print(f"\nReturned {len(result['qm_atoms'])} QM atoms")
print(f"Returned {len(result['active_atoms'])} active atoms")
print(f"Central molecule ID: {result['central_molecule']['id']}")
