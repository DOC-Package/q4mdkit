#!/usr/bin/env python
"""
Select QM atoms from pentacene 4x4x4 supercell GRO file
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from prep.select_qmatoms import select_qmatoms

input_dir = Path(__file__).parent
result = select_qmatoms(
    grofile=str(input_dir / "pentacene.gro"),
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