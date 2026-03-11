#!/usr/bin/env python3
"""
Select molecules near the QM region from ph-btbt-c10 supercell PDB file.
This script finds neighboring molecules around the QM molecule.
"""
from pathlib import Path
from q4mdkit.prep.select_nearby_molecules import (
    select_nearby_molecules,
    read_atom_indices,
    generate_vmd_script
)

input_dir = Path(__file__).parent

# Settings
cutoff = 1.5              # Distance cutoff in nm
n_shells = None           # Number of coordination shells (alternative to cutoff)
use_pbc = True            # Apply periodic boundary conditions
verbose = True            # Show detailed information
generate_vmd = True       # Generate VMD visualization script

result = select_nearby_molecules(
    structure_file=str(input_dir / "ph-btbt-c10.pdb"),
    qm_atoms_file=str(input_dir / "qmatoms"),
    cutoff=cutoff,
    n_shells=n_shells,
    output=str(input_dir / "nearby_atoms"),
    mol_output=str(input_dir / "nearby_molecules"),
    use_pbc=use_pbc,
    verbose=verbose
)

print(f"\nReturned {len(result['nearby_atoms'])} nearby atoms")
print(f"Returned {len(result['nearby_molecules'])} nearby molecules")
print(f"QM molecules: {result['qm_molecules']}")
print(f"Nearby molecules: {result['nearby_molecules']}")

# Generate VMD script
if generate_vmd:
    qm_atoms = read_atom_indices(str(input_dir / "qmatoms"))
    generate_vmd_script(
        structure_file=str(input_dir / "ph-btbt-c10.pdb"),
        qm_atoms=qm_atoms,
        nearby_atoms=result['nearby_atoms'],
        output_script=str(input_dir / "visualize_nearby.tcl"),
        show_all_system=False
    )

