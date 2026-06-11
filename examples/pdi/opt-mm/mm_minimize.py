#!/usr/bin/env python
"""
MM energy minimization before MD.
"""

import os
from ash import Fragment, OpenMM_Opt, OpenMMTheory

# Input/output
input_pdb = "../input/pdi_acetonitrile.pdb"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  MM Energy Minimization  ")
print("="*60)

# Create fragment
frag = Fragment(pdbfile=input_pdb)
print(f"Loaded {frag.numatoms} atoms")

# ASH/OpenMM atom indices are 0-based. These correspond to PDI N1 (PDB atom 5)
# and C21 (PDB atom 27), defining one imide C-N bond to keep constrained.
i_n1 = 4
i_c1 = 26
i_c2 = 28
i_n2 = 5
i_c3 = 27
i_c4 = 29

# Create OpenMM theory WITHOUT constraints for minimization
PBCvectors = [[34.8, 0.0, 0.0], [0.0, 34.8, 0.0], [0.0, 0.0, 34.8]]
omm = OpenMMTheory(
    Amberfiles=True,
    amberprmtopfile="../input/pdi_acetonitrile.prmtop",
    fragment=frag,
    PBCvectors=PBCvectors,
    periodic=True,
    autoconstraints="HBonds", 
    bondconstraints=[[i_n1, i_c1], [i_n1, i_c2], [i_n2, i_c3], [i_n2, i_c4]],
    rigidwater=False,
    periodic_nonbonded_cutoff=9.0,
    numcores=4
)

# Energy minimization
print("\nRunning energy minimization (tolerance=0.1)...")
OpenMM_Opt(fragment=frag, theory=omm, maxiter=5000, tolerance=0.1)

# Save minimized structure
frag.write_pdbfile(f"{output_dir}/minimized.pdb")
#frag.write_xyzfile(f"{output_dir}/minimized.xyz")
print(f"\nMinimized structure saved to {output_dir}/minimized.pdb")
