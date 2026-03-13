#!/usr/bin/env python
"""
MM energy minimization before MD.
"""

import os
from ash import Fragment, OpenMM_Opt, OpenMMTheory
from q4mdkit.qmmm.qmmm_config import get_config

# Load config
qmmm_config = get_config("qmmm_settings.yaml")

# Input/output
input_pdb = "../input/quinacridone_dmso.pdb"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  MM Energy Minimization  ")
print("="*60)

# Create fragment
frag = Fragment(pdbfile=input_pdb)
print(f"Loaded {frag.numatoms} atoms")

# Create OpenMM theory WITHOUT constraints for minimization
PBCvectors = [[34.8, 0.0, 0.0], [0.0, 34.8, 0.0], [0.0, 0.0, 34.8]]
omm = OpenMMTheory(
    Amberfiles=True,
    amberprmtopfile="../input/quinacridone_dmso.prmtop",
    PBCvectors=PBCvectors,
    periodic=True,
    autoconstraints=None, 
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
