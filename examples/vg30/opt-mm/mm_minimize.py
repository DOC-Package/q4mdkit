#!/usr/bin/env python3
"""
MM energy minimization before MD.
"""

import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
sys.path.insert(0, str(SCRIPT_DIR.parent))

from ash import Fragment, OpenMM_Opt, OpenMMTheory

sys.path.insert(0, str(SCRIPT_DIR.parent / "nvt-mm"))
from validate_inputs import load_box_vectors
from workflow_paths import AMBER_PRMTOP, BOX_FILE, MINIMIZED_PDB, PACKED_PDB

# Input/output
input_pdb = str(PACKED_PDB)

print("="*60)
print("  MM Energy Minimization  ")
print("="*60)

# Create fragment
frag = Fragment(pdbfile=input_pdb)
print(f"Loaded {frag.numatoms} atoms")

# Create OpenMM theory WITHOUT constraints for minimization
PBCvectors = load_box_vectors(BOX_FILE)
omm = OpenMMTheory(
    Amberfiles=True,
    amberprmtopfile=str(AMBER_PRMTOP),
    fragment=frag,
    PBCvectors=PBCvectors,
    periodic=True,
    autoconstraints="HBonds", 
    rigidwater=False,
    periodic_nonbonded_cutoff=9.0,
    numcores=4
)

# Energy minimization
print("\nRunning energy minimization (tolerance=0.1)...")
OpenMM_Opt(fragment=frag, theory=omm, maxiter=5000, tolerance=0.1)

# Save minimized structure
frag.write_pdbfile(str(MINIMIZED_PDB))
#frag.write_xyzfile(f"{output_dir}/minimized.xyz")
print(f"\nMinimized structure saved to {MINIMIZED_PDB}")
