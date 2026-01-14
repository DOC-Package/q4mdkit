#!/usr/bin/env python3
"""
Step 1: MM-only Geometry Optimization

This script performs classical MM minimization to relax the initial structure
before QM/MM optimization. This helps avoid SCC convergence issues in DFTB
by providing a more stable starting geometry.
"""
from ash import *
import os
import sys
from qm4d4crystal.qmmm.qmmm_config import get_config

# Load QM/MM configuration
qmmm_config = get_config("qmmm_settings.yaml")

# Input/Output
input_pdb = "../input/spiro-ometad.pdb"
output_dir = "output_mm"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  Step 1: MM-only Geometry Optimization  ")
print("="*60)

# Load fragment
frag = Fragment(pdbfile=input_pdb)
print(f"\nLoaded structure: {frag.numatoms} atoms")

# Create MM theory (OpenMM only)
omm = qmmm_config.create_openmm_theory()

print("\nStarting MM minimization...")
print("  Max iterations: 500")

# Run geometry optimization with MM only
Optimizer(
    fragment=frag,
    theory=omm,
    maxiter=500
)

# Save optimized structure
output_pdb = f"{output_dir}/mm_optimized.pdb"
frag.write_pdbfile(output_pdb)
print(f"\nSaved MM-optimized structure: {output_pdb}")

# Also save XYZ for inspection
output_xyz = f"{output_dir}/mm_optimized.xyz"
frag.write_xyzfile(output_xyz)
print(f"Saved XYZ file: {output_xyz}")
