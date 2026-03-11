#!/usr/bin/env python3
"""
Process trajectory: remove translation/rotation and compute average structure.

Usage:
    python process_trajectory.py
"""

import sys
from q4mdkit.analysis import remove_translation_rotation, compute_average_structure

# Input files
dcd_file = "mol.dcd"
top_file = "mol.pdb"

# Output files
aligned_dcd = "mol_aligned.dcd"
average_pdb = "mol_average.pdb"

# Step 1: Remove translation and rotation
print("=" * 60)
print("Step 1: Remove translation and rotation")
print("=" * 60)
traj = remove_translation_rotation(dcd_file, top_file, aligned_dcd)

# Step 2: Compute average structure from aligned trajectory
print("\n" + "=" * 60)
print("Step 2: Compute average structure")
print("=" * 60)
avg_coords = compute_average_structure(aligned_dcd, top_file, average_pdb)

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Aligned trajectory: {aligned_dcd}")
print(f"Average structure: {average_pdb}")
