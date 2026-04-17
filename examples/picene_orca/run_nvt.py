#!/usr/bin/env python
"""
Picene ORCA QM/MM NVT Simulation

Run with:
    conda activate ash
    python run_nvt.py
"""

from ash import *
import os
import sys
from pathlib import Path
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import MDConfig

# Load configuration
qmmm_config = get_config("qmmm_settings.yaml")
md_config = MDConfig("md_settings.yaml")

# Print configuration
qmmm_config.print_config()
md_config.print_config()

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Load fragment from AMBER topology
frag = Fragment(
    amber_prmtopfile=qmmm_config.amber_prmtop,
    amber_inpcrdfile=qmmm_config.amber_inpcrd
)
print(f"\nLoaded structure: {frag.numatoms} atoms")

# Load QM atom indices
qmatoms = qmmm_config.load_qmatoms()
print(f"QM region: {len(qmatoms)} atoms (Picene C22H14)")

# Create theories
omm = qmmm_config.create_openmm_theory()
qm_theory = qmmm_config.create_qm_theory()  # Uses ORCA based on qm_backend
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_theory)

# Print system info
qmmm_config.print_system_info(frag, qmatoms)

# Run NVT equilibration
print("\n" + "="*60)
print("Running ORCA QM/MM NVT simulation...")
print("="*60)
md_config.run_nvt(frag, qmmm, output_dir=output_dir)

# Save final structure
frag.write_xyzfile(f"{output_dir}/nvt_final.xyz")
print(f"\nNVT complete. Structure saved to {output_dir}/nvt_final.xyz")
