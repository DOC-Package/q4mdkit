"""
Example NPT simulation using pure MM (force field only)

This script demonstrates how to run NPT equilibration using
only classical force fields (OpenMM) without QM/MM.
"""

from ash import *
import os
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config

# Load configurations
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure from NVT
input_pdb = "../nvt/output/nvt_lastframe.pdb"
state_xml = "../nvt/OpenMM_MD_final_state.xml"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  MM NPT Equilibration (300 K, 1 bar)  ")
print("="*60)

# Load fragment
frag = Fragment(pdbfile=input_pdb)

# Create OpenMM theory (no QM region)
omm = mm_config.create_openmm_theory()

# Print system info
mm_config.print_system_info(frag)

print("\nStarting NPT equilibration...")
md_config.run_npt(frag, omm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="npt")

print("\nDone!")
