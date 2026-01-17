"""
Example NVE simulation using pure MM (force field only)

This script demonstrates how to run NVE production using
only classical force fields (OpenMM) without QM/MM.
"""

from ash import *
import os
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config
from q4mdkit.qmmm.xml import remove_montecarlo_params

# Load configurations
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure from NPT
input_pdb = "../npt/output/npt_lastframe.pdb"
state_xml = "../npt/OpenMM_MD_final_state.xml"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Remove MonteCarlo barostat params from state file for NVE
nve_state_xml = "nve_initial_state.xml"
remove_montecarlo_params(state_xml, nve_state_xml)

print("="*60)
print("  MM NVE Production  ")
print("="*60)

# Load fragment
frag = Fragment(pdbfile=input_pdb)

# Create OpenMM theory (no QM region)
omm = mm_config.create_openmm_theory()

# Print system info
mm_config.print_system_info(frag)

print("\nStarting NVE production...")
md_config.run_nve(frag, omm, output_dir=output_dir, statefile=nve_state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nve")

print("\nDone!")
