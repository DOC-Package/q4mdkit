from ash import *
import os
import sys
from ash import Fragment
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config
from q4mdkit.qmmm.xml import remove_montecarlo_xml
import mdtraj as md

# Load configurations from parent directory
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "../npt-mm/output/npt_lastframe.pdb"
state_npt = "../npt-mm/OpenMM_MD_final_state.xml"
# Prepare state XML for NVE (remove MonteCarlo barostat parameters)
state_xml = "nve_initial_state.xml"
remove_montecarlo_xml(state_npt, state_xml)

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NVE Calculation  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)

omm = mm_config.create_openmm_theory()

print("\nStarting NVE calculation...")
md_config.run_nve(frag, omm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nve")
