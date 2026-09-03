from ash import *
import os
import sys
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config
import mdtraj as md

# Load configurations from parent directory
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "../nvt-mm/output/nvt_lastframe.pdb"
state_xml = "../nvt-mm/OpenMM_MD_final_state.xml"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NPT Equilibration (300 K, 1 bar)  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)

omm = mm_config.create_openmm_theory()

print("\nStarting NPT equilibration...")
md_config.run_npt(frag, omm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="npt")

