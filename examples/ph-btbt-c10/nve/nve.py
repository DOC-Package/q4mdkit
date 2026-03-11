from ash import *
import os
import sys
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config
from q4mdkit.qmmm.xml import remove_montecarlo_xml
import mdtraj as md

# Load configurations from parent directory
qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "input/npt_lastframe.pdb"
state_npt = "input/OpenMM_MD_final_state.xml"
# Prepare state XML for NVE (remove MonteCarlo barostat parameters)
state_xml = "input/nve_initial_state.xml"
remove_montecarlo_xml(state_npt, state_xml)

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NVE Calculation  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)
qmatoms = qmmm_config.load_qmatoms()

omm = qmmm_config.create_openmm_theory()
qm_dftb = qmmm_config.create_dftb_theory()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_dftb)

print("\nStarting NVE calculation...")
md_config.run_nve(frag, qmmm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nve")
