from ash import *
import os
import sys
sys.path.insert(0, os.path.abspath('../../../src'))
from qmmm.qmmm_config import get_config
from qmmm.md_config import get_md_config
import mdtraj as md

# Load configurations from parent directory
qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "../nvt/output/nvt_lastframe.pdb"
state_xml = "../nvt/OpenMM_MD_final_state.xml"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NPT Equilibration (300 K, 1 bar)  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)
qmatoms = qmmm_config.load_qmatoms()

omm = qmmm_config.create_openmm_theory()
qm_dftb = qmmm_config.create_dftb_theory()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_dftb)

print("\nStarting NPT equilibration...")
md_config.run_npt(frag, qmmm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="npt")

