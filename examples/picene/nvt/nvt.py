import os
import sys
from ash import Fragment
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config
import mdtraj as md

# Load configurations from parent directory
qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "../input/picene.pdb"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NVT Equilibration (300 K)  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)
qmatoms = qmmm_config.load_qmatoms()

omm = qmmm_config.create_openmm_theory()
qm_dftb = qmmm_config.create_dftb_theory()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_dftb)

print("\nStarting NVT equilibration...")
md_config.run_nvt(frag, qmmm, output_dir=output_dir)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nvt")

