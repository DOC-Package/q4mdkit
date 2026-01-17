import os
from ash import Fragment
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config

# Load configurations
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = "../input/spiro-ometad.pdb"

# output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  MM NVT Equilibration (300 K)  ")
print("="*60)

# Load fragment
frag = Fragment(pdbfile=input_pdb)

# Create OpenMM theory (no QM region)
omm = mm_config.create_openmm_theory()

# Print system info
mm_config.print_system_info(frag)
md_config.print_config()

# Energy minimization before MD
print("\nRunning energy minimization...")
from ash import OpenMM_Opt
OpenMM_Opt(fragment=frag, theory=omm, maxiter=1000, tolerance=10)

print("\nStarting NVT equilibration...")
md_config.run_nvt(frag, omm, output_dir=output_dir)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nvt")

print("\nDone!")

