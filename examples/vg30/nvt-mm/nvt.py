import argparse
import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
sys.path.insert(0, str(SCRIPT_DIR.parent))

from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config
from validate_inputs import validate_vg30_system
from workflow_paths import MINIMIZED_PDB

parser = argparse.ArgumentParser(description="Run VG30 MM NVT equilibration.")
parser.add_argument(
    "--validate-only",
    action="store_true",
    help="validate the prepared system without importing ASH or running MD",
)
args = parser.parse_args()

# Load configurations
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# input: initial structure
input_pdb = str(MINIMIZED_PDB)
mm_config.pdbfile = input_pdb

# Refuse a minimized structure produced with an obsolete or mismatched box.
validate_vg30_system(input_pdb, mm_config.amber_prmtop, mm_config.boxfile)
if args.validate_only:
    print("VG30 NVT inputs are consistent.")
    raise SystemExit(0)

from ash import Fragment

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

print("\nStarting NVT equilibration...")
md_config.run_nvt(frag, omm, output_dir=output_dir)

print("\nDone!")
