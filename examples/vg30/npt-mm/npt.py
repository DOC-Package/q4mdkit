import argparse
import os
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
sys.path.insert(0, str(SCRIPT_DIR.parent / "nvt-mm"))

from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config
from sync_box import sync_box
from validate_inputs import validate_vg30_system


parser = argparse.ArgumentParser(description="Run VG30 MM NPT equilibration.")
parser.add_argument(
    "--validate-only",
    action="store_true",
    help="validate the prepared system without importing ASH or running MD",
)
args = parser.parse_args()

mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")
input_pdb = "../nvt-mm/output/nvt_lastframe.pdb"
state_xml = "../nvt-mm/OpenMM_MD_final_state.xml"
mm_config.pdbfile = input_pdb

validate_vg30_system(input_pdb, mm_config.amber_prmtop, mm_config.boxfile)
if not Path(state_xml).is_file():
    raise FileNotFoundError(f"NVT final state not found: {state_xml}")

if args.validate_only:
    print("VG30 MM NPT inputs are consistent.")
    raise SystemExit(0)

from ash import Fragment


output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("=" * 60)
print("  VG30 MM NPT Equilibration (300 K, 1 bar)  ")
print("=" * 60)

frag = Fragment(pdbfile=input_pdb)
omm = mm_config.create_openmm_theory(fragment=frag)

print("\nStarting NPT equilibration...")
md_config.run_npt(frag, omm, output_dir=output_dir, statefile=state_xml)

print("\nSaving final structure and precise box...")
md_config.save_final_structure(output_dir=output_dir, prefix="npt")
sync_box(
    "OpenMM_MD_final_state.xml",
    os.path.join(output_dir, "npt_lastframe.box"),
)
