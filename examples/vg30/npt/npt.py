import argparse
import os
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
sys.path.insert(0, str(SCRIPT_DIR.parent / "nvt-mm"))

from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config
from validate_inputs import validate_vg30_system


parser = argparse.ArgumentParser(description="Run VG30 OPV-B3LYP QM/MM NPT.")
parser.add_argument(
    "--validate-only",
    action="store_true",
    help="validate the prepared system without importing ASH or running MD",
)
args = parser.parse_args()

qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")
input_pdb = "../npt-mm/output/npt_lastframe.pdb"
state_xml = "../npt-mm/OpenMM_MD_final_state.xml"
qmatoms = qmmm_config.load_qmatoms()

validate_vg30_system(input_pdb, qmmm_config.amber_prmtop, qmmm_config.boxfile)
if qmatoms != list(range(77)):
    raise ValueError("VG30 QM atom list must contain exactly indices 0 through 76")
missing_sk = [path for path in qmmm_config.slater_koster_files.values() if not Path(path).is_file()]
if missing_sk:
    raise FileNotFoundError(f"missing OPV-B3LYP Slater-Koster files: {missing_sk}")
if not Path(state_xml).is_file():
    raise FileNotFoundError(f"NPT-MM final state not found: {state_xml}")

if args.validate_only:
    print("VG30 OPV-B3LYP QM/MM inputs are consistent.")
    print(f"Initial state: {state_xml}")
    raise SystemExit(0)

from ash import Fragment


output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("=" * 60)
print("  VG30 OPV-B3LYP QM/MM NPT (300 K, 1 bar)  ")
print("=" * 60)

frag = Fragment(pdbfile=input_pdb)
omm = qmmm_config.create_openmm_theory(fragment=frag)
qm_theory = qmmm_config.create_qm_theory()
qmmm = qmmm_config.create_qmmm_theory(
    frag,
    qmatoms,
    omm=omm,
    qm_theory=qm_theory,
)

print("\nStarting NPT equilibration...")
md_config.run_npt(
    frag,
    qmmm,
    output_dir=output_dir,
    statefile=state_xml,
)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="npt")
