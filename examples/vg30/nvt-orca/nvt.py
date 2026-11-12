#!/usr/bin/env python3
"""ORCA QM/MM NVT, following ../nvt/nvt.py."""

import argparse
import math
import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
sys.path.insert(0, str(SCRIPT_DIR.parent))

from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config
from q4mdkit.qmmm.xml import remove_montecarlo_xml
from orca_md import read_restart_box
from qmmm_constraints import capture_qm_hbond_constraints, restore_constraints
from validate_inputs import validate_vg30_system


parser = argparse.ArgumentParser(description="Run VG30 ORCA QM/MM NVT.")
parser.add_argument("--validate-only", action="store_true",
                    help="check restart inputs without importing ASH or running ORCA")
parser.add_argument("--simulation-time", type=float, metavar="PS",
                    help="override the simulation duration in ps")
args = parser.parse_args()

qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

input_pdb = qmmm_config.pdbfile
state_xml = "../npt-mm/OpenMM_MD_final_state.xml"
state_nvt = "initial_state.xml"
output_dir = md_config.output_dir

validate_vg30_system(input_pdb, qmmm_config.amber_prmtop, None)
qmatoms = qmmm_config.load_qmatoms()
if qmatoms != list(range(77)):
    raise ValueError("VG30 QM atom list must contain exactly indices 0 through 76")
if qmmm_config.qm_backend != "orca":
    raise ValueError("This example requires qm_backend: orca")
box_vectors = read_restart_box(state_xml)
if args.simulation_time is not None:
    if not math.isfinite(args.simulation_time) or args.simulation_time < md_config.timestep:
        parser.error("--simulation-time must be finite and at least one timestep")
if args.validate_only:
    print("VG30 ORCA QM/MM NVT inputs are consistent.")
    print(f"Initial state: {state_xml}")
    raise SystemExit(0)

existing = [path for path in (
    Path(output_dir) / "nvt.dcd", Path("OpenMM_MD_final_state.xml"),
) if path.exists()]
if existing:
    parser.error("Move existing results before starting a new run: "
                 + ", ".join(str(path) for path in existing))

from ash import Fragment

remove_montecarlo_xml(state_xml, state_nvt)
os.makedirs(output_dir, exist_ok=True)

print("=" * 60)
print(f"  ORCA QM/MM NVT Equilibration ({md_config.temperature:g} K)  ")
print("=" * 60)

frag = Fragment(pdbfile=input_pdb)

# Use the restart State's full precision cell rather than the rounded PDB cell.
omm = qmmm_config.create_openmm_theory(fragment=frag, pbc_vectors=box_vectors)
qm_hbond_constraints = capture_qm_hbond_constraints(
    omm.system,
    qmatoms,
    frag.elems,
)
qm_theory = qmmm_config.create_qm_theory()
qmmm = qmmm_config.create_qmmm_theory(
    frag,
    qmatoms,
    omm=omm,
    qm_theory=qm_theory,
)
restore_constraints(omm.system, qm_hbond_constraints)

qmmm_config.print_config()
md_config.print_config()

print("\nStarting NVT equilibration...")
md_config.run_nvt(
    frag,
    qmmm,
    output_dir=output_dir,
    statefile=state_nvt,
    simulation_time=args.simulation_time,
)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nvt")
