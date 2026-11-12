"""Shared setup for the VG30 ORCA QM/MM constant-volume examples."""

import argparse
import math
import os
from pathlib import Path
import sys

from q4mdkit.qmmm.md_config import get_md_config
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.xml import read_restart_box, remove_montecarlo_xml


VG30_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(VG30_DIR / "nvt-mm"))
sys.path.insert(0, str(VG30_DIR / "nvt"))

from qmmm_constraints import capture_qm_hbond_constraints, restore_constraints
from validate_inputs import EXPECTED_ATOMS, validate_vg30_system


def run(ensemble, stage_dir):
    """Run a stage in its own directory so ASH restart and QM files stay local."""
    parser = argparse.ArgumentParser(description=f"Run VG30 ORCA QM/MM {ensemble.upper()}.")
    parser.add_argument("--validate-only", action="store_true",
                        help="check restart inputs without importing ASH or running ORCA")
    parser.add_argument("--simulation-time", type=float, metavar="PS",
                        help="override the simulation duration in ps")
    args = parser.parse_args()
    os.chdir(stage_dir)

    qmmm_config = get_config("qmmm_settings.yaml")
    md_config = get_md_config("md_settings.yaml")
    predecessor = "npt-mm" if ensemble == "nvt" else "nvt-orca"
    statefile = f"../{predecessor}/OpenMM_MD_final_state.xml"
    input_pdb = qmmm_config.pdbfile
    validate_vg30_system(input_pdb, qmmm_config.amber_prmtop, None)
    qmatoms = qmmm_config.load_qmatoms()
    if qmatoms != list(range(77)):
        raise ValueError("VG30 QM atom list must contain exactly indices 0 through 76")
    if qmmm_config.qm_backend != "orca":
        raise ValueError("These examples require qm_backend: orca")
    box_vectors = read_restart_box(statefile, expected_atoms=EXPECTED_ATOMS)
    if args.simulation_time is not None:
        if not math.isfinite(args.simulation_time) or args.simulation_time < md_config.timestep:
            parser.error("--simulation-time must be finite and at least one timestep")
    if args.validate_only:
        print(f"VG30 ORCA QM/MM {ensemble.upper()} inputs are consistent.")
        print(f"Initial state: {statefile}")
        return

    output_dir = md_config.output_dir
    existing = [path for path in (
        Path(output_dir) / f"{ensemble}.dcd", Path("OpenMM_MD_final_state.xml"),
    ) if path.exists()]
    if existing:
        parser.error("Move existing results before starting a new run: "
                     + ", ".join(str(path) for path in existing))

    from ash import Fragment

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    initial_state = "initial_state.xml"
    remove_montecarlo_xml(statefile, initial_state)
    frag = Fragment(pdbfile=input_pdb)
    # State XML keeps more box precision than the CRYST1 record in the PDB.
    omm = qmmm_config.create_openmm_theory(fragment=frag, pbc_vectors=box_vectors)
    constraints = capture_qm_hbond_constraints(omm.system, qmatoms, frag.elems)
    qm_theory = qmmm_config.create_qm_theory()
    qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm=omm, qm_theory=qm_theory)
    restore_constraints(omm.system, constraints)

    qmmm_config.print_config()
    md_config.print_config()
    runner = md_config.run_nvt if ensemble == "nvt" else md_config.run_nve
    runner(frag, qmmm, output_dir=output_dir, statefile=initial_state,
           simulation_time=args.simulation_time)
    md_config.save_final_structure(output_dir=output_dir, prefix=ensemble)
