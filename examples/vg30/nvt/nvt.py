import os
from ash import Fragment
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config
from q4mdkit.qmmm.xml import remove_montecarlo_xml
from qmmm_constraints import capture_qm_hbond_constraints, restore_constraints

qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

input_pdb = "../npt-mm/output/npt_lastframe.pdb"
state_xml = "../npt-mm/OpenMM_MD_final_state.xml"
state_nve = "nve_initial_state.xml"
remove_montecarlo_xml(state_xml, state_nve)

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

print("="*60)
print("  NVT Equilibration (300 K)  ")
print("="*60)

frag = Fragment(pdbfile=input_pdb)
qmatoms = qmmm_config.load_qmatoms()

omm = qmmm_config.create_openmm_theory(fragment=frag)
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

print("\nStarting NVT equilibration...")
md_config.run_nvt(frag, qmmm, output_dir=output_dir, statefile=state_nve)

print("\nSaving final structure...")
md_config.save_final_structure(output_dir=output_dir, prefix="nvt")
