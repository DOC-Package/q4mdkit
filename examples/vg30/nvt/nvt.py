#!/usr/bin/env python3
"""Run VG30 QM/MM NVT, using the shared q4mdkit setup utilities."""

import os
from pathlib import Path

from q4mdkit.md import prepare_md_run
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config

os.chdir(Path(__file__).resolve().parent)
qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")
statefile = prepare_md_run(
    qmmm_config, md_config, "nvt",
    statefile="../npt-mm/OpenMM_MD_final_state.xml",
    argv=None, description="VG30 QM/MM NVT",
)

from ash import Fragment

frag = Fragment(pdbfile=qmmm_config.pdbfile)
qmatoms = qmmm_config.load_qmatoms()
omm = qmmm_config.create_openmm_theory(fragment=frag)
qm_theory = qmmm_config.create_qm_theory()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm=omm, qm_theory=qm_theory)

md_config.run_nvt(frag, qmmm, statefile=statefile)
md_config.save_final_structure(prefix="nvt")
