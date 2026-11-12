#!/usr/bin/env python3
"""Run VG30 NVT, following the other MM examples."""

import os
from pathlib import Path

from q4mdkit.md import prepare_md_run
from q4mdkit.mm.mm_config import get_config
from q4mdkit.mm.md_config import get_md_config

os.chdir(Path(__file__).resolve().parent)
mm_config = get_config("mm_settings.yaml")
md_config = get_md_config("md_settings.yaml")
statefile = prepare_md_run(
    mm_config, md_config, "nvt",
    argv=None, description="VG30 NVT",
)

from ash import Fragment

frag = Fragment(pdbfile=mm_config.pdbfile)
omm = mm_config.create_openmm_theory(fragment=frag)
md_config.run_nvt(frag, omm, statefile=statefile)
md_config.save_final_structure(prefix="nvt")
