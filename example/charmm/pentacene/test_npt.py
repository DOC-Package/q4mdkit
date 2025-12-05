#!/usr/bin/env python3
"""Test NPT trajectory output."""
import sys
sys.path.insert(0, '/home/takahashi/python/mdcrystal/src')
from engines.pycharmm import run_md

run_md.load_base_topology()
run_md.load_psf('pentacene.psf', 'pentacene.crd', 'pentacene.str')
run_md.setup_pbc_and_pme()
run_md.setup_shake()
print("==> Running NPT...")
run_md.run_npt(nsteps=200, dt_ps=0.002, traj_path='test_npt.dcd', nsavc=50, use_shake=True)
print("==> NPT done")
