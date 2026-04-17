#!/usr/bin/env python
"""
ORCA QM/MM MD Simulation Example

This script demonstrates how to run QM/MM molecular dynamics using
ORCA as the QM backend with ASH and OpenMM.

Prerequisites:
- ASH with OpenMM installed
- ORCA 5.x installed
- AMBER topology files (prmtop, inpcrd)
- QM atom indices file

Usage:
    conda activate ash
    python run_qmmm_md.py

Configuration:
    Edit qmmm_settings.yaml to set:
    - qm_backend: "orca"
    - orca.orcadir: path to ORCA bin directory
    - orca.orcasimpleinput: ORCA method line
    - paths: topology files and atom indices
"""

from ash import *
import os
import sys
from pathlib import Path

# Add q4mdkit to path if not installed
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import MDConfig

def main():
    # Load configuration
    qmmm_config = get_config("qmmm_settings.yaml")
    md_config = MDConfig("md_settings.yaml")
    
    # Print configuration
    qmmm_config.print_config()
    md_config.print_config()
    
    # Create output directory
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load fragment from topology
    frag = Fragment(
        amber_prmtopfile=qmmm_config.amber_prmtop,
        amber_inpcrdfile=qmmm_config.amber_inpcrd
    )
    print(f"\nLoaded structure: {frag.numatoms} atoms")
    
    # Load QM atom indices
    qmatoms = qmmm_config.load_qmatoms()
    print(f"QM region: {len(qmatoms)} atoms")
    
    # Create theories using the configured backend (ORCA in this case)
    omm = qmmm_config.create_openmm_theory()
    qm_theory = qmmm_config.create_qm_theory()  # Uses qm_backend setting
    qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_theory)
    
    # Print system info
    qmmm_config.print_system_info(frag, qmatoms)
    
    # Run NVT equilibration
    print("\n" + "="*60)
    print("Running NVT equilibration...")
    print("="*60)
    md_config.run_nvt(frag, qmmm, output_dir=output_dir)
    
    # Save equilibrated structure
    frag.write_xyzfile(f"{output_dir}/nvt_final.xyz")
    print(f"\nNVT equilibration complete. Structure saved to {output_dir}/nvt_final.xyz")


if __name__ == "__main__":
    main()
