# qm4dcrystal

## Introduction
Wrapper for QM/MM molecular dynamics simulations of molecular crystals using ASH.

### QM/MM MD
ASH employs OpenMM as a MD engine.

## Configuration Files

### qmmm_settings.yaml
QM/MM calculation settings:
- File paths (SK files, topology, coordinates)
- QM region (charge, multiplicity, atom indices)
- DFTB parameters (Slater-Koster files, Hubbard derivatives)
- OpenMM settings (cutoff, platform, etc.)

### md_settings.yaml
MD simulation settings:
- Timestep, trajectory frequency
- NVT: temperature, coupling frequency, integrator
- NPT: pressure, barostat, barostat frequency, integrator
- NVE: integrator
- Simulation times for each ensemble

## Quick Start

```python
from ash import *
from qmmm.qmmm_config import get_config
from qmmm.md_config import get_md_config

# Load configurations
qmmm_config = get_config("qmmm_settings.yaml")
md_config = get_md_config("md_settings.yaml")

# Setup QM/MM
frag = Fragment(grofile="input.gro")
qmatoms = qmmm_config.load_qmatoms()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms)

# Run MD
md_config.run_nvt(frag, qmmm)
md_config.run_npt(frag, qmmm)
md_config.run_nve(frag, qmmm)
```

## Installation

## Usage

After adding `bin/` to your PATH or sourcing `tools.sh`:

```bash
# topology files workflow
cif2gro --cif pentacene.cif --supercell 3 3 3 --out pentacene.gro
buildtop --gro pentacene.gro --itp pentacene.itp --nmol 27 --out topol.top
```

## Molecular Crystal Conversion Tools

This project includes tools for converting molecular crystal structures (CIF format) to simulation-ready formats: