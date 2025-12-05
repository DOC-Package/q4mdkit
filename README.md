# qm4dcrystal

## Introduction
Wrapper for QM/MM molecular dynamics simulations of molecular crystals.

## Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Add bin directory to PATH (optional, for convenience)
export PATH="/home/takahashi/python/mdcrystal/bin:$PATH"

# Or source the tools.sh file to use aliases
source tools.sh
```

## Usage

After adding `bin/` to your PATH or sourcing `tools.sh`:

```bash
# topology files workflow
cif2gro --cif pentacene.cif --supercell 3 3 3 --out pentacene.gro
buildtop --gro pentacene.gro --itp pentacene.itp --nmol 27 --out topol.top
```

## Simulation engines

This project supports both OpenMM and pyCHARMM.

- **OpenMM**: a fast, GPU-enabled MD engine

Select the engine via the `SIMULATION_ENGINE` variable in `src/config.py`.
The common interface allows using the same workflow with either engine.

## Molecular Crystal Conversion Tools

This project includes tools for converting molecular crystal structures (CIF format) to simulation-ready formats:

### Quick Start

```bash
# Convert CIF to GROMACS GRO (3x3x3 supercell)
python src/engines/openmm/cif2gro.py \
    --cif input.cif \
    --super 3 3 3 \
    --out output.gro \
    --template-pdb template.pdb

# Generate GROMACS topology
python src/engines/openmm/build_top.py output.gro --itp molecule.itp

# Convert CIF to CHARMM formats
python src/engines/pycharmm/cif2crd.py \
    --cif input.cif \
    --supercell 3 3 3
```

### Key Features

- **Consistent atom ordering**: All molecule copies have identical internal atom order
- **Graph-based detection**: Uses covalent bonding under periodic boundaries
- **RMSD disambiguation**: Handles symmetric molecules correctly
- **Multiple formats**: Supports GROMACS (GRO/TOP) and CHARMM (CRD/PDB/mol2/PSF)

See [docs/tools.md](docs/tools.md) for comprehensive documentation.

---

## Development
```bash
# Run tests
python -m pytest tests/
```
