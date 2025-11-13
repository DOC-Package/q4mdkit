# mdcrystal

## Introduction
Wrapper for molecular dynamics simulations of molecular crystals.

## Directory structure
```
.
├── src/              # Main source code
│   ├── main.py       # Entry point
│   ├── config.py     # Configuration
│   ├── common/       # Shared modules (OpenMM / pyCHARMM common code)
│   │   ├── crystal_structure.py
│   │   ├── analysis.py
│   │   └── utils.py
│   └── engines/      # Simulation engine implementations
│       ├── base.py
│       ├── openmm/
│       │   └── simulator.py
│       └── pycharmm/
│           └── simulator.py
├── tests/            # Tests
├── data/             # Data files
│   ├── input/
│   └── output/
├── docs/             # Documentation
└── requirements.txt  # Dependencies
```

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

### Method 1: Using convenience commands (recommended)

After adding `bin/` to your PATH or sourcing `tools.sh`:

```bash
cd example/pentacene

# OpenMM workflow
cif2gro --cif pentacene.cif --supercell 3 3 3 --out pentacene.gro
buildtop --gro pentacene.gro --itp pentacene.itp --nmol 27 --out topol.top

# CHARMM workflow
cif2crd --cif pentacene.cif --supercell 3 3 3
buildpsf --pdb pentacene.pdb --crd pentacene.crd --toppar pentacene.str --nmol 27
```

### Method 2: Using PYTHONPATH

```bash
cd example/pentacene
export PYTHONPATH=/home/takahashi/python/mdcrystal/src:$PYTHONPATH

python ../../src/engines/openmm/cif2gro.py --cif pentacene.cif --supercell 3 3 3
```

### Method 3: Running main.py

```bash
# Run the program (default: OpenMM)
python src/main.py

# To switch engine, edit `src/config.py`:
# SIMULATION_ENGINE = "openmm" or "pycharmm"
```

## Simulation engines

This project supports both OpenMM and pyCHARMM.

- **OpenMM**: a fast, GPU-enabled MD engine
- **pyCHARMM**: traditional CHARMM-based MD engine

Select the engine via the `SIMULATION_ENGINE` variable in `src/config.py`.
The common interface allows using the same workflow with either engine.

## Molecular Crystal Conversion Tools

This project includes tools for converting molecular crystal structures (CIF format) to simulation-ready formats:

### Available Tools

| Tool | Description | Location |
|------|-------------|----------|
| **cif2gro.py** | CIF → GROMACS GRO converter | `src/engines/openmm/` |
| **build_top.py** | GROMACS topology generator | `src/engines/openmm/` |
| **cif2crd.py** | CIF → CHARMM CRD/PDB/mol2 converter | `src/engines/pycharmm/` |
| **build_psf.py** | CHARMM PSF generator | `src/engines/pycharmm/` |
| **mol_utils.py** | Core molecular graph utilities | `src/common/` |

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
