# Molecular Crystal Conversion Tools

## Overview

This project includes specialized tools for converting molecular crystal structures from CIF format to simulation-ready formats for GROMACS and CHARMM. These tools ensure **consistent atom ordering** across all molecule copies in the crystal, which is critical for proper topology application.

## Tool Suite

### 1. [mol_utils.py](mol_utils.md) - Core Utilities

Molecular graph and coordinate manipulation functions:
- **`build_graph()`**: Create covalent bond graphs from atomic structures
- **`split_molecules_pbc()`**: Detect individual molecules under periodic boundary conditions
- **`unwrap_to_single_image()`**: Unwrap molecules across periodic boundaries
- Graph isomorphism matching with RMSD-based disambiguation

Used internally by all conversion tools.

---

### 2. [cif2gro.py](cif2gro.md) - CIF → GROMACS GRO Converter

Convert CIF files to GROMACS GRO format with consistent atom ordering.

**Location**: `src/engines/openmm/cif2gro.py`

**Features**:
- Graph-based molecule detection
- Deterministic atom ordering via graph isomorphism
- Optional supercell generation (e.g., 3×3×3)
- Template PDB export for topology generation
- Optional recentering/wrapping

**Quick Start**:
```bash
# Basic conversion
python src/engines/openmm/cif2gro.py --cif input.cif --out output.gro

# 3×3×3 supercell with template
python src/engines/openmm/cif2gro.py \
    --cif pentacene.cif \
    --super 3 3 3 \
    --out pentacene_333.gro \
    --template-pdb pentacene_template.pdb
```

---

### 3. [build_top.py](build_top.md) - GROMACS Topology Generator

Create GROMACS topology files (.top) from GRO structures.

**Location**: `src/engines/openmm/build_top.py`

**Features**:
- Automatic molecule counting
- Forcefield and ITP file inclusion
- System name customization

**Quick Start**:
```bash
# Basic topology
python src/engines/openmm/build_top.py pentacene.gro

# Custom parameters
python src/engines/openmm/build_top.py rubrene.gro \
    --itp rubrene.itp \
    --atoms-per-mol 70 \
    --mol-name RUB \
    -o topol.top
```

---

### 4. [cif2crd.py](cif2crd.md) - CIF → CHARMM CRD/PDB/mol2 Converter

Convert CIF files to CHARMM-compatible formats.

**Location**: `src/engines/pycharmm/cif2crd.py`

**Outputs**:
- **CRD**: CHARMM coordinate file (extended format)
- **PDB**: Visualization format
- **mol2**: SYBYL format for CGenFF parameterization

**Features**:
- All outputs have consistent atom ordering
- SYBYL atom type inference (C.3, C.ar, N.2, etc.)
- Bond information for CGenFF
- Supercell generation

**Quick Start**:
```bash
# All outputs (CRD, PDB, mol2)
python src/engines/pycharmm/cif2crd.py \
    --cif pentacene.cif \
    --supercell 3 3 3

# Custom filenames
python src/engines/pycharmm/cif2crd.py \
    --cif input.cif \
    --supercell 2 2 2 \
    --crd system.crd \
    --pdb system.pdb \
    --mol2 template.mol2
```

---

### 5. [build_psf.py](build_psf.md) - CHARMM PSF Topology Generator

Generate CHARMM PSF topology files using pyCHARMM.

**Location**: `src/engines/pycharmm/build_psf.py`

**Inputs**:
- **PDB**: Single molecule template (defines atom names)
- **CRD**: Full system coordinates (extended format)
- **Stream file**: CGenFF topology/parameters (`.str`)

**Outputs**:
- **PSF**: CHARMM topology (bonds, angles, dihedrals)
- **PDB**: Optional output with proper connectivity

**Features**:
- No external CHARMM executable required
- Automatic connectivity generation
- Crystal structure support
- Extended format for large systems

**Quick Start**:
```bash
# Basic PSF generation
python src/engines/pycharmm/build_psf.py \
    --pdb pentacene.pdb \
    --crd pentacene_333.crd \
    --toppar pentacene.str \
    --nmol 27 \
    --psf system.psf
```

---

## Workflow Examples

### GROMACS Workflow

```bash
# 1. Convert CIF to GRO with template
python src/engines/openmm/cif2gro.py \
    --cif pentacene.cif \
    --super 3 3 3 \
    --out pentacene.gro \
    --template-pdb pentacene_template.pdb

# 2. Generate parameters with LigParGen (external)
# Visit https://traken.chem.yale.edu/ligpargen/
# Upload pentacene_template.pdb → download pentacene.itp

# 3. Create topology
python src/engines/openmm/build_top.py pentacene.gro \
    --itp pentacene.itp \
    --atoms-per-mol 36

# 4. Run GROMACS simulation
gmx grompp -f md.mdp -c pentacene.gro -p topol.top -o run.tpr
gmx mdrun -deffnm run
```

### CHARMM Workflow

```bash
# 1. Convert CIF to CHARMM formats
python src/engines/pycharmm/cif2crd.py \
    --cif pentacene.cif \
    --supercell 3 3 3 \
    --crd pentacene.crd \
    --mol2 pentacene.mol2

# 2. Generate parameters with CGenFF (external)
# Visit https://cgenff.umaryland.edu/
# Upload pentacene.mol2 → download parameters

# 3. Create PSF (using build_psf.py or CHARMM directly)
# ... (to be documented)

# 4. Run CHARMM or pyCHARMM simulation
```

---

## Common Options

All conversion tools support:

| Option | Description |
|--------|-------------|
| `--bond-scale` | Adjust covalent radius cutoff (default: 1.10) |
| `--no-recenter` | Disable automatic recentering |
| `--no-wrap` | Disable wrapping into simulation box |

---

## Algorithm Overview

All tools use the same core algorithm:

1. **Molecule Detection**: Build covalent bond graph under PBC → identify connected components
2. **Template Generation**: Select first molecule, create deterministic atom ordering via BFS
3. **Isomorphism Matching**: Map all molecules to template using graph isomorphism
4. **RMSD Disambiguation**: For symmetric molecules, choose isomorphism with minimum RMSD
5. **Coordinate Processing**: Unwrap → optional recenter → optional wrap
6. **Format-Specific Output**: Write in target format (GRO/CRD/PDB/mol2/TOP)

This ensures **identical atom ordering** for all molecule copies, which is essential for:
- Correct topology application
- Energy conservation during MD
- Meaningful trajectory analysis

---

## Dependencies

- **ASE** (Atomic Simulation Environment): Structure I/O, neighbor lists
- **NetworkX**: Graph algorithms, isomorphism matching
- **NumPy**: Coordinate transformations
- **SciPy**: Kabsch alignment for RMSD minimization

Install with:
```bash
pip install ase networkx numpy scipy
```

---

## See Also

- [mol_utils.md](mol_utils.md) - Detailed API documentation
- [cif2gro.md](cif2gro.md) - GROMACS conversion details
- [cif2crd.md](cif2crd.md) - CHARMM conversion details
- [build_top.md](build_top.md) - GROMACS topology generation
- [architecture.md](architecture.md) - Overall project structure
