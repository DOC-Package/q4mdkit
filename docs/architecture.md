# Architecture

## Overview
This document describes the project architecture and directory layout used in this example Python project.

## Directory structure explanation

### `src/` - Source code
Contains the main source code for the project.

- **`main.py`**: Entry point for the application
- **`config.py`**: Centralized configuration (constants, paths, simulation parameters)
- **`common/`**: Shared modules used by different engines (crystal generation, analysis, utilities)

### `tests/` - Tests
Unit and integration tests live here.

- Test files should be prefixed with `test_`
- Create test files corresponding to modules under `src/`

### `data/` - Data files
Input and output data such as initial structures and simulation results.

- **`input/`**: Input files
- **`output/`**: Generated outputs (simulation results, plots)

### `docs/` - Documentation
Architecture documents, usage guides and API references.

## Module design

### Core Modules

#### 1. `config.py` - Configuration
- Central place for settings used across the project
- Define constants, paths and parameters
- Allow environment-specific overrides

#### 2. `common/crystal_structure.py` - Crystal structure
- Generate and manage crystal lattices (FCC, BCC, etc.)
- Compute atomic positions and box vectors

#### 3. `common/analysis.py` - Analysis
- Statistical processing of results
- Visualization and file output

#### 4. `common/utils.py` - Utilities
- Logging setup
- Common helper functions

#### 5. `common/mol_utils.py` - Molecular graph utilities
- Graph-based molecule detection under PBC
- Graph isomorphism matching
- Coordinate unwrapping and transformations
- Shared by all CIF conversion tools

### Simulation Engines

#### `engines/base.py` - Base simulator interface
- Abstract base class for all simulation engines
- Common methods: setup(), run(), record_data()
- Engine-specific abstract methods: setup_system(), run_dynamics(), etc.

#### `engines/openmm/` - OpenMM engine
- **`simulator.py`**: OpenMM MD simulation implementation
- **`cif2gro.py`**: CIF → GROMACS GRO converter with consistent atom ordering
- **`build_top.py`**: GROMACS topology file generator

#### `engines/pycharmm/` - pyCHARMM engine
- **`simulator.py`**: pyCHARMM MD simulation implementation
- **`cif2crd.py`**: CIF → CHARMM CRD/PDB/mol2 converter
- **`build_psf.py`**: CHARMM PSF topology generator

### Tool Architecture

The CIF conversion tools (`cif2gro.py`, `cif2crd.py`) share a common architecture:

1. **Input**: CIF file with crystal structure
2. **Graph detection**: Build covalent bond graph → identify molecules
3. **Template ordering**: Create deterministic atom order for first molecule
4. **Isomorphism matching**: Map all molecules to template using graph isomorphism
5. **RMSD disambiguation**: Handle symmetric molecules via Kabsch alignment
6. **Coordinate processing**: Unwrap, recenter, wrap as needed
7. **Output**: Write in target format (GRO, CRD, PDB, mol2, etc.)

This ensures **consistent atom ordering** across all molecule copies, which is critical for:
- Proper topology application
- Energy conservation in MD simulations
- Trajectory analysis

## Design principles

### 1. Separation of concerns
Each module has a clear responsibility and can operate independently.

### 2. Centralized configuration
Use `config.py` to centralize settings for easy modification.

### 3. Testability
Modules are designed to be independently testable.

### 4. Readability
- Clear naming conventions
- Adequate documentation
- Consistent code style

### 5. Path handling
Use `Path` objects instead of hard-coded absolute paths for portability.

## Execution flow

1. Load configuration (`config.py`)
2. Ensure directories exist
3. Generate crystal structure (`common/crystal_structure.py`)
4. Run simulation (engine-specific simulator)
5. Analyze results (`common/analysis.py`)
6. Save data and plots

## Extending the project

### Adding new features
1. Add a new module under `src/common/` or `src/engines/`
2. Invoke it from `main.py` as needed
3. Add corresponding tests under `tests/`

### Adding configuration
1. Add new settings to `config.py`
2. Optionally load from environment variables or a config file

### Adding data processing
1. Add new analysis functions to `common/analysis.py` or create a new analysis module
