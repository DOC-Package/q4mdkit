# qm4dcrystal

QM/MM MD simulations of molecular crystals using ASH + OpenMM + DFTB+.


## Features

- **CIF → Simulation files**: Convert crystal structures to PDB with supercell support
- **AMBER topology**: Generate GAFF2 force field parameters via tleap
- **QM atom selection**: Automatic selection of QM region based on geometric center
- **QM/MM MD**: MD simulations with YAML configuration 
- **QM analysis**: Constrained DFTB calculations

## Requirements

- ASH (with OpenMM, DFTB+)
- ASE, pymatgen
- AMBER tools
- pyyaml
- mdtraj
- Python 3.11+