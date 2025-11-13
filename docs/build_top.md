# build_top.py - GROMACS Topology File Generator

## Overview

Creates a GROMACS topology file (`.top`) from a GRO structure file for molecular crystal systems.

## Usage

```bash
python build_top.py [gro_file] [options]
```

## Arguments

### Positional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `gro_file` | `pentacene.gro` | Input GRO structure file |

### Optional Arguments

| Option | Default | Description |
|--------|---------|-------------|
| `--itp` | `pentacene.itp` | Molecule topology ITP file containing atom types, bonds, angles, etc. |
| `--ff`, `--forcefield` | `oplsaa.ff/forcefield.itp` | Forcefield ITP file with parameter definitions |
| `--atoms-per-mol` | `36` | Number of atoms per molecule (used to calculate number of molecules) |
| `--mol-name` | `MOL` | Molecule name as defined in the ITP file's `[ moleculetype ]` section |
| `--system-name` | `pentacene_crystal` | System name for the `[ system ]` section |
| `-o`, `--output` | `topol.top` | Output topology file path |

## Examples

### Basic Usage (Default Parameters)

```bash
python build_top.py pentacene.gro
```

Creates `topol.top` for a pentacene crystal with 36 atoms per molecule.

### Custom Molecule

```bash
python build_top.py rubrene.gro --itp rubrene.itp --atoms-per-mol 70 --mol-name RUB
```

### Specify All Parameters

```bash
python build_top.py system.gro \
    --itp molecule.itp \
    --ff amber99sb.ff/forcefield.itp \
    --atoms-per-mol 50 \
    --mol-name MOLNAME \
    --system-name "My Crystal System" \
    -o topology.top
```

## Output

The script creates a GROMACS topology file with the following structure:

```
; Topology file for <system_name>
; Generated from <gro_file>

; Include forcefield
#include "<forcefield>"

; Include molecule topology
#include "<itp_file>"

[ system ]
<system_name>

[ molecules ]
; Compound        #mols
<mol_name>        <n_molecules>
```

## Validation

- **Atom Count Check**: Warns if total atoms in GRO file is not evenly divisible by `atoms_per_mol`
- **File Requirements**: Requires existing GRO file (ITP and forcefield paths are not validated)

## Python API

```python
from build_top import build_topology_file

info = build_topology_file(
    gro_file="pentacene.gro",
    itp_file="pentacene.itp",
    system_name="pentacene_crystal",
    atoms_per_molecule=36,
    forcefield="oplsaa.ff/forcefield.itp",
    molecule_name="MOL",
    output_file="topol.top"
)

# Returns:
# {
#     'n_atoms': 1944,
#     'n_molecules': 54,
#     'output_file': 'topol.top',
#     'gro_file': 'pentacene.gro',
#     'itp_file': 'pentacene.itp'
# }
```

## Notes

- The GRO file must contain coordinates for all atoms in the crystal
- The ITP file must define the molecule topology with the name specified in `--mol-name`
- The forcefield ITP file must be accessible from the working directory or use absolute paths
- Number of molecules is calculated as: `n_molecules = n_atoms / atoms_per_mol`

## Related Files

- **build_psf_pdb.py**: Creates CHARMM PSF topology files
- **cif2crd.py**: Converts CIF crystal structures to CHARMM CRD format
- **make_gro.py**: Converts structures to GROMACS GRO format
