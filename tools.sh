#!/bin/bash
# Convenience wrapper scripts for mdcrystal tools
# Source this file to add commands to your shell: source tools.sh

# Get the directory where this script is located
MDCRYSTAL_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PYTHONPATH="${MDCRYSTAL_ROOT}/src:${PYTHONPATH}"

# OpenMM tools
alias cif2gro="python ${MDCRYSTAL_ROOT}/src/engines/openmm/cif2gro.py"
alias buildtop="python ${MDCRYSTAL_ROOT}/src/engines/openmm/build_top.py"

# CHARMM tools
alias cif2crd="python ${MDCRYSTAL_ROOT}/src/engines/pycharmm/cif2crd.py"
alias buildpsf="python ${MDCRYSTAL_ROOT}/src/engines/pycharmm/build_psf.py"

echo "mdcrystal tools loaded!"
echo "Available commands:"
echo "  cif2gro   - Convert CIF to GROMACS format"
echo "  buildtop  - Build GROMACS topology"
echo "  cif2crd   - Convert CIF to CHARMM format"
echo "  buildpsf  - Build CHARMM PSF topology"
echo ""
echo "Example usage:"
echo "  cif2gro --cif pentacene.cif --supercell 3 3 3"
