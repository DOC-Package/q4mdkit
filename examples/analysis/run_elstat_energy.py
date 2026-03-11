#!/usr/bin/env python3
"""
Run electrostatic energy analysis.

Calculate the electrostatic interaction energy between QM atoms and MM point charges.

Usage:
    python run_elstat_energy.py config.yaml
    
Or modify the configuration below and run directly.
"""

import sys
from pathlib import Path

# Add parent to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from q4mdkit.analysis.electrostatic_energy import (
    analyze_electrostatic_energy,
    run_from_config_file,
)

# =============================================================================
# Configuration
# =============================================================================
# Input files
DCD_FILE = '../nve/output/nve.dcd'
TOP_FILE = '../input/system.pdb'
QM_INDEX_FILE = '../input/qmatoms'
QM_CHARGES_FILE = 'qm_charges.dat'       # Charges of QM atoms
MM_INDEX_FILE = '../input/mmatoms'       # Indices of MM atoms (specific molecule)

# MM charges: either direct file or PCcharges.dat (auto-extract)
MM_CHARGES_FILE = None                   # Direct charges file (if provided)
PCCHARGES_FILE = '../input/PCcharges.dat'  # ASH/DFTB output (auto-extract)

# Optional: Reference structure (required if FIX_QM=True)
REF_XYZ_FILE = '../input/qm_ref.xyz'

# Fix QM molecule? If True, QM molecule undergoes rigid body motion only
FIX_QM = False

# Output
OUTPUT_PREFIX = 'output_elstat_energy'

# Frame selection
START_FRAME = 0
END_FRAME = None  # None = all frames
STRIDE = 1

# Options
CHARGE_COLUMN = -1  # -1 = last column
VERBOSE = True


# =============================================================================
# Main
# =============================================================================
def main():
    # Check for config file argument
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
        if Path(config_file).exists() and config_file.endswith(('.yaml', '.yml')):
            print(f"Running from config: {config_file}")
            run_from_config_file(config_file)
            return
    
    # Run with hardcoded configuration
    analyze_electrostatic_energy(
        dcd_file=DCD_FILE,
        top_file=TOP_FILE,
        qm_index_file=QM_INDEX_FILE,
        qm_charges_file=QM_CHARGES_FILE,
        mm_index_file=MM_INDEX_FILE,
        mm_charges_file=MM_CHARGES_FILE,
        pccharges_file=PCCHARGES_FILE,
        output_prefix=OUTPUT_PREFIX,
        ref_xyz_file=REF_XYZ_FILE if FIX_QM else None,
        fix_qm=FIX_QM,
        verbose=VERBOSE,
        start_frame=START_FRAME,
        end_frame=END_FRAME,
        stride=STRIDE,
        charge_column=CHARGE_COLUMN,
    )


if __name__ == "__main__":
    main()
