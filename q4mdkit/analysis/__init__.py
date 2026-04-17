# Post-processing and QM analysis module
"""
from .cdftbci import run_cdftbci_analysis
from .dftb import run_dftb_analysis

__all__ = [
    'run_cdftbci_analysis',
    'run_dftb_analysis',
]
"""

from .normal_mode_analysis import (
    NormalModeAnalyzer,
    analyze_trajectory,
    calculate_displacements,
    project_onto_normal_modes,
    load_trajectory_mdtraj,
    read_gen_file,
    read_hessian_eigenvectors,
    read_hessian_eigenvalues,
)

from .dcd2xyz import (
    convert_dcd_to_xyz,
    write_xyz,
)

from .electrostatic_potential import (
    ElectrostaticPotentialAnalyzer,
    ElectrostaticPotentialConfig,
    analyze_electrostatic_potential,
    compute_electrostatic_potential,
    compute_displacement_potential,
    compute_reference_positions,
    load_config as load_electrostatic_potential_config,
    run_potential_analysis,
    run_from_config as run_electrostatic_potential_from_config,
)

from .remove_transrot import (
    remove_translation_rotation,
)

from .average_structure import (
    compute_average_structure,
    compute_average_from_traj,
)

from .electrostatic_energy import (
    ElectrostaticEnergyConfig,
    analyze_electrostatic_energy,
    compute_electrostatic_energy,
    run_electrostatic_energy_analysis,
    run_from_config_file as run_electrostatic_energy_from_config,
)

# TD-DFTB module (requires hsd package)
# from .tddftb import (
#     TDDFTBConfig,
#     Excitation,
#     run_tddftb_analysis,
#     compute_absorption_spectrum,
#     parse_exc_dat,
# )

# ORCA TD-DFT module
# from .tdorca import (
#     TDORCAConfig,
#     run_tdorca_analysis,
# )

# )