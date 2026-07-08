# Post-processing and QM analysis module
#
# Submodules in this package have heavy optional dependencies
# (MDAnalysis, mdtraj, hsd, ...). Import each block defensively so that
# users who only need a subset of the functionality (e.g. the pure-NumPy
# absorption-spectrum utilities) do not have to install everything.

# Always-available, dependency-light modules ---------------------------------
from .absorption import (
    absorption_spectrum,
    dipole_correlation,
    lineshape_function,
    oscillator_strength_to_mu2,
)

# Optional / heavy-dependency modules ---------------------------------------
try:
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
except ImportError:
    pass

try:
    from .dcd2xyz import (
        convert_dcd_to_xyz,
        write_xyz,
    )
except ImportError:
    pass

try:
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
except ImportError:
    pass

try:
    from .remove_transrot import (
        remove_translation_rotation,
    )
except ImportError:
    pass

try:
    from .average_structure import (
        compute_average_structure,
        compute_average_from_traj,
    )
except ImportError:
    pass

try:
    from .electrostatic_energy import (
        ElectrostaticEnergyConfig,
        analyze_electrostatic_energy,
        compute_electrostatic_energy,
        run_electrostatic_energy_analysis,
        run_from_config_file as run_electrostatic_energy_from_config,
    )
except ImportError:
    pass
