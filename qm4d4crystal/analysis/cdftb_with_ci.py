"""
CDFTB with Online CDFTB-CI Calculation Module

This module performs constrained DFT-B calculations on MD trajectories
and computes CDFTB-CI transfer integrals on-the-fly for each frame.

Unlike the separate cdftb.py + cdftbci.py workflow, this module:
- Uses a single working directory that is reused for each frame
- Computes CDFTB-CI immediately after CDFTB, avoiding H/S matrix storage
- Significantly reduces disk usage and file count

Output files:
- energies.dat: CDFTB energies for each fragment
- charges.dat: Mulliken charges
- cdftbci.dat: Transfer integrals and adiabatic energies
- cdftbci_sub.dat: Detailed CI quantities (overlaps, weight matrices, etc.)
"""

import numpy as np
import dftbplus
import hsd
import mdtraj as md
from pathlib import Path
import multiprocessing as mp
import os
import shutil
import yaml
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Tuple

from .cdftb import (
    CDFTBConfig,
    FragmentConfig,
    load_config,
    load_qm_indices,
    load_mm_indices,
    load_pccharges,
    iter_qm_coordinates,
    extract_atom_types_from_hsd,
    save_qm_coords_xyz,
    run_dftb_in_subprocess,
    parse_atom_range,
    ANG_PER_NM,
    BOHR_PER_ANG,
)

from .cdftbci import (
    UnrestrictedOrbitalData,
    UnrestrictedCDFTBCIHamiltonian,
    build_cdftbci_hamiltonian_unrestricted,
    solve_cdftbci_unrestricted,
    compute_transfer_integral_unrestricted,
)

from .cdftb_result_reader import (
    load_spin_polarized_calculation_for_ci,
    get_orbital_info_from_eigenvec,
    get_atom_orbital_map,
    build_fragment_weight_matrix,
)

from .spin import compute_fragment_spin_population


@dataclass
class CDFTBWithCIConfig(CDFTBConfig):
    """Extended configuration for CDFTB with online CI calculation."""
    # CI settings
    ci_enabled: bool = True
    ci_N_A: float = 101.0  # Target population for constraint A
    ci_N_B: float = 101.0  # Target population for constraint B
    
    # Output mode
    output_mode: str = "online_ci"  # "online_ci" or "store_frames"
    work_directory: str = "work"
    
    # CI output files
    ci_output_file: Optional[Path] = None
    ci_sub_file: Optional[Path] = None
    
    # Spin output file
    spin_output_file: Optional[Path] = None
    
    # Use previous frame's charges as initial guess
    use_previous_charges: bool = False


def load_config_with_ci(config_path: Path) -> CDFTBWithCIConfig:
    """Load configuration from YAML file with CI settings."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse input paths (relative to config file) - all required
    input_cfg = data.get('input', {})
    if 'trajectory' not in input_cfg:
        raise ValueError("Missing required input: trajectory")
    if 'topology' not in input_cfg:
        raise ValueError("Missing required input: topology")
    if 'qm_atoms' not in input_cfg:
        raise ValueError("Missing required input: qm_atoms")
    if 'pccharges_template' not in input_cfg:
        raise ValueError("Missing required input: pccharges_template")
    if 'hsd_template' not in input_cfg:
        raise ValueError("Missing required input: hsd_template")
    
    traj_path = base_dir / input_cfg['trajectory']
    topology_path = base_dir / input_cfg['topology']
    qm_atoms_file = base_dir / input_cfg['qm_atoms']
    pccharges_template = base_dir / input_cfg['pccharges_template']
    hsd_template = base_dir / input_cfg['hsd_template']
    
    # Parse output paths
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output')
    energy_filename = output_cfg.get('energy_file', 'energies.dat')
    energy_file = output_dir / energy_filename
    charge_filename = output_cfg.get('charge_file', 'charges.dat')
    charge_file = output_dir / charge_filename
    
    # Output mode settings
    output_mode = output_cfg.get('mode', 'online_ci')
    work_directory = output_cfg.get('work_directory', 'work')
    
    # Parse frame settings
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    n_frames = frames_cfg.get('n_frames', None)
    t0_fs = frames_cfg.get('t0_fs', 0.0)
    dt_fs = frames_cfg.get('dt_fs', 4.0)
    
    # Parse DFTB+ settings
    dftb_cfg = data.get('dftb', {})
    dftb_library_path = dftb_cfg.get('library_path', '/home/takahashi/opt/dftb+/lib/libdftbplus.so')
    num_threads = dftb_cfg.get('num_threads', None)
    timeout = dftb_cfg.get('timeout', 300)
    
    # Parse fragment definitions
    fragments = []
    for frag_data in data.get('fragments', []):
        atom_range = frag_data['atom_range']
        if 'charge_sum_range' in frag_data:
            charge_sum_range = frag_data['charge_sum_range']
        else:
            charge_sum_range = parse_atom_range(atom_range)
        # Parse initial_charges path if provided
        initial_charges = None
        if 'initial_charges' in frag_data:
            initial_charges = base_dir / frag_data['initial_charges']
        frag = FragmentConfig(
            name=frag_data['name'],
            atom_range=atom_range,
            charge_sum_range=charge_sum_range,
            initial_charges=initial_charges
        )
        fragments.append(frag)
    
    # Parse CI settings
    ci_cfg = data.get('cdftb_ci', {})
    ci_enabled = ci_cfg.get('enabled', True)
    ci_N_A = ci_cfg.get('N_A', 101.0)
    ci_N_B = ci_cfg.get('N_B', 101.0)
    
    # Parse SCF settings
    scf_cfg = data.get('scf', {})
    use_previous_charges = scf_cfg.get('use_previous_charges', False)
    
    # CI output files
    ci_output_file = output_dir / output_cfg.get('ci_file', 'cdftbci.dat')
    ci_sub_file = output_dir / output_cfg.get('ci_sub_file', 'cdftbci_sub.dat')
    
    # Spin output file
    spin_output_file = output_dir / output_cfg.get('spin_file', 'spin.dat')
    
    return CDFTBWithCIConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        hsd_template=hsd_template,
        output_dir=output_dir,
        energy_file=energy_file,
        charge_file=charge_file,
        start_frame=start_frame,
        n_frames=n_frames,
        t0_fs=t0_fs,
        dt_fs=dt_fs,
        dftb_library_path=dftb_library_path,
        num_threads=num_threads,
        timeout=timeout,
        fragments=fragments,
        ci_enabled=ci_enabled,
        ci_N_A=ci_N_A,
        ci_N_B=ci_N_B,
        output_mode=output_mode,
        work_directory=work_directory,
        ci_output_file=ci_output_file,
        ci_sub_file=ci_sub_file,
        spin_output_file=spin_output_file,
        use_previous_charges=use_previous_charges,
    )


def setup_work_directory(
    work_dir: Path,
    qm_coords_ang: np.ndarray,
    mm_coords_ang: np.ndarray,
    mm_charges: np.ndarray,
    frame_id: int,
    time_fs: float,
    atom_types: List[str],
) -> None:
    """
    Set up or update the working directory with current frame data.
    
    This function reuses the same directory for each frame, only updating
    the coordinate files.
    
    Parameters
    ----------
    work_dir : Path
        Path to working directory.
    qm_coords_ang : np.ndarray
        QM atom coordinates in Angstrom.
    mm_coords_ang : np.ndarray
        MM atom coordinates in Angstrom.
    mm_charges : np.ndarray
        MM point charges.
    frame_id : int
        Current frame number.
    time_fs : float
        Current time in femtoseconds.
    atom_types : List[str]
        Atom type labels for QM atoms.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Save QM coordinates in XYZ format
    save_qm_coords_xyz(work_dir, frame_id, time_fs, qm_coords_ang, atom_types)
    
    # Save PCcharges.dat with updated MM coordinates
    pc_file = work_dir / "PCcharges.dat"
    pc_data = np.column_stack([mm_coords_ang, mm_charges])
    np.savetxt(pc_file, pc_data, fmt="%20.10f %20.10f %20.10f %10.4f")


def setup_fragment_work_directory(
    work_dir: Path,
    fragment_name: str,
    constrained_atoms: str,
    hsd_template: Path,
    read_initial_charges: bool = False,
) -> Path:
    """
    Set up fragment subdirectory within work directory.
    
    Parameters
    ----------
    work_dir : Path
        Parent working directory.
    fragment_name : str
        Name of fragment (e.g., 'fragment1').
    constrained_atoms : str
        Atom range for constraint, e.g., '1:36'.
    hsd_template : Path
        Path to HSD template file.
    read_initial_charges : bool
        If True, set ReadInitialCharges=Yes to use charges.dat from previous frame.
        
    Returns
    -------
    frag_dir : Path
        Path to fragment directory.
    """
    frag_dir = work_dir / fragment_name
    frag_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy and modify HSD template
    with open(hsd_template, 'r') as f:
        data = hsd.load(f)
    
    # Remove Geometry from data (will be added as text later)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Update PCcharges path (use parent directory's file)
    data['Hamiltonian']['DFTB']['ElectricField']['PointCharges']['CoordsAndCharges']['DirectRead']['File'] = '../PCcharges.dat'
    
    # Modify ElectronicConstraints to constrain specific fragment
    if 'ElectronicConstraints' in data['Hamiltonian']['DFTB']:
        constraints = data['Hamiltonian']['DFTB']['ElectronicConstraints']['Constraints']
        if 'MullikenPopulation' in constraints:
            mulliken = constraints['MullikenPopulation']
            if isinstance(mulliken, list):
                for constraint in mulliken:
                    constraint['Atoms'] = constrained_atoms
            else:
                mulliken['Atoms'] = constrained_atoms
    
    # Modify InitialSpins.Atoms to guide spin localization
    if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
        spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
        if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
            atom_spin = spin_pol['Colinear']['InitialSpins']['AtomSpin']
            atom_spin['Atoms'] = constrained_atoms
    
    # Set ReadInitialCharges if using previous frame's charges
    if read_initial_charges:
        data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'Yes'
        # Remove InitialSpins when using ReadInitialCharges (they conflict)
        if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
            spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
            if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
                del spin_pol['Colinear']['InitialSpins']
    
    # Save modified HSD
    hsd_file = frag_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        # Write Geometry section with file include directive first
        f.write('Geometry {\n')
        f.write('  xyzFormat {\n')
        f.write('    <<< "../qm_coords.xyz"\n')
        f.write('  }\n')
        f.write('}\n\n')
        # Then write the rest of the HSD
        hsd.dump(data, f)
    
    return frag_dir


def compute_cdftbci_for_frame(
    work_dir: Path,
    frag1_name: str,
    frag2_name: str,
    E_A: float,
    E_B: float,
    N_A: float,
    N_B: float,
) -> Tuple[Optional[UnrestrictedCDFTBCIHamiltonian], Optional[float], Optional[float], Optional[str]]:
    """
    Compute CDFTB-CI quantities for current frame using data in work directory.
    
    Parameters
    ----------
    work_dir : Path
        Working directory containing fragment subdirectories.
    frag1_name : str
        Name of first fragment directory.
    frag2_name : str
        Name of second fragment directory.
    E_A : float
        Energy of state A (fragment 1 constrained).
    E_B : float
        Energy of state B (fragment 2 constrained).
    N_A : float
        Target population for constraint A.
    N_B : float
        Target population for constraint B.
        
    Returns
    -------
    ham : UnrestrictedCDFTBCIHamiltonian or None
        Hamiltonian object if successful.
    J_direct : float or None
        Direct transfer integral in Hartree.
    J_lowdin : float or None
        Löwdin-orthogonalized transfer integral in Hartree.
    error : str or None
        Error message if failed.
    """
    frag1_dir = work_dir / frag1_name
    frag2_dir = work_dir / frag2_name
    
    try:
        # Load spin-polarized calculation results (without Hamiltonian matrices)
        data_A = load_spin_polarized_calculation_for_ci(frag1_dir)
        data_B = load_spin_polarized_calculation_for_ci(frag2_dir)
        
        # Create UnrestrictedOrbitalData objects
        orb_A = UnrestrictedOrbitalData(
            C_alpha=data_A['C_alpha'],
            C_beta=data_A['C_beta'],
            occ_alpha=data_A['occ_alpha'],
            occ_beta=data_A['occ_beta'],
            n_alpha=data_A['n_alpha'],
            n_beta=data_A['n_beta']
        )
        orb_B = UnrestrictedOrbitalData(
            C_alpha=data_B['C_alpha'],
            C_beta=data_B['C_beta'],
            occ_alpha=data_B['occ_alpha'],
            occ_beta=data_B['occ_beta'],
            n_alpha=data_B['n_alpha'],
            n_beta=data_B['n_beta']
        )
        
        # Read AO overlap matrix
        S_AO = data_A['S']
        n_orbitals = data_A['n_orbitals']
        
        # Build fragment weight matrices
        eigenvec_file = frag1_dir / "eigenvec.out"
        orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
        atom_to_orbitals = get_atom_orbital_map(orbital_info)
        n_atoms = len(atom_to_orbitals)
        n_atoms_half = n_atoms // 2
        
        frag_A_atoms = list(range(1, n_atoms_half + 1))
        frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
        
        w_A = build_fragment_weight_matrix(
            S_AO, frag_A_atoms, n_atoms, atom_to_orbitals=atom_to_orbitals
        )
        w_B = build_fragment_weight_matrix(
            S_AO, frag_B_atoms, n_atoms, atom_to_orbitals=atom_to_orbitals
        )
        
        # Get constraint potentials
        V_A = data_A['Vc']
        V_B = data_B['Vc']
        
        # Build CDFTB-CI Hamiltonian
        ham = build_cdftbci_hamiltonian_unrestricted(
            orb_A, orb_B, S_AO, w_A, w_B,
            E_A, E_B, V_A, V_B, N_A, N_B
        )
        
        # Compute transfer integrals
        J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
        J_lowdin = compute_transfer_integral_unrestricted(ham.H, ham.S, method="lowdin")
        
        return ham, J_direct, J_lowdin, None
        
    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        return None, None, None, error_msg


def compute_spin_for_fragment(
    frag_dir: Path,
    n_atoms_half: int
) -> Tuple[Optional[Dict[str, Dict[str, float]]], Optional[str]]:
    """
    Compute spin expectation values for both fragments from a CDFTB calculation.
    
    Parameters
    ----------
    frag_dir : Path
        Path to fragment directory containing DFTB+ output files.
    n_atoms_half : int
        Half the number of atoms (atoms per fragment).
        
    Returns
    -------
    spin_data : dict or None
        Dictionary with spin populations for 'frag1' and 'frag2':
        Each contains 'n_alpha_frag', 'n_beta_frag', 'spin_density'.
    error : str or None
        Error message if failed.
    """
    try:
        # Load spin-polarized data
        data = load_spin_polarized_calculation_for_ci(frag_dir)
        
        # Get orbital-to-atom mapping
        eigenvec_file = frag_dir / "eigenvec.out"
        orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
        atom_to_orbitals = get_atom_orbital_map(orbital_info)
        n_atoms = len(atom_to_orbitals)
        
        # Define fragment atoms (1-indexed)
        frag_A_atoms = list(range(1, n_atoms_half + 1))
        frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
        
        # Compute spin populations
        spin_frag1 = compute_fragment_spin_population(
            data['C_alpha'], data['C_beta'], data['S'],
            data['n_alpha'], data['n_beta'],
            frag_A_atoms, atom_to_orbitals
        )
        spin_frag2 = compute_fragment_spin_population(
            data['C_alpha'], data['C_beta'], data['S'],
            data['n_alpha'], data['n_beta'],
            frag_B_atoms, atom_to_orbitals
        )
        
        return {'frag1': spin_frag1, 'frag2': spin_frag2}, None
        
    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        return None, error_msg


def run_cdftb_with_ci(config_path: Path) -> None:
    """
    Run CDFTB analysis with online CDFTB-CI calculation.
    
    This is the main entry point for the combined CDFTB + CI workflow.
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file.
    """
    # Load configuration
    config = load_config_with_ci(config_path)
    
    print("=" * 70)
    print("CDFTB with Online CDFTB-CI Analysis")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Trajectory: {config.traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output directory: {config.output_dir}")
    print(f"Output mode: {config.output_mode}")
    print(f"CI enabled: {config.ci_enabled}")
    print("=" * 70)
    
    # Load QM/MM indices
    qm_indices = load_qm_indices(config.qm_atoms_file)
    print(f"Loaded {len(qm_indices)} QM atom indices from {config.qm_atoms_file}")
    
    traj_info = md.load(str(config.traj_path), top=str(config.topology_path), frame=0)
    total_atoms = traj_info.n_atoms
    mm_indices = load_mm_indices(qm_indices, total_atoms)
    print(f"MM atoms: {len(mm_indices)}")
    
    # Load point charges
    mm_charges = load_pccharges(config.pccharges_template)
    print(f"Loaded {len(mm_charges)} point charges from {config.pccharges_template}")
    
    # Extract atom types from HSD template
    atom_types = extract_atom_types_from_hsd(config.hsd_template)
    print(f"Extracted {len(atom_types)} atom types from HSD template")
    
    # Print fragment information
    print(f"\nFragments ({len(config.fragments)}):")
    for frag in config.fragments:
        print(f"  - {frag.name}: atoms {frag.atom_range}, charge sum range {frag.charge_sum_range}")
    
    if config.ci_enabled:
        print(f"\nCDFTB-CI settings:")
        print(f"  N_A = {config.ci_N_A}")
        print(f"  N_B = {config.ci_N_B}")
    
    if config.use_previous_charges:
        print(f"\nSCF settings:")
        print(f"  Use previous charges: Yes")
    print()
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up work directory path
    work_dir = config.output_dir / config.work_directory
    
    # Build headers for output files
    energy_header_parts = ["# Frame", "Time[fs]"]
    for frag in config.fragments:
        energy_header_parts.append(f"Energy_{frag.name}[a.u.]")
    energy_header = "  ".join(energy_header_parts) + "\n"
    
    charge_header_lines = []
    charge_header_lines.append("# Mulliken charges from constrained DFT calculations")
    charge_header_lines.append("# Qn|calcm = charge on fragment n when constraint applied to fragment m")
    charge_header_parts = [f"{'# Frame':>7s}", f"{'Time[fs]':>12s}"]
    for i, frag in enumerate(config.fragments, 1):
        for j, frag2 in enumerate(config.fragments, 1):
            charge_header_parts.append(f"{'Q' + str(j) + '|calc' + str(i):>12s}")
    charge_header_lines.append("  ".join(charge_header_parts))
    charge_header = "\n".join(charge_header_lines) + "\n"
    
    # Open output files
    energy_file = open(config.energy_file, "w")
    energy_file.write(energy_header)
    
    charge_file = open(config.charge_file, "w")
    charge_file.write(charge_header)
    
    ci_file = None
    ci_sub_file = None
    spin_file = None
    
    if config.ci_enabled:
        ci_file = open(config.ci_output_file, "w")
        ci_file.write("# CDFTB-CI Results (Online Calculation)\n")
        ci_file.write("# Frame  Time(fs)   J_lowdin(meV)        E1(Ha)        E2(Ha)        dE(eV)\n")
        
        ci_sub_file = open(config.ci_sub_file, "w")
        ci_sub_file.write("# CDFTB-CI Sub Values\n")
        ci_sub_file.write("# Frame  Time(fs)      E_A(Ha)      E_B(Ha)        H_AB(Ha)        J_direct(meV)        "
                         "V_A(Ha)        V_B(Ha)      N_A      N_B        S_AB        "
                         "S_AB_alpha      S_AB_beta        W_BA        W_BA_alpha      W_BA_beta        "
                         "W_AB        W_AB_alpha      W_AB_beta\n")
    
    # Always open spin output file (useful even without CI)
    spin_file = open(config.spin_output_file, "w")
    spin_file.write("# Spin Populations from CDFTB Calculations\n")
    spin_file.write("# Spin density = N_α - N_β on each fragment\n")
    spin_file.write("# For each constraint state (A or B):\n")
    spin_file.write("#   frag1 = first fragment atoms, frag2 = second fragment atoms\n")
    frag1_name = config.fragments[0].name if len(config.fragments) > 0 else "frag1"
    frag2_name = config.fragments[1].name if len(config.fragments) > 1 else "frag2"
    spin_file.write(f"# Frame  Time(fs)  "
                   f"{frag1_name}_A_nalpha  {frag1_name}_A_nbeta  {frag1_name}_A_spin  "
                   f"{frag2_name}_A_nalpha  {frag2_name}_A_nbeta  {frag2_name}_A_spin  "
                   f"{frag1_name}_B_nalpha  {frag1_name}_B_nbeta  {frag1_name}_B_spin  "
                   f"{frag2_name}_B_nalpha  {frag2_name}_B_nbeta  {frag2_name}_B_spin\n")
    
    try:
        processed_count = 0
        
        for frame_id, time_fs, qm_coords_bohr, mm_coords_ang in iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices, mm_indices
        ):
            # Skip frames before start_frame
            if frame_id < config.start_frame:
                continue
            
            # Stop if we've processed enough frames
            if config.n_frames is not None and processed_count >= config.n_frames:
                break
            
            # Convert QM coords to Angstrom
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Calculate time using user-defined t0 and dt
            time_val = config.t0_fs + frame_id * config.dt_fs
            time_str = f"t = {time_val:.3f} fs"
            
            print(f"Frame {frame_id:05d} ({time_str})")
            
            # Set up work directory with current coordinates
            setup_work_directory(
                work_dir, qm_coords_ang, mm_coords_ang, mm_charges,
                frame_id, time_fs, atom_types
            )
            
            # Process each fragment
            energies = []
            charges = []
            cdftb_success = True
            
            for frag in config.fragments:
                frag_dir = work_dir / frag.name
                charges_file = frag_dir / "charges.dat"
                
                # Determine if we should use initial charges
                use_initial_charges = False
                
                if config.use_previous_charges:
                    if processed_count == 0:
                        # First frame: use initial_charges from config if provided
                        if frag.initial_charges is not None and frag.initial_charges.exists():
                            use_initial_charges = True
                            # Copy initial charges to work directory
                            frag_dir.mkdir(parents=True, exist_ok=True)
                            shutil.copy(frag.initial_charges, charges_file)
                    else:
                        # Subsequent frames: use previous frame's charges if they exist
                        if charges_file.exists():
                            use_initial_charges = True
                
                # Set up fragment directory
                frag_dir = setup_fragment_work_directory(
                    work_dir, frag.name, frag.atom_range, config.hsd_template,
                    read_initial_charges=use_initial_charges
                )
                
                # Run DFTB+ calculation
                energy, mcharge, error = run_dftb_in_subprocess(
                    qm_coords_bohr, frag_dir, write_hs=False,
                    dftb_library_path=config.dftb_library_path,
                    num_threads=config.num_threads,
                    timeout=config.timeout
                )
                
                if error:
                    print(f"  {frag.name}: ERROR - {error}")
                    frag_charge = float("nan")
                    energy = float("nan")
                    all_frag_charges = [float("nan")] * len(config.fragments)
                    cdftb_success = False
                else:
                    start_idx, end_idx = frag.charge_sum_range
                    frag_charge = float(np.sum(mcharge[start_idx:end_idx]))
                    all_frag_charges = []
                    for f in config.fragments:
                        s, e = f.charge_sum_range
                        all_frag_charges.append(float(np.sum(mcharge[s:e])))
                    charges_str = ", ".join([f"Q_{f.name}={q:+.6f}" for f, q in zip(config.fragments, all_frag_charges)])
                    print(f"  {frag.name}: E = {energy:.10f} a.u., {charges_str}")
                    
                    # Run WriteHS calculation to output oversqr.dat for CDFTB-CI
                    # Note: DFTB+ crashes after WriteHS output, which is expected behavior
                    if config.ci_enabled:
                        _, _, error_hs = run_dftb_in_subprocess(
                            qm_coords_bohr, frag_dir, write_hs=True,
                            dftb_library_path=config.dftb_library_path,
                            num_threads=config.num_threads,
                            timeout=config.timeout
                        )
                        # Ignore WriteHS errors as DFTB+ crashes after output (expected)
                
                energies.append(energy)
                charges.append(all_frag_charges)
            
            # Save CDFTB results
            energy_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for e in energies:
                energy_parts.append(f"{e:18.10f}")
            energy_file.write("  ".join(energy_parts) + "\n")
            energy_file.flush()
            
            charge_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for qs in charges:
                for q in qs:
                    charge_parts.append(f"{q:+12.6f}")
            charge_file.write("  ".join(charge_parts) + "\n")
            charge_file.flush()
            
            # Compute CDFTB-CI if enabled and CDFTB was successful
            if config.ci_enabled and cdftb_success and len(config.fragments) == 2:
                E_A = energies[0]
                E_B = energies[1]
                
                ham, J_direct, J_lowdin, ci_error = compute_cdftbci_for_frame(
                    work_dir,
                    config.fragments[0].name,
                    config.fragments[1].name,
                    E_A, E_B,
                    config.ci_N_A, config.ci_N_B
                )
                
                if ci_error:
                    print(f"  CDFTB-CI: ERROR - {ci_error}")
                    # Write NaN values
                    ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {'nan':>12s}  "
                                  f"{'nan':>16s}  {'nan':>16s}  {'nan':>10s}\n")
                    ci_sub_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 17) + "\n")
                else:
                    # Solve eigenvalue problem
                    eigenvalues, _ = solve_cdftbci_unrestricted(ham.H, ham.S)
                    
                    J_direct_meV = J_direct * 27211.386
                    J_lowdin_meV = J_lowdin * 27211.386
                    dE_eV = (eigenvalues[1] - eigenvalues[0]) * 27.211386
                    
                    print(f"  CDFTB-CI: S_AB={ham.S_AB:.6f}, J={J_lowdin_meV:.2f} meV, ΔE={dE_eV:.4f} eV")
                    
                    # Write CI results
                    ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {J_lowdin_meV:12.4f}  "
                                  f"{eigenvalues[0]:16.10f}  {eigenvalues[1]:16.10f}  {dE_eV:10.6f}\n")
                    ci_file.flush()
                    
                    # Write CI sub values
                    ci_sub_file.write(
                        f"{frame_id:5d}  {time_val:8.2f}  {E_A:14.10f}  {E_B:14.10f}  "
                        f"{ham.H_AB:14.10f}  {J_direct_meV:12.4f}  "
                        f"{ham.V_A:14.10f}  {ham.V_B:14.10f}  "
                        f"{ham.N_A:6.1f}  {ham.N_B:6.1f}  {ham.S_AB:14.10f}  "
                        f"{ham.S_AB_alpha:14.10f}  {ham.S_AB_beta:14.10f}  "
                        f"{ham.W_BA:14.10f}  {ham.W_BA_alpha:14.10f}  {ham.W_BA_beta:14.10f}  "
                        f"{ham.W_AB:14.10f}  {ham.W_AB_alpha:14.10f}  {ham.W_AB_beta:14.10f}\n"
                    )
                    ci_sub_file.flush()
            
            # Compute spin populations for both constraint states
            if cdftb_success and len(config.fragments) == 2:
                # Get number of atoms per fragment from atom_types
                n_atoms = len(atom_types)
                n_atoms_half = n_atoms // 2
                
                # State A: constraint on fragment 1
                spin_A, spin_A_error = compute_spin_for_fragment(
                    work_dir / config.fragments[0].name,
                    n_atoms_half
                )
                
                # State B: constraint on fragment 2
                spin_B, spin_B_error = compute_spin_for_fragment(
                    work_dir / config.fragments[1].name,
                    n_atoms_half
                )
                
                if spin_A_error or spin_B_error:
                    if spin_A_error:
                        print(f"  Spin (state A): ERROR - {spin_A_error[:100]}")
                    if spin_B_error:
                        print(f"  Spin (state B): ERROR - {spin_B_error[:100]}")
                    spin_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 12) + "\n")
                else:
                    # Print spin populations
                    print(f"  Spin (state A): frag1={spin_A['frag1']['spin_density']:+.4f}, "
                          f"frag2={spin_A['frag2']['spin_density']:+.4f}")
                    print(f"  Spin (state B): frag1={spin_B['frag1']['spin_density']:+.4f}, "
                          f"frag2={spin_B['frag2']['spin_density']:+.4f}")
                    
                    # Write spin values
                    spin_file.write(
                        f"{frame_id:5d}  {time_val:8.2f}  "
                        f"{spin_A['frag1']['n_alpha_frag']:14.8f}  "
                        f"{spin_A['frag1']['n_beta_frag']:14.8f}  "
                        f"{spin_A['frag1']['spin_density']:14.8f}  "
                        f"{spin_A['frag2']['n_alpha_frag']:14.8f}  "
                        f"{spin_A['frag2']['n_beta_frag']:14.8f}  "
                        f"{spin_A['frag2']['spin_density']:14.8f}  "
                        f"{spin_B['frag1']['n_alpha_frag']:14.8f}  "
                        f"{spin_B['frag1']['n_beta_frag']:14.8f}  "
                        f"{spin_B['frag1']['spin_density']:14.8f}  "
                        f"{spin_B['frag2']['n_alpha_frag']:14.8f}  "
                        f"{spin_B['frag2']['n_beta_frag']:14.8f}  "
                        f"{spin_B['frag2']['spin_density']:14.8f}\n"
                    )
                    spin_file.flush()
            
            processed_count += 1
            print()
    
    finally:
        energy_file.close()
        charge_file.close()
        if ci_file:
            ci_file.close()
        if ci_sub_file:
            ci_sub_file.close()
        if spin_file:
            spin_file.close()
        
        print("=" * 70)
        print(f"Energies saved to {config.energy_file}")
        print(f"Charges saved to {config.charge_file}")
        print(f"Spin populations saved to {config.spin_output_file}")
        if config.ci_enabled:
            print(f"CDFTB-CI results saved to {config.ci_output_file}")
            print(f"CDFTB-CI sub values saved to {config.ci_sub_file}")
        print("=" * 70)


def main():
    """Main entry point for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run CDFTB analysis with online CDFTB-CI calculation"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftb_settings.yaml",
        help="Path to YAML configuration file (default: cdftb_settings.yaml)"
    )
    args = parser.parse_args()
    
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        return 1
    
    run_cdftb_with_ci(config_path)
    return 0


if __name__ == "__main__":
    exit(main())
