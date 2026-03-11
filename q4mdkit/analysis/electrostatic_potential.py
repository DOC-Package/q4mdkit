#!/usr/bin/env python3
"""
Electrostatic Potential Analysis for QM/MM Systems

This module calculates the displacement-induced electrostatic potential
fluctuations at QM atoms due to MM point charges.

The method calculates:
    φ_a^disp(t) = φ_a(t) - φ_a^env(t)

where:
    φ_a(t) = Σ_{b∈MM} Q_b / ||r_a(t) - R_b(t)||
    φ_a^env(t) = Σ_{b∈MM} Q_b / ||r_a^ref(t) - R_b(t)||

    r_a^ref(t) = r_cm(t) + R(t)^T @ r̄_a^(0)
    r_cm(t) = (1/M) Σ_a m_a r_a(t)

φ_a(t) is the electrostatic potential at QM atom a from MM atoms.
φ_a^env(t) is the "envelope" potential - the potential if the molecule
            only underwent rigid body motion (no internal deformation).
φ_a^disp(t) is the displacement-induced fluctuation due to internal motion.

Units:
- Coordinates: Angstrom
- Charges: Elementary charge (e)
- Potential: V (Volt) = e / (4πε₀ Å) ≈ 14.3996 V

Author: q4mdkit
"""

import numpy as np
import yaml
import mdtraj as md
from typing import Tuple, Optional, List, Union
from pathlib import Path
from dataclasses import dataclass

# Import from normal_mode_analysis for consistency
from .normal_mode_analysis import (
    read_gen_file,
    get_masses,
    center_of_mass,
    remove_center_of_mass,
    kabsch_rotation,
    ATOMIC_MASSES,
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
)


def load_dcd_trajectory(
    dcd_file: str, 
    topology_file: str, 
    unit: str = 'angstrom'
) -> Tuple[np.ndarray, int]:
    """
    Load DCD trajectory file using MDTraj.
    
    Parameters
    ----------
    dcd_file : str
        Path to DCD trajectory file
    topology_file : str
        Path to topology file (PDB, GRO, etc.)
    unit : str
        Output unit: 'angstrom' or 'bohr'
        
    Returns
    -------
    coords : np.ndarray
        Coordinates array of shape (n_frames, n_atoms, 3)
    n_frames : int
        Number of frames
    """
    traj = md.load(dcd_file, top=topology_file)
    # MDTraj uses nm, convert to angstrom
    coords = traj.xyz * 10.0  # nm -> angstrom
    
    if unit.lower() == 'bohr':
        coords = coords * ANGSTROM_TO_BOHR
    
    n_frames = coords.shape[0]
    return coords, n_frames


# Unit conversion constants
# Coulomb constant k = 1/(4πε₀) in appropriate units
# e^2 / (4πε₀) = 14.3996 eV·Å
# So potential φ = Q / r has units of (e / Å) * 14.3996 V
COULOMB_CONSTANT_EV_ANGSTROM = 14.3996  # eV·Å / e^2
HARTREE_TO_EV = 27.211386245988


def read_index_file(filename: str) -> np.ndarray:
    """
    Read atom indices from an index file.
    
    Supports formats:
    - Simple list of integers (one per line or space-separated)
    - GROMACS-style index file (reads first group after header)
    
    Parameters
    ----------
    filename : str
        Path to index file
    
    Returns
    -------
    indices : np.ndarray
        0-based atom indices
    """
    indices = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
            # Skip GROMACS group headers like "[ System ]"
            if line.startswith('['):
                continue
            # Parse integers
            for token in line.split():
                try:
                    idx = int(token)
                    indices.append(idx)
                except ValueError:
                    continue
    
    indices = np.array(indices, dtype=np.int64)
    
    # Check if 1-indexed (GROMACS style) and convert to 0-indexed
    if len(indices) > 0 and indices.min() >= 1:
        indices = indices - 1
    
    return indices


def read_charges_file(filename: str, n_atoms: Optional[int] = None,
                      charge_column: int = -1) -> np.ndarray:
    """
    Read atomic charges from a text file.
    
    Supports formats:
    - Simple list of floats (one per line)
    - Multi-column format (e.g., x y z charge)
    - PCcharges.dat format (x y z charge) - 4 columns, charge in column 4
    
    Parameters
    ----------
    filename : str
        Path to charges file
    n_atoms : int, optional
        Expected number of atoms (for validation)
    charge_column : int
        Column index for charge (0-indexed). Default -1 means last column.
        For PCcharges.dat format (x y z charge), use charge_column=3 or -1.
    
    Returns
    -------
    charges : np.ndarray
        Atomic charges in elementary charge units
    """
    charges = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
            parts = line.split()
            if len(parts) == 1:
                # Single column: just charge values
                charges.append(float(parts[0]))
            elif len(parts) >= 2:
                # Multi-column format: use specified column (default: last column)
                charges.append(float(parts[charge_column]))
    
    charges = np.array(charges)
    
    if n_atoms is not None and len(charges) != n_atoms:
        raise ValueError(f"Expected {n_atoms} charges, got {len(charges)}")
    
    return charges


def compute_electrostatic_potential(
    qm_coords: np.ndarray,
    mm_coords: np.ndarray,
    mm_charges: np.ndarray
) -> np.ndarray:
    """
    Compute electrostatic potential at QM atoms from MM point charges.
    
    φ_a = Σ_{b∈MM} Q_b / ||r_a - R_b||
    
    Parameters
    ----------
    qm_coords : np.ndarray
        QM atom coordinates (N_qm, 3) in Angstrom
    mm_coords : np.ndarray
        MM atom coordinates (N_mm, 3) in Angstrom
    mm_charges : np.ndarray
        MM atom charges (N_mm,) in elementary charge units
    
    Returns
    -------
    potential : np.ndarray
        Electrostatic potential at each QM atom (N_qm,) in Volt
    """
    n_qm = qm_coords.shape[0]
    n_mm = mm_coords.shape[0]
    
    # Compute distance matrix (N_qm, N_mm)
    # r_ab = ||r_a - R_b||
    diff = qm_coords[:, np.newaxis, :] - mm_coords[np.newaxis, :, :]  # (N_qm, N_mm, 3)
    dist = np.linalg.norm(diff, axis=2)  # (N_qm, N_mm)
    
    # Avoid division by zero (shouldn't happen for QM-MM pairs)
    dist = np.maximum(dist, 1e-10)
    
    # Potential: φ_a = Σ_b Q_b / r_ab
    # Unit: e / Å, converted to V by Coulomb constant
    potential = np.sum(mm_charges[np.newaxis, :] / dist, axis=1)  # (N_qm,)
    
    # Convert to Volt: multiply by Coulomb constant
    potential = potential * COULOMB_CONSTANT_EV_ANGSTROM
    
    return potential


def compute_reference_positions(
    qm_coords: np.ndarray,
    ref_coords_centered: np.ndarray,
    masses: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute reference positions for QM atoms.
    
    r_a^ref(t) = r_cm(t) + R(t)^T @ r̄_a^(0)
    
    Parameters
    ----------
    qm_coords : np.ndarray
        QM atom coordinates at time t (N_qm, 3)
    ref_coords_centered : np.ndarray
        Centered reference coordinates r̄_a^(0) (N_qm, 3)
    masses : np.ndarray
        Atomic masses (N_qm,)
    
    Returns
    -------
    ref_positions : np.ndarray
        Reference positions r_a^ref(t) (N_qm, 3)
    r_cm : np.ndarray
        Center of mass (3,)
    R : np.ndarray
        Rotation matrix (3, 3)
    """
    # Compute center of mass of current frame
    r_cm = center_of_mass(qm_coords, masses)
    
    # Center the current coordinates
    qm_centered = qm_coords - r_cm
    
    # Find rotation matrix R that aligns reference to current:
    # R @ r̄_a^(0) ≈ qm_centered
    # So R^T @ qm_centered ≈ r̄_a^(0)
    # We want r_a^ref = r_cm + R^T @ r̄_a^(0)
    # 
    # kabsch_rotation(P, Q) finds R such that R @ P ≈ Q
    # So kabsch_rotation(ref, current) gives R such that R @ ref ≈ current
    # Therefore R^T @ current ≈ ref, which means R^T maps current → ref
    R = kabsch_rotation(ref_coords_centered, qm_centered, weights=masses)
    
    # Compute reference positions in lab frame
    # r_a^ref = r_cm + R^T @ r̄_a^(0)
    # Note: kabsch_rotation returns R such that aligned = coords @ R
    # So we need R (not R.T) here because of the row-vector convention
    ref_positions = r_cm + ref_coords_centered @ R
    
    return ref_positions, r_cm, R


def compute_displacement_potential(
    qm_coords: np.ndarray,
    mm_coords: np.ndarray,
    mm_charges: np.ndarray,
    ref_coords_centered: np.ndarray,
    masses: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute displacement-induced electrostatic potential.
    
    φ_a^disp(t) = φ_a(t) - φ_a^env(t)
    
    Parameters
    ----------
    qm_coords : np.ndarray
        QM atom coordinates (N_qm, 3)
    mm_coords : np.ndarray
        MM atom coordinates (N_mm, 3)
    mm_charges : np.ndarray
        MM atom charges (N_mm,)
    ref_coords_centered : np.ndarray
        Centered reference coordinates (N_qm, 3)
    masses : np.ndarray
        Atomic masses (N_qm,)
    
    Returns
    -------
    phi_disp : np.ndarray
        Displacement potential φ_a^disp (N_qm,) in Volt
    phi : np.ndarray
        Total potential φ_a (N_qm,) in Volt
    phi_env : np.ndarray
        Envelope potential φ_a^env (N_qm,) in Volt
    """
    # Compute reference positions
    ref_positions, r_cm, R = compute_reference_positions(
        qm_coords, ref_coords_centered, masses
    )
    
    # Compute potentials
    phi = compute_electrostatic_potential(qm_coords, mm_coords, mm_charges)
    phi_env = compute_electrostatic_potential(ref_positions, mm_coords, mm_charges)
    
    # Displacement potential
    phi_disp = phi - phi_env
    
    return phi_disp, phi, phi_env


class ElectrostaticPotentialAnalyzer:
    """
    Analyzer for computing displacement-induced electrostatic potentials
    from QM/MM molecular dynamics trajectories.
    
    Attributes
    ----------
    ref_coords : np.ndarray
        Reference structure coordinates (N_qm, 3)
    ref_coords_centered : np.ndarray
        Reference structure centered at COM (N_qm, 3)
    atom_types : list
        Element symbols for QM atoms
    masses : np.ndarray
        Atomic masses (N_qm,)
    qm_indices : np.ndarray
        Indices of QM atoms in trajectory
    mm_indices : np.ndarray
        Indices of MM atoms in trajectory
    mm_charges : np.ndarray
        MM atom charges
    """
    
    def __init__(
        self,
        ref_gen_file: str,
        qm_index_file: str,
        mm_charges_file: str,
        total_atoms: int,
        charge_column: int = -1
    ):
        """
        Initialize the analyzer.
        
        Parameters
        ----------
        ref_gen_file : str
            Path to reference structure (.gen file)
        qm_index_file : str
            Path to QM atom index file
        mm_charges_file : str
            Path to MM charges file
        total_atoms : int
            Total number of atoms in trajectory
        charge_column : int
            Column index for charge in mm_charges file (-1 = last column).
            For PCcharges.dat format (x y z charge), use -1 or 3.
        """
        # Load reference structure
        self.ref_coords, self.atom_types = read_gen_file(ref_gen_file, unit='angstrom')
        self.masses = get_masses(self.atom_types)
        
        # Center reference structure
        self.ref_coords_centered = remove_center_of_mass(self.ref_coords, self.masses)
        
        # Load QM indices
        self.qm_indices = read_index_file(qm_index_file)
        
        # Validate QM atom count
        if len(self.qm_indices) != len(self.atom_types):
            raise ValueError(
                f"QM index count ({len(self.qm_indices)}) doesn't match "
                f"reference structure atoms ({len(self.atom_types)})"
            )
        
        # Compute MM indices (all atoms not in QM)
        all_indices = np.arange(total_atoms)
        self.mm_indices = np.setdiff1d(all_indices, self.qm_indices)
        
        # Load MM charges
        self.mm_charges = read_charges_file(
            mm_charges_file, n_atoms=len(self.mm_indices), charge_column=charge_column
        )
        
        # Store dimensions
        self.n_qm = len(self.qm_indices)
        self.n_mm = len(self.mm_indices)
        self.total_atoms = total_atoms
    
    def process_frame(
        self,
        coords: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Process a single frame and compute potentials.
        
        Parameters
        ----------
        coords : np.ndarray
            All atom coordinates for this frame (N_total, 3)
        
        Returns
        -------
        phi_disp : np.ndarray
            Displacement potential (N_qm,)
        phi : np.ndarray
            Total potential (N_qm,)
        phi_env : np.ndarray
            Envelope potential (N_qm,)
        r_cm : np.ndarray
            Center of mass (3,)
        R : np.ndarray
            Rotation matrix (3, 3)
        """
        # Extract QM and MM coordinates
        qm_coords = coords[self.qm_indices]
        mm_coords = coords[self.mm_indices]
        
        # Compute reference positions
        ref_positions, r_cm, R = compute_reference_positions(
            qm_coords, self.ref_coords_centered, self.masses
        )
        
        # Compute potentials
        phi = compute_electrostatic_potential(qm_coords, mm_coords, self.mm_charges)
        phi_env = compute_electrostatic_potential(ref_positions, mm_coords, self.mm_charges)
        phi_disp = phi - phi_env
        
        return phi_disp, phi, phi_env, r_cm, R
    
    def process_trajectory(
        self,
        trajectory: np.ndarray,
        verbose: bool = True
    ) -> dict:
        """
        Process entire trajectory and compute potentials for all frames.
        
        Parameters
        ----------
        trajectory : np.ndarray
            Trajectory coordinates (n_frames, N_total, 3) in Angstrom
        verbose : bool
            Print progress information
        
        Returns
        -------
        results : dict
            Dictionary containing:
            - 'phi_disp': Displacement potentials (n_frames, N_qm)
            - 'phi': Total potentials (n_frames, N_qm)
            - 'phi_env': Envelope potentials (n_frames, N_qm)
            - 'r_cm': Center of mass trajectory (n_frames, 3)
            - 'rotation_angles': Rotation angles in degrees (n_frames,)
        """
        n_frames = trajectory.shape[0]
        
        # Allocate output arrays
        phi_disp = np.zeros((n_frames, self.n_qm))
        phi = np.zeros((n_frames, self.n_qm))
        phi_env = np.zeros((n_frames, self.n_qm))
        r_cm = np.zeros((n_frames, 3))
        rotation_angles = np.zeros(n_frames)
        
        if verbose:
            print(f"Processing {n_frames} frames...")
        
        for i in range(n_frames):
            if verbose and (i + 1) % 1000 == 0:
                print(f"  Frame {i + 1}/{n_frames}")
            
            phi_disp[i], phi[i], phi_env[i], r_cm[i], R = self.process_frame(
                trajectory[i]
            )
            
            # Calculate rotation angle from rotation matrix
            # θ = arccos((tr(R) - 1) / 2)
            trace = np.trace(R)
            cos_theta = np.clip((trace - 1) / 2, -1, 1)
            rotation_angles[i] = np.degrees(np.arccos(cos_theta))
        
        if verbose:
            print("Done.")
        
        return {
            'phi_disp': phi_disp,
            'phi': phi,
            'phi_env': phi_env,
            'r_cm': r_cm,
            'rotation_angles': rotation_angles,
        }
    
    def save_results(
        self,
        results: dict,
        output_prefix: str,
        save_individual: bool = True
    ):
        """
        Save analysis results to files.
        
        Parameters
        ----------
        results : dict
            Results from process_trajectory()
        output_prefix : str
            Prefix for output files
        save_individual : bool
            Save individual atom potentials (can be large)
        """
        output_dir = Path(output_prefix).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        prefix = Path(output_prefix).name
        
        # Save displacement potential (primary output)
        np.savetxt(
            output_dir / f"{prefix}_phi_disp.dat",
            results['phi_disp'],
            header=f"Displacement potential φ_a^disp(t) [V] for {self.n_qm} QM atoms\n"
                   f"Columns: atom 0, atom 1, ..., atom {self.n_qm - 1}",
            fmt='%.8e'
        )
        
        # Save sum/mean over all QM atoms
        phi_disp_sum = np.sum(results['phi_disp'], axis=1)
        phi_disp_mean = np.mean(results['phi_disp'], axis=1)
        np.savetxt(
            output_dir / f"{prefix}_phi_disp_summary.dat",
            np.column_stack([phi_disp_sum, phi_disp_mean]),
            header="Displacement potential summary\n"
                   "Column 1: sum over QM atoms [V]\n"
                   "Column 2: mean over QM atoms [V]",
            fmt='%.8e'
        )
        
        if save_individual:
            # Save total potential
            np.savetxt(
                output_dir / f"{prefix}_phi.dat",
                results['phi'],
                header=f"Total potential φ_a(t) [V] for {self.n_qm} QM atoms",
                fmt='%.8e'
            )
            
            # Save envelope potential
            np.savetxt(
                output_dir / f"{prefix}_phi_env.dat",
                results['phi_env'],
                header=f"Envelope potential φ_a^env(t) [V] for {self.n_qm} QM atoms",
                fmt='%.8e'
            )
        
        # Save center of mass trajectory
        np.savetxt(
            output_dir / f"{prefix}_r_cm.dat",
            results['r_cm'],
            header="Center of mass trajectory [Å]\nColumns: x, y, z",
            fmt='%.6f'
        )
        
        # Save rotation angles
        np.savetxt(
            output_dir / f"{prefix}_rotation.dat",
            results['rotation_angles'],
            header="Rotation angle from reference [degrees]",
            fmt='%.6f'
        )
        
        print(f"Results saved with prefix: {output_prefix}")


@dataclass
class ElectrostaticPotentialConfig:
    """Configuration for electrostatic potential analysis."""
    # Input paths
    trajectory: Path
    topology: Path
    reference: Path
    qm_indices: Path
    mm_charges: Path
    
    # Output settings
    output_dir: Path = Path("output_potential")
    output_prefix: str = "potential"
    save_individual: bool = True
    
    # Frame selection
    start_frame: int = 0
    end_frame: Optional[int] = None
    stride: int = 1
    
    # Analysis settings
    verbose: bool = True
    
    # File format options
    charge_column: int = -1  # Column index for charge in mm_charges file (-1 = last column)


def load_config(config_path: str) -> ElectrostaticPotentialConfig:
    """
    Load configuration from YAML file.
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file
    
    Returns
    -------
    config : ElectrostaticPotentialConfig
        Configuration object
    """
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse input paths (relative to config file)
    input_cfg = data.get('input', {})
    
    required_inputs = ['trajectory', 'topology', 'reference', 'qm_indices', 'mm_charges']
    for key in required_inputs:
        if key not in input_cfg:
            raise ValueError(f"Missing required input: {key}")
    
    trajectory = base_dir / input_cfg['trajectory']
    topology = base_dir / input_cfg['topology']
    reference = base_dir / input_cfg['reference']
    qm_indices = base_dir / input_cfg['qm_indices']
    mm_charges = base_dir / input_cfg['mm_charges']
    
    # Parse output settings
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output_potential')
    output_prefix = output_cfg.get('prefix', 'potential')
    save_individual = output_cfg.get('save_individual', True)
    
    # Parse frame selection
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    end_frame = frames_cfg.get('end', None)
    stride = frames_cfg.get('stride', 1)
    
    # Parse analysis settings
    analysis_cfg = data.get('analysis', {})
    verbose = analysis_cfg.get('verbose', True)
    charge_column = input_cfg.get('charge_column', -1)
    
    return ElectrostaticPotentialConfig(
        trajectory=trajectory,
        topology=topology,
        reference=reference,
        qm_indices=qm_indices,
        mm_charges=mm_charges,
        output_dir=output_dir,
        output_prefix=output_prefix,
        save_individual=save_individual,
        start_frame=start_frame,
        end_frame=end_frame,
        stride=stride,
        verbose=verbose,
        charge_column=charge_column,
    )


def run_from_config(config: ElectrostaticPotentialConfig) -> dict:
    """
    Run electrostatic potential analysis from configuration.
    
    Parameters
    ----------
    config : ElectrostaticPotentialConfig
        Configuration object
    
    Returns
    -------
    n_frames : int
        Number of frames processed
    """
    # Ensure output directory exists first
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    if config.verbose:
        print(f"Output directory: {config.output_dir}")
        print(f"Loading trajectory: {config.trajectory}")
        print(f"Topology: {config.topology}")
    
    # Load trajectory
    trajectory, n_frames = load_dcd_trajectory(
        str(config.trajectory), str(config.topology), unit='angstrom'
    )
    total_atoms = trajectory.shape[1]
    
    if config.verbose:
        print(f"Loaded {n_frames} frames, {total_atoms} atoms")
    
    # Apply frame selection
    start = config.start_frame
    end = config.end_frame if config.end_frame is not None else n_frames
    stride = config.stride
    
    trajectory = trajectory[start:end:stride]
    n_frames_selected = trajectory.shape[0]
    
    if config.verbose:
        print(f"Selected {n_frames_selected} frames (start={start}, end={end}, stride={stride})")
    
    # Create analyzer
    if config.verbose:
        print(f"Reference structure: {config.reference}")
        print(f"QM indices: {config.qm_indices}")
        print(f"MM charges: {config.mm_charges}")
    
    analyzer = ElectrostaticPotentialAnalyzer(
        str(config.reference),
        str(config.qm_indices),
        str(config.mm_charges),
        total_atoms,
        charge_column=config.charge_column
    )
    
    if config.verbose:
        print(f"QM atoms: {analyzer.n_qm}")
        print(f"MM atoms: {analyzer.n_mm}")
    
    # Prepare output files
    prefix = config.output_prefix
    output_dir = config.output_dir
    
    # Open all output files for streaming write
    f_phi_disp = open(output_dir / f"{prefix}_phi_disp.dat", 'w')
    f_summary = open(output_dir / f"{prefix}_phi_disp_summary.dat", 'w')
    f_r_cm = open(output_dir / f"{prefix}_r_cm.dat", 'w')
    f_rotation = open(output_dir / f"{prefix}_rotation.dat", 'w')
    
    if config.save_individual:
        f_phi = open(output_dir / f"{prefix}_phi.dat", 'w')
        f_phi_env = open(output_dir / f"{prefix}_phi_env.dat", 'w')
    else:
        f_phi = None
        f_phi_env = None
    
    try:
        # Write headers
        f_phi_disp.write(f"# Displacement potential φ_a^disp(t) [V] for {analyzer.n_qm} QM atoms\n")
        f_phi_disp.write(f"# Columns: atom 0, atom 1, ..., atom {analyzer.n_qm - 1}\n")
        
        f_summary.write("# Displacement potential summary\n")
        f_summary.write("# Column 1: sum over QM atoms [V]\n")
        f_summary.write("# Column 2: mean over QM atoms [V]\n")
        
        f_r_cm.write("# Center of mass trajectory [Å]\n")
        f_r_cm.write("# Columns: x, y, z\n")
        
        f_rotation.write("# Rotation angle from reference [degrees]\n")
        
        if f_phi is not None:
            f_phi.write(f"# Total potential φ_a(t) [V] for {analyzer.n_qm} QM atoms\n")
        if f_phi_env is not None:
            f_phi_env.write(f"# Envelope potential φ_a^env(t) [V] for {analyzer.n_qm} QM atoms\n")
        
        # Process trajectory frame by frame
        if config.verbose:
            print(f"Processing {n_frames_selected} frames...")
        
        for i in range(n_frames_selected):
            if config.verbose and (i + 1) % 1000 == 0:
                print(f"  Frame {i + 1}/{n_frames_selected}")
            
            # Process frame
            phi_disp, phi, phi_env, r_cm, R = analyzer.process_frame(trajectory[i])
            
            # Calculate rotation angle
            trace = np.trace(R)
            cos_theta = np.clip((trace - 1) / 2, -1, 1)
            rotation_angle = np.degrees(np.arccos(cos_theta))
            
            # Write results immediately
            f_phi_disp.write(" ".join(f"{v:.8e}" for v in phi_disp) + "\n")
            
            phi_disp_sum = np.sum(phi_disp)
            phi_disp_mean = np.mean(phi_disp)
            f_summary.write(f"{phi_disp_sum:.8e} {phi_disp_mean:.8e}\n")
            
            f_r_cm.write(f"{r_cm[0]:.6f} {r_cm[1]:.6f} {r_cm[2]:.6f}\n")
            f_rotation.write(f"{rotation_angle:.6f}\n")
            
            if f_phi is not None:
                f_phi.write(" ".join(f"{v:.8e}" for v in phi) + "\n")
            if f_phi_env is not None:
                f_phi_env.write(" ".join(f"{v:.8e}" for v in phi_env) + "\n")
            
            # Flush periodically to ensure data is written
            if (i + 1) % 100 == 0:
                f_phi_disp.flush()
                f_summary.flush()
                f_r_cm.flush()
                f_rotation.flush()
                if f_phi is not None:
                    f_phi.flush()
                if f_phi_env is not None:
                    f_phi_env.flush()
        
        if config.verbose:
            print("Done.")
    
    finally:
        # Close all files
        f_phi_disp.close()
        f_summary.close()
        f_r_cm.close()
        f_rotation.close()
        if f_phi is not None:
            f_phi.close()
        if f_phi_env is not None:
            f_phi_env.close()
    
    print(f"Results saved to: {output_dir}")
    
    return n_frames_selected


def analyze_electrostatic_potential(
    dcd_file: str,
    top_file: str,
    ref_gen_file: str,
    qm_index_file: str,
    mm_charges_file: str,
    output_prefix: str = "potential",
    save_individual: bool = True,
    verbose: bool = True,
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    stride: int = 1,
    charge_column: int = -1
) -> int:
    """
    Main function to analyze electrostatic potentials from trajectory.
    
    Parameters
    ----------
    dcd_file : str
        Path to DCD trajectory file
    top_file : str
        Path to topology file (PDB)
    ref_gen_file : str
        Path to reference structure (.gen file)
    qm_index_file : str
        Path to QM atom index file
    mm_charges_file : str
        Path to MM charges file
    output_prefix : str
        Prefix for output files
    save_individual : bool
        Save individual atom potentials
    verbose : bool
        Print progress information
    start_frame : int
        First frame to process (0-indexed)
    end_frame : int, optional
        Last frame to process (exclusive), None for all
    stride : int
        Frame stride
    charge_column : int
        Column index for charge in mm_charges file (-1 = last column)
    
    Returns
    -------
    n_frames : int
        Number of frames processed
    """
    # Ensure output directory exists
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = Path(output_prefix).name
    
    if verbose:
        print(f"Output directory: {output_dir}")
        print(f"Loading trajectory: {dcd_file}")
        print(f"Topology: {top_file}")
    
    # Load trajectory
    trajectory, n_frames = load_dcd_trajectory(dcd_file, top_file, unit='angstrom')
    total_atoms = trajectory.shape[1]
    
    if verbose:
        print(f"Loaded {n_frames} frames, {total_atoms} atoms")
    
    # Apply frame selection
    end = end_frame if end_frame is not None else n_frames
    trajectory = trajectory[start_frame:end:stride]
    n_frames_selected = trajectory.shape[0]
    
    if verbose:
        print(f"Selected {n_frames_selected} frames (start={start_frame}, end={end}, stride={stride})")
    
    # Create analyzer
    if verbose:
        print(f"Reference structure: {ref_gen_file}")
        print(f"QM indices: {qm_index_file}")
        print(f"MM charges: {mm_charges_file}")
    
    analyzer = ElectrostaticPotentialAnalyzer(
        ref_gen_file, qm_index_file, mm_charges_file, total_atoms,
        charge_column=charge_column
    )
    
    if verbose:
        print(f"QM atoms: {analyzer.n_qm}")
        print(f"MM atoms: {analyzer.n_mm}")
    
    # Open all output files for streaming write
    f_phi_disp = open(output_dir / f"{prefix}_phi_disp.dat", 'w')
    f_summary = open(output_dir / f"{prefix}_phi_disp_summary.dat", 'w')
    f_r_cm = open(output_dir / f"{prefix}_r_cm.dat", 'w')
    f_rotation = open(output_dir / f"{prefix}_rotation.dat", 'w')
    
    if save_individual:
        f_phi = open(output_dir / f"{prefix}_phi.dat", 'w')
        f_phi_env = open(output_dir / f"{prefix}_phi_env.dat", 'w')
    else:
        f_phi = None
        f_phi_env = None
    
    try:
        # Write headers
        f_phi_disp.write(f"# Displacement potential φ_a^disp(t) [V] for {analyzer.n_qm} QM atoms\n")
        f_phi_disp.write(f"# Columns: atom 0, atom 1, ..., atom {analyzer.n_qm - 1}\n")
        
        f_summary.write("# Displacement potential summary\n")
        f_summary.write("# Column 1: sum over QM atoms [V]\n")
        f_summary.write("# Column 2: mean over QM atoms [V]\n")
        
        f_r_cm.write("# Center of mass trajectory [Å]\n")
        f_r_cm.write("# Columns: x, y, z\n")
        
        f_rotation.write("# Rotation angle from reference [degrees]\n")
        
        if f_phi is not None:
            f_phi.write(f"# Total potential φ_a(t) [V] for {analyzer.n_qm} QM atoms\n")
        if f_phi_env is not None:
            f_phi_env.write(f"# Envelope potential φ_a^env(t) [V] for {analyzer.n_qm} QM atoms\n")
        
        # Process trajectory frame by frame
        if verbose:
            print(f"Processing {n_frames_selected} frames...")
        
        for i in range(n_frames_selected):
            if verbose and (i + 1) % 1000 == 0:
                print(f"  Frame {i + 1}/{n_frames_selected}")
            
            # Process frame
            phi_disp, phi, phi_env, r_cm, R = analyzer.process_frame(trajectory[i])
            
            # Calculate rotation angle
            trace = np.trace(R)
            cos_theta = np.clip((trace - 1) / 2, -1, 1)
            rotation_angle = np.degrees(np.arccos(cos_theta))
            
            # Write results immediately
            f_phi_disp.write(" ".join(f"{v:.8e}" for v in phi_disp) + "\n")
            
            phi_disp_sum = np.sum(phi_disp)
            phi_disp_mean = np.mean(phi_disp)
            f_summary.write(f"{phi_disp_sum:.8e} {phi_disp_mean:.8e}\n")
            
            f_r_cm.write(f"{r_cm[0]:.6f} {r_cm[1]:.6f} {r_cm[2]:.6f}\n")
            f_rotation.write(f"{rotation_angle:.6f}\n")
            
            if f_phi is not None:
                f_phi.write(" ".join(f"{v:.8e}" for v in phi) + "\n")
            if f_phi_env is not None:
                f_phi_env.write(" ".join(f"{v:.8e}" for v in phi_env) + "\n")
            
            # Flush periodically
            if (i + 1) % 100 == 0:
                f_phi_disp.flush()
                f_summary.flush()
                f_r_cm.flush()
                f_rotation.flush()
                if f_phi is not None:
                    f_phi.flush()
                if f_phi_env is not None:
                    f_phi_env.flush()
        
        if verbose:
            print("Done.")
    
    finally:
        # Close all files
        f_phi_disp.close()
        f_summary.close()
        f_r_cm.close()
        f_rotation.close()
        if f_phi is not None:
            f_phi.close()
        if f_phi_env is not None:
            f_phi_env.close()
    
    print(f"Results saved with prefix: {output_prefix}")
    
    return n_frames_selected


def run_potential_analysis(config_path: str) -> int:
    """
    Run electrostatic potential analysis from YAML config file.
    
    This is the main entry point for running analysis from a configuration file.
    Results are written to files frame by frame (streaming mode).
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file
    
    Returns
    -------
    n_frames : int
        Number of frames processed
    
    Example
    -------
    >>> from q4mdkit.analysis.electrostatic_potential import run_potential_analysis
    >>> n_frames = run_potential_analysis("potential_settings.yaml")
    """
    config = load_config(config_path)
    return run_from_config(config)
