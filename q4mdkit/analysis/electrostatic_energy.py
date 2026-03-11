#!/usr/bin/env python3
"""
Electrostatic Energy Analysis for QM/MM Systems

This module calculates the electrostatic interaction energy between
QM atoms and MM point charges for each snapshot.

The method calculates:
    E_a(t) = q_a * φ_a(t) = q_a * Σ_{b∈MM} Q_b / ||r_a(t) - R_b(t)||
    
    E_total(t) = Σ_{a∈QM} E_a(t)

With optional QM molecule fixing:
    E_a^fixed(t) = q_a * Σ_{b∈MM} Q_b / ||r_a^ref(t) - R_b(t)||
    
    where r_a^ref(t) = r_cm(t) + R(t)^T @ r̄_a^(0)  (rigid body motion only)

Units:
- Coordinates: Angstrom
- Charges: Elementary charge (e)
- Energy: eV

Author: q4mdkit
"""

import numpy as np
import yaml
import mdtraj as md
from typing import Tuple, Optional, List
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
)

# Unit conversion constants
COULOMB_CONSTANT_EV_ANGSTROM = 14.3996  # eV·Å / e^2


def read_xyz_file(filename: str) -> Tuple[np.ndarray, List[str]]:
    """
    Read XYZ file format.
    
    Parameters
    ----------
    filename : str
        Path to XYZ file
        
    Returns
    -------
    coords : np.ndarray
        Coordinates (N, 3) in Angstrom
    atom_types : list
        Element symbols
    """
    coords = []
    atom_types = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    n_atoms = int(lines[0].strip())
    # Line 1 is comment, lines 2+ are atoms
    
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        atom_types.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return np.array(coords), atom_types


# Bohr per Angstrom conversion
BOHR_PER_ANG = 1.8897259886


def read_structure_file(
    filename: str,
    output_unit: str = "angstrom",
) -> Tuple[np.ndarray, List[str]]:
    """
    Read structure file (XYZ or GEN format).
    
    Parameters
    ----------
    filename : str
        Path to structure file (.xyz or .gen)
    output_unit : str
        Output unit: "angstrom" or "bohr"
    
    Returns
    -------
    coords : np.ndarray
        Coordinates (N, 3)
    atom_types : list
        Element symbols
    """
    filepath = Path(filename)
    suffix = filepath.suffix.lower()
    
    if suffix == '.xyz':
        coords_ang, atom_types = read_xyz_file(filename)
    elif suffix == '.gen':
        coords_ang, atom_types = read_gen_file(filename)
    else:
        raise ValueError(f"Unsupported file format: {suffix}. Use .xyz or .gen")
    
    if output_unit == "bohr":
        return coords_ang * BOHR_PER_ANG, atom_types
    elif output_unit == "angstrom":
        return coords_ang, atom_types
    else:
        raise ValueError(f"Unknown unit: {output_unit}. Use 'angstrom' or 'bohr'")


def load_dcd_trajectory(
    dcd_file: str, 
    topology_file: str
) -> Tuple[np.ndarray, int]:
    """
    Load DCD trajectory file using MDTraj.
    
    Returns coordinates in Angstrom.
    """
    traj = md.load(dcd_file, top=topology_file)
    coords = traj.xyz * 10.0  # nm -> angstrom
    n_frames = coords.shape[0]
    return coords, n_frames


def read_index_file(filename: str) -> np.ndarray:
    """
    Read atom indices from an index file.
    
    Supports:
    - Simple list of integers (one per line or space-separated)
    - GROMACS-style index file (reads first group)
    
    Returns 0-based indices.
    """
    indices = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
            if line.startswith('['):
                continue
            for token in line.split():
                try:
                    idx = int(token)
                    indices.append(idx)
                except ValueError:
                    continue
    
    indices = np.array(indices, dtype=np.int64)
    
    # Convert to 0-indexed if 1-indexed
    if len(indices) > 0 and indices.min() >= 1:
        indices = indices - 1
    
    return indices


def read_charges_file(filename: str, charge_column: int = -1) -> np.ndarray:
    """
    Read atomic charges from a text file.
    
    Supports:
    - Single column (just charges)
    - Multi-column (e.g., x y z charge) - use charge_column to specify
    
    Parameters
    ----------
    filename : str
        Path to charges file
    charge_column : int
        Column index for charge (0-indexed). Default -1 = last column.
    
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
                charges.append(float(parts[0]))
            elif len(parts) >= 2:
                charges.append(float(parts[charge_column]))
    
    return np.array(charges)


def extract_mm_charges_from_pccharges(
    pccharges_file: str,
    qm_indices: np.ndarray,
    mm_indices: np.ndarray,
    total_atoms: int,
    charge_column: int = -1
) -> np.ndarray:
    """
    Extract MM charges from PCcharges.dat for specific MM atoms.
    
    PCcharges.dat contains charges for all MM atoms (excluding QM atoms),
    ordered by system index. This function extracts charges for a specific
    subset of MM atoms (e.g., a single MM molecule).
    
    Parameters
    ----------
    pccharges_file : str
        Path to PCcharges.dat file
    qm_indices : np.ndarray
        QM atom indices in system (0-based)
    mm_indices : np.ndarray
        MM atom indices to extract (0-based, system indices)
    total_atoms : int
        Total number of atoms in system
    charge_column : int
        Column index for charge (-1 = last column)
    
    Returns
    -------
    charges : np.ndarray
        Charges for the specified MM atoms
    """
    # Read all charges from PCcharges.dat
    all_mm_charges = read_charges_file(pccharges_file, charge_column)
    
    # All MM indices in system (sorted)
    all_indices = np.arange(total_atoms)
    all_mm_indices = np.setdiff1d(all_indices, qm_indices)
    
    if len(all_mm_charges) != len(all_mm_indices):
        raise ValueError(
            f"PCcharges.dat has {len(all_mm_charges)} charges, "
            f"but system has {len(all_mm_indices)} MM atoms"
        )
    
    # Create mapping: system index -> PCcharges row
    system_to_pcrow = {idx: row for row, idx in enumerate(all_mm_indices)}
    
    # Extract charges for requested MM atoms
    charges = np.zeros(len(mm_indices))
    for i, mm_idx in enumerate(mm_indices):
        if mm_idx not in system_to_pcrow:
            raise ValueError(
                f"MM index {mm_idx} is not in MM atom list "
                f"(it may be a QM atom or out of range)"
            )
        row = system_to_pcrow[mm_idx]
        charges[i] = all_mm_charges[row]
    
    return charges


def compute_electrostatic_energy(
    qm_coords: np.ndarray,
    qm_charges: np.ndarray,
    mm_coords: np.ndarray,
    mm_charges: np.ndarray
) -> Tuple[np.ndarray, float]:
    """
    Compute electrostatic interaction energy between QM and MM atoms.
    
    E_a = q_a * Σ_{b∈MM} Q_b / ||r_a - R_b||
    
    Parameters
    ----------
    qm_coords : np.ndarray
        QM atom coordinates (N_qm, 3) in Angstrom
    qm_charges : np.ndarray
        QM atom charges (N_qm,) in elementary charge units
    mm_coords : np.ndarray
        MM atom coordinates (N_mm, 3) in Angstrom
    mm_charges : np.ndarray
        MM atom charges (N_mm,) in elementary charge units
    
    Returns
    -------
    energy_per_atom : np.ndarray
        Electrostatic energy at each QM atom (N_qm,) in eV
    total_energy : float
        Total electrostatic energy in eV
    """
    n_qm = qm_coords.shape[0]
    
    # Compute distance matrix (N_qm, N_mm)
    diff = qm_coords[:, np.newaxis, :] - mm_coords[np.newaxis, :, :]  # (N_qm, N_mm, 3)
    dist = np.linalg.norm(diff, axis=2)  # (N_qm, N_mm)
    
    # Avoid division by zero
    dist = np.maximum(dist, 1e-10)
    
    # Potential at each QM atom: φ_a = Σ_b Q_b / r_ab
    potential = np.sum(mm_charges[np.newaxis, :] / dist, axis=1)  # (N_qm,)
    
    # Energy: E_a = q_a * φ_a * k (convert to eV)
    energy_per_atom = qm_charges * potential * COULOMB_CONSTANT_EV_ANGSTROM
    
    total_energy = np.sum(energy_per_atom)
    
    return energy_per_atom, total_energy


def compute_reference_positions(
    qm_coords: np.ndarray,
    ref_coords_centered: np.ndarray,
    masses: np.ndarray
) -> np.ndarray:
    """
    Compute reference positions for QM atoms (rigid body motion only).
    
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
    """
    # Compute center of mass of current frame
    r_cm = center_of_mass(qm_coords, masses)
    
    # Center the current coordinates
    qm_centered = qm_coords - r_cm
    
    # Find rotation matrix
    R = kabsch_rotation(ref_coords_centered, qm_centered, weights=masses)
    
    # Compute reference positions in lab frame
    ref_positions = r_cm + ref_coords_centered @ R
    
    return ref_positions


@dataclass
class ElectrostaticEnergyConfig:
    """Configuration for electrostatic energy analysis."""
    # Input paths
    trajectory: Path
    topology: Path
    qm_indices: Path
    qm_charges: Path
    mm_indices: Path           # MM molecule indices (required)
    
    # MM charges: either direct file or extract from PCcharges.dat
    mm_charges: Optional[Path] = None      # Direct MM charges file
    pccharges_file: Optional[Path] = None  # Full PCcharges.dat (auto-extract)
    
    # Optional: reference structure for fixing QM molecule
    reference: Optional[Path] = None  # XYZ file for reference structure
    fix_qm: bool = False
    
    # Output settings
    output_dir: Path = Path("output_elstat_energy")
    output_prefix: str = "elstat_energy"
    
    # Frame selection
    start_frame: int = 0
    end_frame: Optional[int] = None
    stride: int = 1
    
    # Analysis settings
    verbose: bool = True
    charge_column: int = -1


def load_config(config_path: str) -> ElectrostaticEnergyConfig:
    """Load configuration from YAML file."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    
    input_cfg = data.get('input', {})
    
    trajectory = base_dir / input_cfg['trajectory']
    topology = base_dir / input_cfg['topology']
    qm_indices = base_dir / input_cfg['qm_indices']
    qm_charges = base_dir / input_cfg['qm_charges']
    mm_indices = base_dir / input_cfg['mm_indices']
    
    # MM charges: either direct file or PCcharges.dat
    mm_charges = None
    pccharges_file = None
    if 'mm_charges' in input_cfg:
        mm_charges = base_dir / input_cfg['mm_charges']
    if 'pccharges_file' in input_cfg:
        pccharges_file = base_dir / input_cfg['pccharges_file']
    
    if mm_charges is None and pccharges_file is None:
        raise ValueError("Either 'mm_charges' or 'pccharges_file' must be specified")
    
    reference = None
    if 'reference' in input_cfg:
        reference = base_dir / input_cfg['reference']
    
    fix_qm = input_cfg.get('fix_qm', False)
    charge_column = input_cfg.get('charge_column', -1)
    
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output_elstat_energy')
    output_prefix = output_cfg.get('prefix', 'elstat_energy')
    
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    end_frame = frames_cfg.get('end', None)
    stride = frames_cfg.get('stride', 1)
    
    analysis_cfg = data.get('analysis', {})
    verbose = analysis_cfg.get('verbose', True)
    
    return ElectrostaticEnergyConfig(
        trajectory=trajectory,
        topology=topology,
        qm_indices=qm_indices,
        qm_charges=qm_charges,
        mm_indices=mm_indices,
        mm_charges=mm_charges,
        pccharges_file=pccharges_file,
        reference=reference,
        fix_qm=fix_qm,
        output_dir=output_dir,
        output_prefix=output_prefix,
        start_frame=start_frame,
        end_frame=end_frame,
        stride=stride,
        verbose=verbose,
        charge_column=charge_column,
    )


def run_electrostatic_energy_analysis(config: ElectrostaticEnergyConfig) -> int:
    """
    Run electrostatic energy analysis.
    
    Parameters
    ----------
    config : ElectrostaticEnergyConfig
        Configuration object
    
    Returns
    -------
    n_frames : int
        Number of frames processed
    """
    # Ensure output directory exists
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    if config.verbose:
        print("=" * 70)
        print("Electrostatic Energy Analysis")
        print("=" * 70)
        print(f"Output directory: {config.output_dir}")
        print(f"Loading trajectory: {config.trajectory}")
    
    # Load trajectory
    trajectory, n_frames = load_dcd_trajectory(
        str(config.trajectory), str(config.topology)
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
    
    # Load indices
    qm_indices = read_index_file(str(config.qm_indices))
    n_qm = len(qm_indices)
    
    # Load MM indices from file
    mm_indices = read_index_file(str(config.mm_indices))
    n_mm = len(mm_indices)
    
    if config.verbose:
        print(f"QM atoms: {n_qm}")
        print(f"MM atoms: {n_mm}")
    
    # Load charges
    qm_charges = read_charges_file(str(config.qm_charges), config.charge_column)
    
    # Load MM charges - either from direct file or extract from PCcharges.dat
    if config.pccharges_file is not None:
        mm_charges = extract_mm_charges_from_pccharges(
            str(config.pccharges_file),
            mm_indices,
            qm_indices,
        )
        if config.verbose:
            print(f"Extracted MM charges from PCcharges.dat: {config.pccharges_file}")
    elif config.mm_charges is not None:
        mm_charges = read_charges_file(str(config.mm_charges), config.charge_column)
    else:
        raise ValueError("Either mm_charges or pccharges_file must be specified")
    
    if len(qm_charges) != n_qm:
        raise ValueError(f"QM charges count ({len(qm_charges)}) != QM atoms ({n_qm})")
    if len(mm_charges) != n_mm:
        raise ValueError(f"MM charges count ({len(mm_charges)}) != MM atoms ({n_mm})")
    
    if config.verbose:
        print(f"QM total charge: {np.sum(qm_charges):.4f} e")
        print(f"MM total charge: {np.sum(mm_charges):.4f} e")
    
    # Load reference structure if fixing QM
    ref_coords_centered = None
    masses = None
    
    if config.fix_qm:
        if config.reference is None:
            raise ValueError("Reference structure (XYZ file) required when fix_qm=True")
        
        ref_coords, atom_types = read_xyz_file(str(config.reference))
        masses = get_masses(atom_types)
        
        if len(ref_coords) != n_qm:
            raise ValueError(f"Reference atoms ({len(ref_coords)}) != QM atoms ({n_qm})")
        
        ref_coords_centered = remove_center_of_mass(ref_coords, masses)
        
        if config.verbose:
            print(f"Reference structure: {config.reference}")
            print("QM molecule will be FIXED (rigid body motion only)")
    
    # Open output files
    prefix = config.output_prefix
    output_dir = config.output_dir
    
    f_energy = open(output_dir / f"{prefix}_per_atom.dat", 'w')
    f_total = open(output_dir / f"{prefix}_total.dat", 'w')
    
    try:
        # Write headers
        f_energy.write(f"# Electrostatic energy per QM atom [eV]\n")
        f_energy.write(f"# {'Fixed QM' if config.fix_qm else 'Dynamic QM'}\n")
        f_energy.write(f"# Columns: Frame, Time(fs), E_1, E_2, ..., E_N\n")
        
        f_total.write(f"# Total electrostatic energy [eV]\n")
        f_total.write(f"# {'Fixed QM' if config.fix_qm else 'Dynamic QM'}\n")
        f_total.write(f"# Columns: Frame, Time(fs), E_total\n")
        
        if config.verbose:
            print(f"\nProcessing {n_frames_selected} frames...")
        
        dt = stride  # Assume 1 fs base timestep
        
        for i in range(n_frames_selected):
            frame_id = start + i * stride
            time_val = frame_id * 1.0  # Assuming 1 fs timestep
            
            if config.verbose and (i + 1) % 1000 == 0:
                print(f"  Frame {i + 1}/{n_frames_selected}")
            
            # Extract coordinates
            coords = trajectory[i]
            qm_coords = coords[qm_indices]
            mm_coords = coords[mm_indices]
            
            # If fixing QM, use reference positions
            if config.fix_qm:
                qm_coords = compute_reference_positions(
                    qm_coords, ref_coords_centered, masses
                )
            
            # Compute energy
            energy_per_atom, total_energy = compute_electrostatic_energy(
                qm_coords, qm_charges, mm_coords, mm_charges
            )
            
            # Write results
            energy_str = "  ".join(f"{e:18.10f}" for e in energy_per_atom)
            f_energy.write(f"{frame_id:5d}  {time_val:12.3f}  {energy_str}\n")
            f_total.write(f"{frame_id:5d}  {time_val:12.3f}  {total_energy:18.10f}\n")
            
            # Flush periodically
            if (i + 1) % 100 == 0:
                f_energy.flush()
                f_total.flush()
        
        if config.verbose:
            print("Done.")
    
    finally:
        f_energy.close()
        f_total.close()
    
    print(f"\nResults saved to: {output_dir}")
    print(f"  {prefix}_per_atom.dat - Energy per QM atom")
    print(f"  {prefix}_total.dat - Total electrostatic energy")
    
    return n_frames_selected


def analyze_electrostatic_energy(
    dcd_file: str,
    top_file: str,
    qm_index_file: str,
    qm_charges_file: str,
    mm_index_file: str,
    mm_charges_file: Optional[str] = None,
    pccharges_file: Optional[str] = None,
    output_prefix: str = "elstat_energy",
    ref_xyz_file: Optional[str] = None,
    fix_qm: bool = False,
    verbose: bool = True,
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    stride: int = 1,
    charge_column: int = -1
) -> int:
    """
    Analyze electrostatic energy from trajectory.
    
    Parameters
    ----------
    dcd_file : str
        Path to DCD trajectory file
    top_file : str
        Path to topology file (PDB)
    qm_index_file : str
        Path to QM atom index file
    qm_charges_file : str
        Path to QM atom charges file
    mm_index_file : str
        Path to MM atom index file (specific MM molecule)
    mm_charges_file : str, optional
        Path to MM atom charges file (pre-extracted)
    pccharges_file : str, optional
        Path to PCcharges.dat file (ASH/DFTB output, auto-extract charges)
    output_prefix : str
        Prefix for output files
    ref_xyz_file : str, optional
        Path to reference structure (XYZ file) for fixing QM
    fix_qm : bool
        If True, fix QM molecule (rigid body motion only)
    verbose : bool
        Print progress
    start_frame : int
        First frame to process
    end_frame : int, optional
        Last frame to process
    stride : int
        Frame stride
    charge_column : int
        Column index for charge in charge files
    
    Returns
    -------
    n_frames : int
        Number of frames processed
    """
    if mm_charges_file is None and pccharges_file is None:
        raise ValueError("Either mm_charges_file or pccharges_file must be specified")
    
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    config = ElectrostaticEnergyConfig(
        trajectory=Path(dcd_file),
        topology=Path(top_file),
        qm_indices=Path(qm_index_file),
        qm_charges=Path(qm_charges_file),
        mm_indices=Path(mm_index_file),
        mm_charges=Path(mm_charges_file) if mm_charges_file else None,
        pccharges_file=Path(pccharges_file) if pccharges_file else None,
        reference=Path(ref_xyz_file) if ref_xyz_file else None,
        fix_qm=fix_qm,
        output_dir=output_dir,
        output_prefix=Path(output_prefix).name,
        start_frame=start_frame,
        end_frame=end_frame,
        stride=stride,
        verbose=verbose,
        charge_column=charge_column,
    )
    
    return run_electrostatic_energy_analysis(config)


def run_from_config_file(config_path: str) -> int:
    """Run analysis from YAML config file."""
    config = load_config(config_path)
    return run_electrostatic_energy_analysis(config)
