"""
Run standard DFTB calculations along a trajectory.

This module runs closed-shell, non-constrained DFTB calculations
for each frame in a trajectory file. It outputs energies.

Usage:
    from q4mdkit.analysis.dftb import run_dftb_analysis
    run_dftb_analysis("dftb_settings.yaml")
"""

import numpy as np
import hsd
import mdtraj as md
from pathlib import Path
import shutil
import yaml
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

# Import utility functions from cdftb module
from .cdftb import (
    ANG_PER_NM, BOHR_PER_ANG,
    load_qm_indices,
    load_mm_indices,
    load_pccharges,
    save_pccharges,
    iter_qm_coordinates,
    extract_atom_types_from_hsd,
    save_qm_coords_xyz,
)
from .electrostatic_energy import (
    COULOMB_CONSTANT_EV_ANGSTROM,
    read_index_file,
    read_structure_file,
)
import dftbplus
import multiprocessing as mp
import os


@dataclass
class DFTBConfig:
    """Configuration for standard DFTB trajectory analysis."""
    # Input paths
    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    pccharges_template: Path
    hsd_template: Path
    
    # Optional: separate QM trajectory file (DCD format)
    qm_traj_path: Optional[Path] = None
    qm_topology_path: Optional[Path] = None
    
    # Output paths
    output_dir: Path = None
    energy_file: Path = None
    work_directory: str = "work"
    
    # Frame selection
    start_frame: int = 0
    n_frames: Optional[int] = None
    t0_fs: float = 0.0
    dt_fs: float = 4.0
    
    # DFTB+ settings
    dftb_library_path: str = "/home/takahashi/opt/dftb+/lib/libdftbplus.so"
    num_threads: Optional[int] = None
    timeout: int = 300
    
    # SCC settings
    use_previous_charges: bool = False
    
    # Output settings
    output_detailed_energy: bool = False  # Output Electronic, Repulsive, Total, PointCharges
    output_mulliken_charges: bool = False  # Output Mulliken charges for each atom
    
    # Fixed geometry mode: only vary point charges, keep QM geometry fixed
    fixed_geometry: bool = False
    fixed_geometry_file: Optional[Path] = None  # Path to fixed QM structure (xyz/gen)
    
    # Electrostatic energy calculation with specific MM molecule
    output_elstat_energy: bool = False  # Compute E_elstat = Σ q_a * Σ Q_b / r_ab
    mm_molecule_indices_file: Optional[Path] = None  # Indices of target MM molecule (system indices)



def load_dftb_config(config_path: Path) -> DFTBConfig:
    """Load configuration from YAML file."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse input paths (relative to config file)
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
    
    # Optional: separate QM trajectory file (DCD format)
    qm_traj_path = None
    qm_topology_path = None
    if 'qm_trajectory' in input_cfg:
        qm_traj_path = base_dir / input_cfg['qm_trajectory']
        if 'qm_topology' in input_cfg:
            qm_topology_path = base_dir / input_cfg['qm_topology']
        else:
            raise ValueError("qm_topology is required when qm_trajectory is specified")
    
    # Parse output paths
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output')
    energy_file = output_cfg.get('energy_file', 'energies.dat')
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
    
    # Parse SCC settings
    scc_cfg = data.get('scc', {})
    use_previous_charges = scc_cfg.get('use_previous_charges', False)
    
    # Output detail settings
    output_detailed_energy = output_cfg.get('detailed_energy', False)
    output_mulliken_charges = output_cfg.get('mulliken_charges', False)
    
    # Fixed geometry mode settings
    fixed_cfg = data.get('fixed_geometry', {})
    fixed_geometry = fixed_cfg.get('enabled', False)
    fixed_geometry_file = None
    if fixed_geometry and 'structure_file' in fixed_cfg:
        fixed_geometry_file = base_dir / fixed_cfg['structure_file']
    
    # Electrostatic energy settings
    output_elstat_energy = output_cfg.get('elstat_energy', False)
    mm_molecule_indices_file = None
    if output_elstat_energy:
        if 'mm_molecule_indices' not in input_cfg:
            raise ValueError("mm_molecule_indices required when output.elstat_energy=true")
        mm_molecule_indices_file = base_dir / input_cfg['mm_molecule_indices']
    
    return DFTBConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        hsd_template=hsd_template,
        qm_traj_path=qm_traj_path,
        qm_topology_path=qm_topology_path,
        output_dir=output_dir,
        energy_file=output_dir / energy_file,
        work_directory=work_directory,
        start_frame=start_frame,
        n_frames=n_frames,
        t0_fs=t0_fs,
        dt_fs=dt_fs,
        dftb_library_path=dftb_library_path,
        num_threads=num_threads,
        timeout=timeout,
        use_previous_charges=use_previous_charges,
        output_detailed_energy=output_detailed_energy,
        output_mulliken_charges=output_mulliken_charges,
        fixed_geometry=fixed_geometry,
        fixed_geometry_file=fixed_geometry_file,
        output_elstat_energy=output_elstat_energy,
        mm_molecule_indices_file=mm_molecule_indices_file,
    )


def iter_qm_coordinates_separate(
    qm_traj_path: Path,
    qm_top_path: Path,
    traj_path: Path,
    top_path: Path,
    mm_indices: np.ndarray,
    chunk_size: int = 10,
):
    """
    Yield (global_frame_idx, time_fs, qm_coords_bohr, mm_coords_ang) for each trajectory frame,
    reading QM coordinates from a separate trajectory file (DCD/XTC etc.).
    
    Parameters
    ----------
    qm_traj_path : Path
        Path to QM trajectory file (DCD/XTC etc.) containing only QM atoms.
    qm_top_path : Path
        Path to topology file for QM trajectory.
    traj_path : Path
        Path to main trajectory file (DCD/XTC etc.) for MM coordinates.
    top_path : Path
        Path to topology file for main trajectory.
    mm_indices : np.ndarray
        Indices of MM atoms in the main trajectory.
    chunk_size : int
        Number of frames to load at a time.
    """
    if not qm_traj_path.exists():
        raise FileNotFoundError(f"QM trajectory not found: {qm_traj_path}")
    if not qm_top_path.exists():
        raise FileNotFoundError(f"QM topology not found: {qm_top_path}")
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory not found: {traj_path}")
    if not top_path.exists():
        raise FileNotFoundError(f"Topology not found: {top_path}")
    
    # Load QM trajectory info
    qm_traj_info = md.load(str(qm_traj_path), top=str(qm_top_path), frame=0)
    print(f"QM trajectory atoms: {qm_traj_info.n_atoms}")
    
    # Create iterators for both trajectories
    qm_iter = md.iterload(str(qm_traj_path), top=str(qm_top_path), chunk=chunk_size)
    mm_iter = md.iterload(str(traj_path), top=str(top_path), chunk=chunk_size)
    
    frame_counter = 0
    
    for qm_chunk, mm_chunk in zip(qm_iter, mm_iter):
        # Check frame counts match
        if qm_chunk.n_frames != mm_chunk.n_frames:
            print(f"Warning: Frame count mismatch at chunk starting at frame {frame_counter}")
            min_frames = min(qm_chunk.n_frames, mm_chunk.n_frames)
        else:
            min_frames = qm_chunk.n_frames
        
        # QM coords: nm -> Å -> Bohr (all atoms in QM trajectory)
        qm_coords_chunk = qm_chunk.xyz[:min_frames, :, :] * ANG_PER_NM * BOHR_PER_ANG
        # MM coords: nm -> Å (PCcharges.dat uses Angstrom)
        mm_coords_chunk = mm_chunk.xyz[:min_frames, mm_indices, :] * ANG_PER_NM
        times = mm_chunk.time[:min_frames] if mm_chunk.time is not None else [None] * min_frames
        
        for local_idx in range(min_frames):
            qm_coords_bohr = np.asarray(qm_coords_chunk[local_idx], dtype=np.float64)
            mm_coords_ang = np.asarray(mm_coords_chunk[local_idx], dtype=np.float64)
            time_fs = times[local_idx]
            yield frame_counter + local_idx, time_fs, qm_coords_bohr, mm_coords_ang
        
        frame_counter += min_frames


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
    Set up work directory with current frame's coordinates.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Save QM coordinates
    save_qm_coords_xyz(work_dir, frame_id, time_fs, qm_coords_ang, atom_types)
    
    # Save PCcharges.dat with updated MM coordinates
    pc_file = work_dir / "PCcharges.dat"
    save_pccharges(pc_file, mm_coords_ang, mm_charges)


def setup_dftb_hsd(
    work_dir: Path,
    hsd_template: Path,
    n_pc_records: Optional[int] = None,
    use_previous_charges: bool = False,
) -> Path:
    """
    Set up HSD file by updating coordinates and PCcharges paths.
    Uses hsd library like cdftbci.
    """
    # Load template HSD
    data = hsd.load(str(hsd_template))
    
    # Remove Geometry section (will be added separately)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Update PCcharges.dat path and Records count
    if 'Hamiltonian' in data and 'DFTB' in data['Hamiltonian']:
        dftb = data['Hamiltonian']['DFTB']
        if 'ElectricField' in dftb and 'PointCharges' in dftb['ElectricField']:
            pc = dftb['ElectricField']['PointCharges']
            if 'CoordsAndCharges' in pc and 'DirectRead' in pc['CoordsAndCharges']:
                pc['CoordsAndCharges']['DirectRead']['File'] = 'PCcharges.dat'
                if n_pc_records is not None:
                    pc['CoordsAndCharges']['DirectRead']['Records'] = n_pc_records
    
    # Set ReadInitialCharges if using previous frame's charges
    if use_previous_charges and 'Hamiltonian' in data and 'DFTB' in data['Hamiltonian']:
        data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'Yes'
        # Need ReadChargesAsText to read text format charges.dat
        if 'Options' not in data:
            data['Options'] = {}
        data['Options']['ReadChargesAsText'] = 'Yes'
        # Remove InitialSpins when using ReadInitialCharges (they conflict)
        if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
            spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
            if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
                del spin_pol['Colinear']['InitialSpins']
    
    # Save modified HSD
    hsd_file = work_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        # Write Geometry section with file include directive first
        f.write('Geometry {\n')
        f.write('  xyzFormat {\n')
        f.write('    <<< "qm_coords.xyz"\n')
        f.write('  }\n')
        f.write('}\n\n')
        # Then write the rest of the HSD
        hsd.dump(data, f)
    
    return hsd_file


def run_dftb_calculation(qm_coords_bohr, work_dir, result_queue, dftb_library_path, num_threads):
    """Run DFTB+ calculation in a separate process."""
    try:
        if num_threads is not None:
            os.environ['OMP_NUM_THREADS'] = str(num_threads)
        
        original_dir = os.getcwd()
        os.chdir(work_dir)
        
        cdftb = dftbplus.DftbPlus(
            libpath=dftb_library_path,
            hsdpath="dftb_in.hsd",
            logfile="dftb.log",
        )
        cdftb.set_geometry(qm_coords_bohr)
        energy = cdftb.get_energy()
        mcharge = cdftb.get_gross_charges()
        cdftb.close()
        
        os.chdir(original_dir)
        result_queue.put(("success", energy, mcharge))
    except Exception as e:
        try:
            os.chdir(original_dir)
        except:
            pass
        result_queue.put(("error", str(e), None))


def run_dftb_in_subprocess(qm_coords_bohr, work_dir,
                           dftb_library_path="/home/takahashi/opt/dftb+/lib/libdftbplus.so",
                           num_threads=None, timeout=300):
    """Run DFTB+ in subprocess without modifying HSD."""
    old_omp_threads = os.environ.get('OMP_NUM_THREADS')
    if num_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(num_threads)
    
    result_queue = mp.Queue()
    proc = mp.Process(target=run_dftb_calculation, 
                     args=(qm_coords_bohr, work_dir, result_queue, dftb_library_path, num_threads))
    proc.start()
    proc.join(timeout=timeout)
    
    if old_omp_threads is not None:
        os.environ['OMP_NUM_THREADS'] = old_omp_threads
    elif num_threads is not None:
        del os.environ['OMP_NUM_THREADS']
    
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return None, None, "Timeout"
    
    if proc.exitcode != 0:
        return None, None, f"Process crashed with exit code {proc.exitcode}"
    
    if result_queue.empty():
        return None, None, "No result returned"
    
    status, energy_or_error, mcharge = result_queue.get()
    if status == "success":
        return energy_or_error, mcharge, None
    else:
        return None, None, energy_or_error


def read_detailed_energy(detailed_path):
    """
    Read energy components from DFTB+ detailed.out file.
    
    Args:
        detailed_path: Path to detailed.out file.
    
    Returns:
        dict: {'electronic': float, 'repulsive': float, 'total': float, 'point_charges': float}
              or None if not found.
    """
    p = Path(detailed_path)
    if not p.exists():
        return None
    txt = p.read_text(errors="ignore")
    
    # Parse energy lines from detailed.out (Hartree only)
    pat_elec = re.compile(r"Total Electronic energy:\s+([-+]?\d+\.\d+)\s+H")
    pat_rep = re.compile(r"Repulsive energy:\s+([-+]?\d+\.\d+)\s+H")
    pat_pc = re.compile(r"Energy point charges:\s+([-+]?\d+\.\d+)\s+H")
    pat_total = re.compile(r"Total energy:\s+([-+]?\d+\.\d+)\s+H")
    
    m_elec = pat_elec.search(txt)
    m_rep = pat_rep.search(txt)
    m_pc = pat_pc.search(txt)
    m_total = pat_total.search(txt)
    
    if not (m_elec and m_rep and m_total):
        return None
    
    return {
        'electronic': float(m_elec.group(1)),
        'repulsive': float(m_rep.group(1)),
        'total': float(m_total.group(1)),
        'point_charges': float(m_pc.group(1)) if m_pc else 0.0,
    }


def load_fixed_geometry(filepath: Path) -> np.ndarray:
    """
    Load fixed QM geometry from XYZ or GEN file.
    
    Parameters
    ----------
    filepath : Path
        Path to structure file (.xyz or .gen format)
    
    Returns
    -------
    coords_bohr : np.ndarray
        Atomic coordinates in Bohr (N, 3)
    """
    coords_bohr, _ = read_structure_file(str(filepath), output_unit="bohr")
    return coords_bohr


def compute_mm_molecule_elstat_energy(
    qm_coords_ang: np.ndarray,
    qm_charges: np.ndarray,
    mm_coords_ang: np.ndarray,
    mm_charges: np.ndarray,
    mm_molecule_pc_indices: np.ndarray,
) -> float:
    """
    Compute electrostatic energy between QM atoms and a specific MM molecule.
    
    E = Σ_a q_a * Σ_b Q_b / r_ab
    
    Parameters
    ----------
    qm_coords_ang : np.ndarray
        QM atom coordinates in Angstrom (n_qm, 3)
    qm_charges : np.ndarray
        Mulliken charges of QM atoms (n_qm,)
    mm_coords_ang : np.ndarray
        All MM atom coordinates in Angstrom (n_mm, 3)
    mm_charges : np.ndarray
        All MM atom charges (n_mm,)
    mm_molecule_pc_indices : np.ndarray
        Indices of target MM molecule in PCcharges order (subset of 0..n_mm-1)
    
    Returns
    -------
    energy_ev : float
        Electrostatic interaction energy in eV
    """
    # Extract target MM molecule coords and charges
    mol_coords = mm_coords_ang[mm_molecule_pc_indices]  # (n_mol, 3)
    mol_charges = mm_charges[mm_molecule_pc_indices]    # (n_mol,)
    
    # Compute pairwise distances (n_qm, n_mol)
    diff = qm_coords_ang[:, np.newaxis, :] - mol_coords[np.newaxis, :, :]  # (n_qm, n_mol, 3)
    distances = np.linalg.norm(diff, axis=2)  # (n_qm, n_mol)
    
    # E = Σ_a q_a * Σ_b Q_b / r_ab
    # Note: we sum over all qm-mm pairs
    energy_ev = COULOMB_CONSTANT_EV_ANGSTROM * np.sum(
        qm_charges[:, np.newaxis] * mol_charges[np.newaxis, :] / distances
    )
    
    return energy_ev


def map_system_indices_to_pccharges(
    system_indices: np.ndarray,
    qm_indices: np.ndarray,
) -> np.ndarray:
    """
    Map system-wide atom indices to PCcharges.dat row indices.
    
    PCcharges.dat contains only MM atoms, so we need to subtract
    the number of QM atoms that come before each MM atom in the system.
    
    Parameters
    ----------
    system_indices : np.ndarray
        System-wide atom indices (0-based)
    qm_indices : np.ndarray
        QM atom indices in system
    
    Returns
    -------
    pc_indices : np.ndarray
        Corresponding row indices in PCcharges.dat
    """
    qm_set = set(qm_indices)
    pc_indices = []
    
    for sys_idx in system_indices:
        if sys_idx in qm_set:
            raise ValueError(f"System index {sys_idx} is a QM atom, not MM")
        # Count how many QM atoms are before this system index
        n_qm_before = np.sum(qm_indices < sys_idx)
        pc_idx = sys_idx - n_qm_before
        pc_indices.append(pc_idx)
    
    return np.array(pc_indices, dtype=np.int32)


def read_nearby_molecules_file(
    filename: str,
    qm_indices: np.ndarray,
) -> List[Tuple[int, np.ndarray]]:
    """
    Read nearby_molecules file format.
    
    Format: each line is one MM molecule
        mol_id  atom_idx1  atom_idx2  ...  atom_idxN
    
    Parameters
    ----------
    filename : str
        Path to nearby_molecules file
    qm_indices : np.ndarray
        QM atom indices (for mapping to PCcharges)
    
    Returns
    -------
    molecules : list of (mol_id, pc_indices)
        List of tuples containing molecule ID and PCcharges indices
    """
    molecules = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            mol_id = int(parts[0])
            atom_indices = np.array([int(x) for x in parts[1:]], dtype=np.int32)
            
            # Map to PCcharges indices
            pc_indices = map_system_indices_to_pccharges(atom_indices, qm_indices)
            molecules.append((mol_id, pc_indices))
    
    return molecules


def run_dftb_analysis(config_path: str):
    """
    Run standard DFTB analysis on trajectory frames.
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file.
    """
    # Load configuration
    config = load_dftb_config(config_path)
    
    print("=" * 70)
    print("DFTB Trajectory Analysis (Non-constrained)")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Trajectory: {config.traj_path}")
    if config.qm_traj_path is not None:
        print(f"QM Trajectory: {config.qm_traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output directory: {config.output_dir}")
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
    
    if config.use_previous_charges:
        print(f"\nSCC settings:")
        print(f"  Use previous charges: Yes")
    
    # Fixed geometry mode
    fixed_qm_coords_bohr = None
    if config.fixed_geometry:
        print(f"\nFixed geometry mode: ENABLED")
        if config.fixed_geometry_file is not None:
            print(f"  Structure file: {config.fixed_geometry_file}")
            fixed_qm_coords_bohr = load_fixed_geometry(config.fixed_geometry_file)
        else:
            print(f"  Using first frame as fixed structure")
        print(f"  Only point charges will vary along trajectory")
    
    # Load MM molecule indices for electrostatic energy calculation
    mm_molecules = None
    if config.output_elstat_energy:
        mm_molecules = read_nearby_molecules_file(
            str(config.mm_molecule_indices_file), qm_indices
        )
        print(f"\nElectrostatic energy calculation: ENABLED")
        print(f"  Number of MM molecules: {len(mm_molecules)}")
        for mol_id, pc_indices in mm_molecules:
            print(f"    Molecule {mol_id}: {len(pc_indices)} atoms")
    print()
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up work directory path
    work_dir = config.output_dir / config.work_directory
    
    # Open output file
    energy_file = open(config.energy_file, "w")
    energy_file.write("# DFTB Trajectory Analysis (Closed-shell, Non-constrained)\n")
    if config.output_detailed_energy:
        energy_file.write("# Frame  Time(fs)        Electronic(a.u.)      Repulsive(a.u.)       Total(a.u.)           PointCharges(a.u.)\n")
    else:
        energy_file.write("# Frame  Time(fs)        Energy(a.u.)\n")
    
    # Open Mulliken charges file if requested
    mulliken_file = None
    if config.output_mulliken_charges:
        mulliken_path = config.output_dir / "mulliken_charges.dat"
        mulliken_file = open(mulliken_path, "w")
        mulliken_file.write("# Mulliken charges for each frame\n")
        mulliken_file.write(f"# Atom types: {' '.join(atom_types)}\n")
        mulliken_file.write(f"# Columns: Frame, Time(fs), Q_1, Q_2, ..., Q_N\n")
        mulliken_file.write(f"# Units: elementary charge (e)\n")
    
    # Open electrostatic energy file if requested
    elstat_file = None
    if config.output_elstat_energy:
        elstat_path = config.output_dir / "elstat_energy.dat"
        elstat_file = open(elstat_path, "w")
        elstat_file.write("# Electrostatic energy between QM and MM molecules\n")
        elstat_file.write("# E = Σ_a q_a * Σ_b Q_b / r_ab (using Mulliken charges)\n")
        mol_ids = [mol_id for mol_id, _ in mm_molecules]
        elstat_file.write(f"# Molecule IDs: {' '.join(str(m) for m in mol_ids)}\n")
        header = "# Frame  Time(fs)        " + "  ".join(f"E_mol{mol_id:03d}(eV)" for mol_id in mol_ids) + "\n"
        elstat_file.write(header)
    
    try:
        processed_count = 0
        
        # Choose iterator based on whether separate QM trajectory is provided
        if config.qm_traj_path is not None:
            coord_iterator = iter_qm_coordinates_separate(
                config.qm_traj_path, config.qm_topology_path,
                config.traj_path, config.topology_path, mm_indices
            )
        else:
            coord_iterator = iter_qm_coordinates(
                config.traj_path, config.topology_path, qm_indices, mm_indices
            )
        
        for frame_id, time_fs, qm_coords_bohr, mm_coords_ang in coord_iterator:
            # Skip frames before start_frame
            if frame_id < config.start_frame:
                continue
            
            # Stop if we've processed enough frames
            if config.n_frames is not None and processed_count >= config.n_frames:
                break
            
            # Fixed geometry mode: use fixed QM structure
            if config.fixed_geometry:
                if fixed_qm_coords_bohr is None:
                    # Use first frame as fixed structure
                    fixed_qm_coords_bohr = qm_coords_bohr.copy()
                    print(f"  Using frame {frame_id} as fixed QM structure")
                qm_coords_bohr = fixed_qm_coords_bohr
            
            # Convert QM coords to Angstrom
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Calculate time using user-defined t0 and dt
            time_val = config.t0_fs + frame_id * config.dt_fs
            
            print(f"Frame {frame_id:05d} (t = {time_val:.3f} fs)")
            
            # Set up work directory with current coordinates
            setup_work_directory(
                work_dir, qm_coords_ang, mm_coords_ang, mm_charges,
                frame_id, time_fs, atom_types
            )
            
            # Set up HSD file
            use_prev_charges = config.use_previous_charges and processed_count > 0
            setup_dftb_hsd(
                work_dir,
                config.hsd_template,
                n_pc_records=len(mm_charges),
                use_previous_charges=use_prev_charges,
            )
            
            # Run DFTB+ calculation
            energy, mcharge, error = run_dftb_in_subprocess(
                qm_coords_bohr, work_dir,
                dftb_library_path=config.dftb_library_path,
                num_threads=config.num_threads,
                timeout=config.timeout
            )
            
            if error:
                print(f"  ERROR: {error}")
                energy = float("nan")
                mcharge = None
                detailed_energies = None
            else:
                print(f"  E = {energy:.10f} a.u.")
                # Read detailed energies if requested
                if config.output_detailed_energy:
                    detailed_energies = read_detailed_energy(work_dir / "detailed.out")
                else:
                    detailed_energies = None
            
            # Save results
            if config.output_detailed_energy and detailed_energies is not None:
                energy_file.write(f"{frame_id:5d}  {time_val:12.3f}  {detailed_energies['electronic']:18.10f}  {detailed_energies['repulsive']:18.10f}  {detailed_energies['total']:18.10f}  {detailed_energies['point_charges']:18.10f}\n")
            else:
                energy_file.write(f"{frame_id:5d}  {time_val:12.3f}  {energy:18.10f}\n")
            energy_file.flush()
            
            # Save Mulliken charges if requested
            if config.output_mulliken_charges and mulliken_file is not None:
                if mcharge is not None:
                    charges_str = "  ".join(f"{q:12.8f}" for q in mcharge)
                    mulliken_file.write(f"{frame_id:5d}  {time_val:12.3f}  {charges_str}\n")
                else:
                    # Write NaN for error frames
                    nan_str = "  ".join("         nan" for _ in atom_types)
                    mulliken_file.write(f"{frame_id:5d}  {time_val:12.3f}  {nan_str}\n")
                mulliken_file.flush()
            
            # Save electrostatic energy if requested
            if config.output_elstat_energy and elstat_file is not None:
                if mcharge is not None:
                    energies = []
                    for mol_id, pc_indices in mm_molecules:
                        e = compute_mm_molecule_elstat_energy(
                            qm_coords_ang, mcharge, mm_coords_ang, mm_charges,
                            pc_indices
                        )
                        energies.append(e)
                    energies_str = "  ".join(f"{e:18.10f}" for e in energies)
                    elstat_file.write(f"{frame_id:5d}  {time_val:12.3f}  {energies_str}\n")
                    total_e = sum(energies)
                    print(f"  E_elstat = {total_e:.6f} eV (total of {len(mm_molecules)} molecules)")
                else:
                    nan_str = "  ".join("               nan" for _ in mm_molecules)
                    elstat_file.write(f"{frame_id:5d}  {time_val:12.3f}  {nan_str}\n")
                elstat_file.flush()
            
            processed_count += 1
    
    finally:
        energy_file.close()
        if mulliken_file is not None:
            mulliken_file.close()
        if elstat_file is not None:
            elstat_file.close()
    
    print("\n" + "=" * 70)
    print(f"Energies saved to {config.energy_file}")
    if config.output_mulliken_charges:
        print(f"Mulliken charges saved to {config.output_dir / 'mulliken_charges.dat'}")
    if config.output_elstat_energy:
        print(f"Electrostatic energy saved to {config.output_dir / 'elstat_energy.dat'}")
    print("=" * 70)


def main():
    """Command-line entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Run DFTB analysis on trajectory")
    parser.add_argument("config", help="Path to YAML configuration file")
    args = parser.parse_args()
    
    run_dftb_analysis(args.config)


if __name__ == "__main__":
    main()
