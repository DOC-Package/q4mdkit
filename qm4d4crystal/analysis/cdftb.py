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
from typing import List, Optional

# Unit conversion constants
ANG_PER_NM = 10.0  # Convert nm (MDtraj default) to Å
BOHR_PER_ANG = 1.0 / 0.529177  # Convert Å to Bohr for DFTB+ set_geometry


def parse_atom_range(atom_range: str) -> List[int]:
    """
    Parse DFTB+ style atom range (1-indexed) to Python-style [start, end) (0-indexed).
    
    Examples:
        "1:36" -> [0, 36]
        "37:72" -> [36, 72]
    """
    parts = atom_range.split(':')
    if len(parts) != 2:
        raise ValueError(f"Invalid atom range format: {atom_range}. Expected 'start:end'")
    start = int(parts[0]) - 1  # Convert 1-indexed to 0-indexed
    end = int(parts[1])        # End is exclusive in Python, so no -1
    return [start, end]


@dataclass
class FragmentConfig:
    """Configuration for a single fragment in constrained DFT."""
    name: str
    atom_range: str  # DFTB+ style, e.g., "1:36"
    charge_sum_range: List[int]  # Python-style [start, end)
    initial_charges: Optional[Path] = None  # Path to initial charges.dat for first frame


@dataclass
class CDFTBConfig:
    """Configuration for CDFTB analysis."""
    # Input paths
    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    pccharges_template: Path
    hsd_template: Path
    
    # Output paths
    output_dir: Path
    energy_file: Path
    charge_file: Path
    
    # Frame selection
    start_frame: int = 0
    n_frames: Optional[int] = None
    t0_fs: float = 0.0
    dt_fs: float = 4.0
    
    # DFTB+ settings
    dftb_library_path: str = "/home/takahashi/opt/dftb+/lib/libdftbplus.so"
    num_threads: Optional[int] = None
    timeout: int = 300
    
    # Fragment definitions
    fragments: List[FragmentConfig] = field(default_factory=list)


def load_config(config_path: Path) -> CDFTBConfig:
    """Load configuration from YAML file."""
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
        # Auto-calculate charge_sum_range from atom_range if not provided
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
    
    return CDFTBConfig(
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
        fragments=fragments
    )


# Global config object (set by run_cdftb_analysis)
_config: Optional[CDFTBConfig] = None

def load_qm_indices(path):
    """Return 0-indexed NumPy array of QM atoms (file is already 0-indexed)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"QM atom file not found: {path}")
    tokens = path.read_text().split()
    if not tokens:
        raise ValueError(f"QM atom file is empty: {path}")
    indices = np.array(sorted(int(tok) for tok in tokens), dtype=int)
    return indices


def load_mm_indices(qm_indices: np.ndarray, total_atoms: int):
    """Return 0-indexed NumPy array of MM atoms (all atoms except QM)."""
    all_atoms = set(range(total_atoms))
    mm_atoms = sorted(all_atoms - set(qm_indices))
    return np.array(mm_atoms, dtype=int)


def load_pccharges(path: Path):
    """Load point charges from file. Returns charges (4th column)."""
    if not path.exists():
        raise FileNotFoundError(f"PCcharges file not found: {path}")
    data = np.loadtxt(path)
    # columns: x, y, z, charge (in Angstrom)
    return data[:, 3]  # return only charges


def save_pccharges(path: Path, coords_ang: np.ndarray, charges: np.ndarray):
    """Save updated point charges file with new coordinates."""
    data = np.column_stack([coords_ang, charges])
    np.savetxt(path, data, fmt="%20.10f %20.10f %20.10f %10.4f")


def iter_qm_coordinates(
    traj_path: Path,
    top_path: Path,
    qm_indices: np.ndarray,
    mm_indices: np.ndarray,
    chunk_size: int = 10,
):
    """Yield (global_frame_idx, time_fs, qm_coords_bohr, mm_coords_ang) for each trajectory frame."""
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory not found: {traj_path}")
    if not top_path.exists():
        raise FileNotFoundError(f"Topology not found: {top_path}")

    frame_counter = 0
    for chunk in md.iterload(str(traj_path), top=str(top_path), chunk=chunk_size):
        times = chunk.time if chunk.time is not None else [None] * chunk.n_frames
        # QM coords: nm -> Å -> Bohr
        qm_coords_chunk = chunk.xyz[:, qm_indices, :] * ANG_PER_NM * BOHR_PER_ANG
        # MM coords: nm -> Å (PCcharges.dat uses Angstrom)
        mm_coords_chunk = chunk.xyz[:, mm_indices, :] * ANG_PER_NM
        for local_idx in range(chunk.n_frames):
            qm_coords_bohr = np.asarray(qm_coords_chunk[local_idx], dtype=np.float64)
            mm_coords_ang = np.asarray(mm_coords_chunk[local_idx], dtype=np.float64)
            time_fs = times[local_idx]
            yield frame_counter + local_idx, time_fs, qm_coords_bohr, mm_coords_ang
        frame_counter += chunk.n_frames


def extract_atom_types_from_hsd(hsd_path: Path) -> List[str]:
    """Extract atom type labels from HSD template's Geometry xyzFormat section."""
    # Read raw file to parse xyzFormat lines directly
    with open(hsd_path, 'r') as f:
        content = f.read()
    
    # Find xyzFormat block
    import re
    match = re.search(r'xyzFormat\s*\{([^}]+)\}', content, re.DOTALL)
    if not match:
        raise ValueError("HSD template must have Geometry with xyzFormat")
    
    xyz_block = match.group(1).strip()
    lines = xyz_block.split('\n')
    
    # First line: number of atoms
    natoms = int(lines[0].strip())
    # Second line: title (skip)
    # Remaining lines: atom type x y z
    atom_types = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        atom_types.append(parts[0])
    
    return atom_types


def save_qm_coords_xyz(frame_dir: Path, frame_id, time_fs, coords_ang: np.ndarray,
                       atom_types: List[str]):
    """Write QM coordinates in XYZ format (Angstrom)."""
    outfile = frame_dir / "qm_coords.xyz"
    natoms = len(atom_types)
    time_label = "NaN" if time_fs is None else f"{time_fs:.6f}"
    
    with open(outfile, 'w') as f:
        f.write(f"{natoms}\n")
        f.write(f"frame={frame_id} time_fs={time_label}\n")
        for atom, coord in zip(atom_types, coords_ang):
            f.write(f"{atom} {coord[0]:20.10f} {coord[1]:20.10f} {coord[2]:20.10f}\n")


def setup_frame_directory(frame_id: int, qm_coords_ang: np.ndarray, mm_coords_ang: np.ndarray, 
                          mm_charges: np.ndarray, time_fs, config: CDFTBConfig,
                          atom_types: List[str]) -> Path:
    """Create a directory for this frame and set up input files."""
    frame_dir = config.output_dir / f"frame_{frame_id:05d}"
    frame_dir.mkdir(parents=True, exist_ok=True)
    
    # Save QM coordinates in XYZ format
    save_qm_coords_xyz(frame_dir, frame_id, time_fs, qm_coords_ang, atom_types)
    
    # Save PCcharges.dat with updated MM coordinates
    pc_file = frame_dir / "PCcharges.dat"
    pc_data = np.column_stack([mm_coords_ang, mm_charges])
    np.savetxt(pc_file, pc_data, fmt="%20.10f %20.10f %20.10f %10.4f")
    
    # Copy and modify HSD template
    with open(config.hsd_template, 'r') as f:
        data = hsd.load(f)
    
    # Remove Geometry from data (will be added as text later)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Update paths to be relative to frame directory
    # PCcharges path
    data['Hamiltonian']['DFTB']['ElectricField']['PointCharges']['CoordsAndCharges']['DirectRead']['File'] = 'PCcharges.dat'
    
    # Save modified HSD
    hsd_file = frame_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        # Write Geometry section with file include directive first
        f.write('Geometry {\n')
        f.write('  xyzFormat {\n')
        f.write('    <<< "qm_coords.xyz"\n')
        f.write('  }\n')
        f.write('}\n\n')
        # Then write the rest of the HSD
        hsd.dump(data, f)
    
    return frame_dir


def setup_fragment_directory(frame_dir: Path, fragment_name: str, constrained_atoms: str, 
                              hsd_template: Path) -> Path:
    """
    Create a fragment subdirectory with constraint on specific atoms.
    
    Parameters:
    -----------
    frame_dir : Path
        Parent frame directory
    fragment_name : str
        'fragment1' or 'fragment2'
    constrained_atoms : str
        Atom range for constraint, e.g., '1:36' or '37:72'
    hsd_template : Path
        Path to HSD template file
    """
    frag_dir = frame_dir / fragment_name
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
            # Handle both single constraint (dict) and multiple constraints (list)
            if isinstance(mulliken, list):
                for constraint in mulliken:
                    constraint['Atoms'] = constrained_atoms
            else:
                mulliken['Atoms'] = constrained_atoms
    
    # Modify InitialSpins.Atoms to guide spin localization on the constrained fragment
    if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
        spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
        if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
            atom_spin = spin_pol['Colinear']['InitialSpins']['AtomSpin']
            atom_spin['Atoms'] = constrained_atoms
    
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


def run_dftb_calculation(qm_coords_bohr, frame_dir, write_hs, result_queue, dftb_library_path, num_threads):
    """Run DFTB+ calculation in a separate process within frame directory.
    
    Parameters:
    -----------
    qm_coords_bohr : np.ndarray
        QM atom coordinates in Bohr units for set_geometry()
    frame_dir : Path
        Directory containing dftb_in.hsd
    write_hs : bool
        Whether to write Hamiltonian/Overlap matrices
    result_queue : mp.Queue
        Queue to return results
    dftb_library_path : str
        Path to DFTB+ library
    num_threads : int or None
        Number of OpenMP threads (None = use system default)
    """
    try:
        # Set OpenMP threads if specified
        if num_threads is not None:
            os.environ['OMP_NUM_THREADS'] = str(num_threads)
        
        # Change to frame directory
        original_dir = os.getcwd()
        os.chdir(frame_dir)
        
        # Modify HSD for this calculation
        with open("dftb_in.hsd", 'r') as f:
            data = hsd.load(f)
        
        if write_hs:
            data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'Yes'
            data['Options']['WriteHS'] = 'Yes'
            data['Options']['ReadChargesAsText'] = 'Yes'
            # Disable unnecessary output for WriteHS calculation
            data['Analysis']['WriteBandOut'] = 'No'
            data['Options']['WriteResultsTag'] = 'No'
            # Remove InitialSpins when using ReadInitialCharges
            if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
                spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
                if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
                    del spin_pol['Colinear']['InitialSpins']
        else:
            # SCC calculation: output charges as text for subsequent WriteHS calculation
            data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'No'
            data['Options']['WriteHS'] = 'No'
            data['Options']['WriteChargesAsText'] = 'Yes'  # Output charges.dat for WriteHS
            if 'ReadChargesAsText' in data['Options']:
                del data['Options']['ReadChargesAsText']
        
        with open("dftb_in.hsd", 'w') as f:
            hsd.dump(data, f)
        
        cdftb = dftbplus.DftbPlus(
            libpath=dftb_library_path,
            hsdpath="dftb_in.hsd",
            logfile="dftb.log",
        )
        # Set geometry triggers the SCC calculation
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


def run_dftb_in_subprocess(qm_coords_bohr, frame_dir, write_hs=False,
                           dftb_library_path="/home/takahashi/opt/dftb+/lib/libdftbplus.so",
                           num_threads=None, timeout=300):
    """Run DFTB+ in a subprocess that can be killed if it hangs or crashes.
    
    Parameters:
    -----------
    qm_coords_bohr : np.ndarray
        QM atom coordinates in Bohr units for set_geometry()
    frame_dir : Path
        Directory containing dftb_in.hsd
    write_hs : bool
        Whether to write Hamiltonian/Overlap matrices
    dftb_library_path : str
        Path to DFTB+ library
    num_threads : int or None
        Number of OpenMP threads (None = use system default)
    timeout : int
        Timeout in seconds
    """
    # Set OpenMP threads in parent process before forking
    # This ensures child process inherits the correct environment
    old_omp_threads = os.environ.get('OMP_NUM_THREADS')
    if num_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(num_threads)
    
    result_queue = mp.Queue()
    proc = mp.Process(target=run_dftb_calculation, args=(qm_coords_bohr, frame_dir, write_hs, result_queue, dftb_library_path, num_threads))
    proc.start()
    proc.join(timeout=timeout)
    
    # Restore original OMP_NUM_THREADS
    if old_omp_threads is not None:
        os.environ['OMP_NUM_THREADS'] = old_omp_threads
    elif num_threads is not None:
        del os.environ['OMP_NUM_THREADS']
    
    if proc.is_alive():
        # Timeout - kill the process
        proc.terminate()
        proc.join()
        return None, None, "Timeout"
    
    if proc.exitcode != 0:
        # Process crashed (e.g., Fortran ERROR STOP)
        return None, None, f"Process crashed with exit code {proc.exitcode}"
    
    if result_queue.empty():
        return None, None, "No result returned"
    
    status, energy_or_error, mcharge = result_queue.get()
    if status == "success":
        return energy_or_error, mcharge, None
    else:
        return None, None, energy_or_error