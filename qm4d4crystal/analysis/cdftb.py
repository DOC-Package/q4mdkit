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
        frag = FragmentConfig(
            name=frag_data['name'],
            atom_range=atom_range,
            charge_sum_range=charge_sum_range
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

def load_qm_indices(path: Path):
    """Return 0-indexed NumPy array of QM atoms (file is already 0-indexed)."""
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
            data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'No'
            data['Options']['WriteHS'] = 'No'
            if 'ReadChargesAsText' in data['Options']:
                del data['Options']['ReadChargesAsText']
        
        with open("dftb_in.hsd", 'w') as f:
            hsd.dump(data, f)
        
        cdftb = dftbplus.DftbPlus(
            libpath=dftb_library_path,
            hsdpath="dftb_in.hsd",
            logfile="dftb.log",
        )
        # Set geometry triggers the SCF calculation
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
    result_queue = mp.Queue()
    proc = mp.Process(target=run_dftb_calculation, args=(qm_coords_bohr, frame_dir, write_hs, result_queue, dftb_library_path, num_threads))
    proc.start()
    proc.join(timeout=timeout)
    
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


def run_cdftb_analysis(config_path: Path):
    """
    Run constrained DFT-B analysis using settings from YAML config file.
    
    Parameters:
    -----------
    config_path : Path
        Path to YAML configuration file
    """
    # Load configuration
    config = load_config(config_path)
    
    print("=" * 70)
    print("CDFTB Analysis")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Trajectory: {config.traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output directory: {config.output_dir}")
    print("=" * 70)
    
    qm_indices = load_qm_indices(config.qm_atoms_file)
    print(f"Loaded {len(qm_indices)} QM atom indices from {config.qm_atoms_file}")
    print(f"QM atoms (0-indexed): {qm_indices}")

    # Load total atoms from topology to get MM indices
    traj_info = md.load(str(config.traj_path), top=str(config.topology_path), frame=0)
    total_atoms = traj_info.n_atoms
    mm_indices = load_mm_indices(qm_indices, total_atoms)
    print(f"MM atoms: {len(mm_indices)}")

    # Load original charges from PCcharges.dat
    mm_charges = load_pccharges(config.pccharges_template)
    print(f"Loaded {len(mm_charges)} point charges from {config.pccharges_template}")
    
    # Extract atom types from HSD template for xyz file generation
    atom_types = extract_atom_types_from_hsd(config.hsd_template)
    print(f"Extracted {len(atom_types)} atom types from HSD template")
    
    # Print fragment information
    print(f"\nFragments ({len(config.fragments)}):")
    for frag in config.fragments:
        print(f"  - {frag.name}: atoms {frag.atom_range}, charge sum range {frag.charge_sum_range}")
    print()

    # Create output directory and open energy file
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build header for energy file
    energy_header_parts = ["# Frame", "Time[fs]"]
    for frag in config.fragments:
        energy_header_parts.append(f"Energy_{frag.name}[a.u.]")
    energy_header = "  ".join(energy_header_parts) + "\n"
    
    # Build header for charge file
    # Format: Qn|calcm = charge on fragment n when constraint applied to fragment m
    charge_header_lines = []
    charge_header_lines.append("# Mulliken charges from constrained DFT calculations")
    charge_header_lines.append("# Qn|calcm = charge on fragment n when constraint applied to fragment m")
    charge_header_parts = [f"{'# Frame':>7s}", f"{'Time[fs]':>12s}"]
    for i, frag in enumerate(config.fragments, 1):
        for j, frag2 in enumerate(config.fragments, 1):
            charge_header_parts.append(f"{'Q' + str(j) + '|calc' + str(i):>12s}")
    charge_header_lines.append("  ".join(charge_header_parts))
    charge_header = "\n".join(charge_header_lines) + "\n"
    
    energy_file = open(config.energy_file, "w")
    energy_file.write(energy_header)
    charge_file = open(config.charge_file, "w")
    charge_file.write(charge_header)

    try:
        # Track frame count for n_frames limit
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
            
            # Convert QM coords to Angstrom for saving
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Setup frame directory with input files
            frame_dir = setup_frame_directory(frame_id, qm_coords_ang, mm_coords_ang, mm_charges, time_fs, config, atom_types)
            
            # Calculate time using user-defined t0 and dt
            time_val = config.t0_fs + frame_id * config.dt_fs
            time_str = f"t = {time_val:.3f} fs"
            
            print(f"Frame {frame_id:05d} ({time_str})")
            
            # Process each fragment from config
            energies = []
            charges = []
            
            for frag in config.fragments:
                frag_dir = setup_fragment_directory(
                    frame_dir, frag.name, frag.atom_range, config.hsd_template
                )
                
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
                else:
                    start_idx, end_idx = frag.charge_sum_range
                    frag_charge = float(np.sum(mcharge[start_idx:end_idx]))
                    # Calculate charges for all fragments
                    all_frag_charges = []
                    for f in config.fragments:
                        s, e = f.charge_sum_range
                        all_frag_charges.append(float(np.sum(mcharge[s:e])))
                    charges_str = ", ".join([f"Q_{f.name}={q:+.6f}" for f, q in zip(config.fragments, all_frag_charges)])
                    print(f"  {frag.name}: E = {energy:.10f} a.u., {charges_str}")
                    
                    # Rename detailed.out and dftb.log to preserve them before WriteHS calculation
                    detailed_out = frag_dir / "detailed.out"
                    detailed_out_saved = frag_dir / "detailed_scc.out"
                    if detailed_out.exists():
                        shutil.move(detailed_out, detailed_out_saved)
                    
                    dftb_log = frag_dir / "dftb.log"
                    dftb_log_saved = frag_dir / "dftb_scc.log"
                    if dftb_log.exists():
                        shutil.move(dftb_log, dftb_log_saved)
                    
                    # WriteHS calculation
                    energy_hs, _, error_hs = run_dftb_in_subprocess(
                        qm_coords_bohr, frag_dir, write_hs=True,
                        dftb_library_path=config.dftb_library_path,
                        num_threads=config.num_threads,
                        timeout=config.timeout
                    )
                    if error_hs:
                        print(f"    WriteHS: {error_hs}")
                    else:
                        print(f"    WriteHS: OK")
                
                energies.append(energy)
                charges.append(all_frag_charges)
            
            # Save to energy file
            energy_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for e in energies:
                energy_parts.append(f"{e:18.10f}")
            energy_file.write("  ".join(energy_parts) + "\n")
            energy_file.flush()
            
            # Save to charge file
            charge_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for qs in charges:
                for q in qs:
                    charge_parts.append(f"{q:+12.6f}")
            charge_file.write("  ".join(charge_parts) + "\n")
            charge_file.flush()
            
            processed_count += 1
            print()

    finally:
        energy_file.close()
        charge_file.close()
        print(f"Energies saved to {config.energy_file}")
        print(f"Charges saved to {config.charge_file}")


def main():
    """Main entry point for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run constrained DFT-B analysis on MD trajectory"
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
    
    run_cdftb_analysis(config_path)
    return 0


if __name__ == "__main__":
    exit(main())