import numpy as np
import dftbplus
import hsd
import mdtraj as md
from pathlib import Path
import multiprocessing as mp
import os
import shutil

BASE_DIR = Path(__file__).resolve().parent
TRAJ_PATH = BASE_DIR / "input" / "prod.dcd"
TOPOLOGY_PATH = BASE_DIR.parent / "pentacene.gro"
QM_ATOMS_FILE = BASE_DIR.parent / "qmatoms"
PCCHARGES_TEMPLATE = BASE_DIR / "input" / "PCcharges.dat"
HSD_TEMPLATE = BASE_DIR / "input" / "dftb_in.hsd"
OUTPUT_DIR = BASE_DIR / "output"
ENERGY_FILE = OUTPUT_DIR / "energies.dat"
ANG_PER_NM = 10.0  # Convert nm (MDtraj default) to Å
BOHR_PER_ANG = 1.0 / 0.529177  # Convert Å to Bohr for DFTB+ set_geometry

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


def save_qm_coords(frame_dir: Path, frame_id, time_fs, coords_ang):
    """Write QM coordinates for debugging in plain text."""
    outfile = frame_dir / "qm_coords.txt"
    time_label = "NaN" if time_fs is None else f"{time_fs:.6f}"
    header = f"frame={frame_id}, time_fs={time_label}, coordinates in Angstrom"
    np.savetxt(outfile, coords_ang, header=header, fmt="%20.10f")


def setup_frame_directory(frame_id: int, qm_coords_ang: np.ndarray, mm_coords_ang: np.ndarray, 
                          mm_charges: np.ndarray, time_fs) -> Path:
    """Create a directory for this frame and set up input files."""
    frame_dir = OUTPUT_DIR / f"frame_{frame_id:05d}"
    frame_dir.mkdir(parents=True, exist_ok=True)
    
    # Create input subdirectory
    input_dir = frame_dir / "input"
    input_dir.mkdir(exist_ok=True)
    
    # Save QM coordinates (for reference)
    save_qm_coords(frame_dir, frame_id, time_fs, qm_coords_ang)
    
    # Save PCcharges.dat with updated MM coordinates
    pc_file = input_dir / "PCcharges.dat"
    pc_data = np.column_stack([mm_coords_ang, mm_charges])
    np.savetxt(pc_file, pc_data, fmt="%20.10f %20.10f %20.10f %10.4f")
    
    # Copy and modify HSD template
    with open(HSD_TEMPLATE, 'r') as f:
        data = hsd.load(f)
    
    # Update paths to be relative to frame directory
    # PCcharges path
    data['Hamiltonian']['DFTB']['ElectricField']['PointCharges']['CoordsAndCharges']['DirectRead']['File'] = 'input/PCcharges.dat'
    
    # Save modified HSD
    hsd_file = frame_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        hsd.dump(data, f)
    
    return frame_dir


def setup_fragment_directory(frame_dir: Path, fragment_name: str, constrained_atoms: str, 
                              mm_coords_ang: np.ndarray, mm_charges: np.ndarray) -> Path:
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
    mm_coords_ang : np.ndarray
        MM coordinates in Angstrom
    mm_charges : np.ndarray
        MM point charges
    """
    frag_dir = frame_dir / fragment_name
    frag_dir.mkdir(parents=True, exist_ok=True)
    
    # Create input subdirectory
    input_dir = frag_dir / "input"
    input_dir.mkdir(exist_ok=True)
    
    # Save PCcharges.dat
    pc_file = input_dir / "PCcharges.dat"
    pc_data = np.column_stack([mm_coords_ang, mm_charges])
    np.savetxt(pc_file, pc_data, fmt="%20.10f %20.10f %20.10f %10.4f")
    
    # Copy and modify HSD template
    with open(HSD_TEMPLATE, 'r') as f:
        data = hsd.load(f)
    
    # Update PCcharges path
    data['Hamiltonian']['DFTB']['ElectricField']['PointCharges']['CoordsAndCharges']['DirectRead']['File'] = 'input/PCcharges.dat'
    
    # Modify ElectronicConstraints to constrain specific fragment
    if 'ElectronicConstraints' in data['Hamiltonian']['DFTB']:
        constraints = data['Hamiltonian']['DFTB']['ElectronicConstraints']['Constraints']
        if 'MullikenPopulation' in constraints:
            constraints['MullikenPopulation']['Atoms'] = constrained_atoms
    
    # Save modified HSD
    hsd_file = frag_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        hsd.dump(data, f)
    
    return frag_dir


def run_dftb_calculation(qm_coords_bohr, frame_dir, write_hs, result_queue):
    """Run DFTB+ calculation in a separate process within frame directory."""
    try:
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
        else:
            data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'No'
            data['Options']['WriteHS'] = 'No'
            if 'ReadChargesAsText' in data['Options']:
                del data['Options']['ReadChargesAsText']
        
        with open("dftb_in.hsd", 'w') as f:
            hsd.dump(data, f)
        
        cdftb = dftbplus.DftbPlus(
            libpath="/home/takahashi/opt/dftb+/lib/libdftbplus.so",
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


def run_dftb_in_subprocess(qm_coords_bohr, frame_dir, write_hs=False, timeout=300):
    """Run DFTB+ in a subprocess that can be killed if it hangs or crashes."""
    result_queue = mp.Queue()
    proc = mp.Process(target=run_dftb_calculation, args=(qm_coords_bohr, frame_dir, write_hs, result_queue))
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


def main():
    # =========================================================================
    # User parameters - modify these as needed
    # =========================================================================
    start_frame = 0       # First frame to process (0-indexed)
    n_frames = None       # Number of frames to process (None = all remaining frames)
    t0_fs = 0.0           # Initial time (fs) for output
    dt_fs = 4.0           # Time step (fs) between frames
    # =========================================================================
    
    qm_indices = load_qm_indices(QM_ATOMS_FILE)
    print(f"Loaded {len(qm_indices)} QM atom indices from {QM_ATOMS_FILE}")
    print(f"QM atoms (0-indexed): {qm_indices}")

    # Load total atoms from topology to get MM indices
    traj_info = md.load(str(TRAJ_PATH), top=str(TOPOLOGY_PATH), frame=0)
    total_atoms = traj_info.n_atoms
    mm_indices = load_mm_indices(qm_indices, total_atoms)
    print(f"MM atoms: {len(mm_indices)}")

    # Load original charges from PCcharges.dat
    mm_charges = load_pccharges(PCCHARGES_TEMPLATE)
    print(f"Loaded {len(mm_charges)} point charges from {PCCHARGES_TEMPLATE}")

    # Create output directory and open energy file
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    energy_file = open(ENERGY_FILE, "w")
    energy_file.write("# Frame  Time[fs]  Energy_frag1[a.u.]  Energy_frag2[a.u.]  Charge_frag1[e]  Charge_frag2[e]\n")

    try:
        # Track frame count for n_frames limit
        processed_count = 0
        
        for frame_id, time_fs, qm_coords_bohr, mm_coords_ang in iter_qm_coordinates(
            TRAJ_PATH, TOPOLOGY_PATH, qm_indices, mm_indices
        ):
            # Skip frames before start_frame
            if frame_id < start_frame:
                continue
            
            # Stop if we've processed enough frames
            if n_frames is not None and processed_count >= n_frames:
                break
            
            # Convert QM coords to Angstrom for saving
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Setup frame directory with input files
            frame_dir = setup_frame_directory(frame_id, qm_coords_ang, mm_coords_ang, mm_charges, time_fs)
            
            # Calculate time using user-defined t0 and dt
            time_val = t0_fs + frame_id * dt_fs
            time_str = f"t = {time_val:.3f} fs"
            
            print(f"Frame {frame_id:05d} ({time_str})")
            
            # =====================================================================
            # Fragment 1: Constrain atoms 1:36 to have charge +1
            # =====================================================================
            frag1_dir = setup_fragment_directory(frame_dir, "fragment1", "1:36", mm_coords_ang, mm_charges)
            
            energy1, mcharge1, error1 = run_dftb_in_subprocess(qm_coords_bohr, frag1_dir, write_hs=False)
            
            if error1:
                print(f"  Fragment1: ERROR - {error1}")
                frag1_charge = float("nan")
                energy1 = float("nan")
            else:
                frag1_charge = float(np.sum(mcharge1[:36]))
                print(f"  Fragment1: E = {energy1:.10f} a.u., Q1 = {frag1_charge:+.6f}")
                
                # Rename detailed.out and dftb.log to preserve them before WriteHS calculation
                detailed_out = frag1_dir / "detailed.out"
                detailed_out_saved = frag1_dir / "detailed_scc.out"
                if detailed_out.exists():
                    shutil.move(detailed_out, detailed_out_saved)
                
                dftb_log = frag1_dir / "dftb.log"
                dftb_log_saved = frag1_dir / "dftb_scc.log"
                if dftb_log.exists():
                    shutil.move(dftb_log, dftb_log_saved)
                
                # WriteHS calculation
                energy1_hs, _, error1_hs = run_dftb_in_subprocess(qm_coords_bohr, frag1_dir, write_hs=True)
                if error1_hs:
                    print(f"    WriteHS: {error1_hs}")
                else:
                    print(f"    WriteHS: OK")
            
            # =====================================================================
            # Fragment 2: Constrain atoms 37:72 to have charge +1
            # =====================================================================
            frag2_dir = setup_fragment_directory(frame_dir, "fragment2", "37:72", mm_coords_ang, mm_charges)
            
            energy2, mcharge2, error2 = run_dftb_in_subprocess(qm_coords_bohr, frag2_dir, write_hs=False)
            
            if error2:
                print(f"  Fragment2: ERROR - {error2}")
                frag2_charge = float("nan")
                energy2 = float("nan")
            else:
                frag2_charge = float(np.sum(mcharge2[36:72]))
                print(f"  Fragment2: E = {energy2:.10f} a.u., Q2 = {frag2_charge:+.6f}")
                
                # Rename detailed.out and dftb.log to preserve them before WriteHS calculation
                detailed_out = frag2_dir / "detailed.out"
                detailed_out_saved = frag2_dir / "detailed_scc.out"
                if detailed_out.exists():
                    shutil.move(detailed_out, detailed_out_saved)
                
                dftb_log = frag2_dir / "dftb.log"
                dftb_log_saved = frag2_dir / "dftb_scc.log"
                if dftb_log.exists():
                    shutil.move(dftb_log, dftb_log_saved)
                
                # WriteHS calculation
                energy2_hs, _, error2_hs = run_dftb_in_subprocess(qm_coords_bohr, frag2_dir, write_hs=True)
                if error2_hs:
                    print(f"    WriteHS: {error2_hs}")
                else:
                    print(f"    WriteHS: OK")
            
            # Save to energy file
            energy_file.write(f"{frame_id:5d}  {time_val:12.3f}  {energy1:18.10f}  {energy2:18.10f}  {frag1_charge:+12.6f}  {frag2_charge:+12.6f}\n")
            energy_file.flush()
            
            processed_count += 1
            print()

    finally:
        energy_file.close()
        print(f"Energies saved to {ENERGY_FILE}")

if __name__ == "__main__":
    main()