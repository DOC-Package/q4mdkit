"""
Run TD-DFTB calculations along a trajectory and extract excitations.

This module runs TD-DFTB (Linear Response) calculations for each frame
in a trajectory file and writes excitation energies and oscillator strengths.

Usage:
    from q4mdkit.analysis.tddftb import run_tddftb_analysis
    run_tddftb_analysis("tddftb_settings.yaml")
"""

import numpy as np
try:
    import hsd
except ModuleNotFoundError:  # pragma: no cover - optional runtime dependency
    hsd = None

try:
    import mdtraj as md
except ModuleNotFoundError:  # pragma: no cover - optional runtime dependency
    md = None
from pathlib import Path
import shutil
import yaml
import re
import subprocess
import os
from dataclasses import dataclass
from typing import List, Optional, Dict

# Physical constants
BOHR_PER_ANG = 1.0 / 0.529177


def write_xyz_file(coords_ang: np.ndarray, atom_types: List[str], filepath: Path):
    """Write coordinates to XYZ file."""
    natoms = len(atom_types)
    with open(filepath, 'w') as f:
        f.write(f"{natoms}\n")
        f.write("TD-DFTB calculation\n")
        for atom, coord in zip(atom_types, coords_ang):
            f.write(f"{atom} {coord[0]:20.10f} {coord[1]:20.10f} {coord[2]:20.10f}\n")


@dataclass
class TDDFTBConfig:
    """Configuration for TD-DFTB trajectory analysis."""
    # Input paths
    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    pccharges_template: Optional[Path]  # Optional: may not use point charges
    hsd_template: Path
    
    # Output paths
    output_dir: Path = None
    excitations_file: Path = None
    work_directory: str = "work"
    
    # Frame selection
    start_frame: int = 0
    n_frames: Optional[int] = None
    t0_fs: float = 0.0
    dt_fs: float = 4.0
    
    # DFTB+ settings
    dftb_binary: str = "dftb+"
    library_path: Optional[str] = None  # LD_LIBRARY_PATH for DFTB+ shared libs
    num_threads: Optional[int] = None
    timeout: int = 600
    
    # TD-DFTB settings
    n_excitations: int = 10


def load_tddftb_config(config_path: Path) -> TDDFTBConfig:
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
    if 'hsd_template' not in input_cfg:
        raise ValueError("Missing required input: hsd_template")
    
    traj_path = base_dir / input_cfg['trajectory']
    topology_path = base_dir / input_cfg['topology']
    qm_atoms_file = base_dir / input_cfg['qm_atoms']
    hsd_template = base_dir / input_cfg['hsd_template']
    
    # Optional point charges
    pccharges_template = None
    if 'pccharges_template' in input_cfg:
        pccharges_template = base_dir / input_cfg['pccharges_template']
    
    # Parse output paths
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output_tddftb')
    excitations_file = output_cfg.get('excitations_file', 'excitations.dat')
    work_directory = output_cfg.get('work_directory', 'work')
    
    # Parse frame settings
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    n_frames = frames_cfg.get('n_frames', None)
    t0_fs = frames_cfg.get('t0_fs', 0.0)
    dt_fs = frames_cfg.get('dt_fs', 4.0)
    
    # Parse DFTB+ settings
    dftb_cfg = data.get('dftb', {})
    dftb_binary = dftb_cfg.get('binary', 'dftb+')
    library_path = dftb_cfg.get('library_path', None)
    num_threads = dftb_cfg.get('num_threads', None)
    timeout = dftb_cfg.get('timeout', 600)
    
    # Parse TD-DFTB settings
    tddftb_cfg = data.get('tddftb', {})
    n_excitations = tddftb_cfg.get('n_excitations', 10)
    
    return TDDFTBConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        hsd_template=hsd_template,
        output_dir=output_dir,
        excitations_file=output_dir / excitations_file,
        work_directory=work_directory,
        start_frame=start_frame,
        n_frames=n_frames,
        t0_fs=t0_fs,
        dt_fs=dt_fs,
        dftb_binary=dftb_binary,
        library_path=library_path,
        num_threads=num_threads,
        timeout=timeout,
        n_excitations=n_excitations,
    )


@dataclass
class Excitation:
    """Single excitation data."""
    energy_ev: float
    oscillator_strength: float
    transition_from: Optional[int] = None
    transition_to: Optional[int] = None
    transition_weight: Optional[float] = None
    ks_energy_ev: Optional[float] = None


def parse_exc_dat(exc_dat_path: Path) -> List[Excitation]:
    """
    Parse EXC.DAT file from TD-DFTB calculation.
    
    Returns list of Excitation objects with energy and oscillator strength.
    
    EXC.DAT format:
        w [eV]       Osc.Str.         Transition         Weight      KS [eV]    Sym.
     ===============================================================================
          1.130        0.00015988       102   ->   103        1.000       1.128      S
    """
    excitations = []
    
    if not exc_dat_path.exists():
        return excitations
    
    with open(exc_dat_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header, empty lines, and separator lines
            if not line or line.startswith('#') or line.startswith('w [') or '===' in line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # Format: Energy(eV)  OscStrength  occ -> virt  weight  KS(eV)  sym
                    energy_ev = float(parts[0])
                    osc_strength = float(parts[1])
                    transition_from = None
                    transition_to = None
                    transition_weight = None
                    ks_energy_ev = None

                    if len(parts) >= 7 and parts[3] == '->':
                        transition_from = int(parts[2])
                        transition_to = int(parts[4])
                        transition_weight = float(parts[5])
                        ks_energy_ev = float(parts[6])

                    excitations.append(Excitation(
                        energy_ev=energy_ev,
                        oscillator_strength=osc_strength,
                        transition_from=transition_from,
                        transition_to=transition_to,
                        transition_weight=transition_weight,
                        ks_energy_ev=ks_energy_ev,
                    ))
                except (ValueError, IndexError):
                    continue
    
    return excitations


def parse_detailed_out(detailed_out_path: Path) -> List[Excitation]:
    """
    Parse detailed.out file from TD-DFTB calculation.
    
    Alternative parser if EXC.DAT is not available.
    """
    excitations = []
    
    if not detailed_out_path.exists():
        return excitations
    
    with open(detailed_out_path, 'r') as f:
        content = f.read()
    
    # Look for excitation energies section
    # Format varies, common patterns:
    # "Excitation Energy:" or table format
    
    # Try to find excitation table
    pattern = r'(\d+)\s+(\d+\.\d+)\s+eV\s+(\d+\.\d+)'
    matches = re.findall(pattern, content)
    
    for match in matches:
        try:
            energy_ev = float(match[1])
            osc_strength = float(match[2])
            excitations.append(Excitation(
                energy_ev=energy_ev,
                oscillator_strength=osc_strength
            ))
        except (ValueError, IndexError):
            continue
    
    return excitations


def initialize_excitation_output(output_path: Path, n_excitations: int) -> None:
    """Create the horizontal excitation energy and oscillator-strength output."""
    with open(output_path, 'w') as f:
        f.write("# TD-DFTB Excitation Data\n")
        columns = ["# Frame", "Time(fs)"]
        for state_index in range(1, n_excitations + 1):
            columns.extend([f"E{state_index}(eV)", f"f{state_index}"])
        f.write("  ".join(columns) + "\n")


def append_excitation_frame(
    output_path: Path,
    frame_data: Dict,
    n_excitations: int,
) -> None:
    """Append one frame as horizontal energy and oscillator-strength pairs."""
    fields = [f"{frame_data['frame']:6d}", f"{frame_data['time_fs']:.2f}"]
    for excitation in frame_data['excitations'][:n_excitations]:
        fields.extend([
            f"{excitation.energy_ev:.4f}",
            f"{excitation.oscillator_strength:.6f}",
        ])
    with open(output_path, 'a') as f:
        f.write("  ".join(fields) + "\n")


def prepare_tddftb_hsd(
    hsd_template: Path,
    work_dir: Path,
    n_excitations: int,
    n_pc_records: Optional[int] = None,
) -> Path:
    """
    Prepare TD-DFTB HSD input file.
    
    Ensures LinearResponse section is properly configured.
    """
    if hsd is None:
        raise ModuleNotFoundError("hsd is required to prepare TD-DFTB input files")

    with open(hsd_template, 'r') as f:
        data = hsd.load(f)
    
    # Remove Geometry section (will be added separately)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Handle point charges
    if 'Hamiltonian' in data and 'DFTB' in data['Hamiltonian']:
        dftb = data['Hamiltonian']['DFTB']
        if n_pc_records is not None:
            # Update point charges path if present
            if 'ElectricField' in dftb and 'PointCharges' in dftb['ElectricField']:
                pc = dftb['ElectricField']['PointCharges']
                if 'CoordsAndCharges' in pc and 'DirectRead' in pc['CoordsAndCharges']:
                    pc['CoordsAndCharges']['DirectRead']['File'] = 'PCcharges.dat'
                    pc['CoordsAndCharges']['DirectRead']['Records'] = n_pc_records
        else:
            # Remove ElectricField section if no point charges
            if 'ElectricField' in dftb:
                del dftb['ElectricField']
    
    # Ensure LinearResponse section exists
    if 'ExcitedState' not in data:
        data['ExcitedState'] = {}
    
    # Add/update Casida section for TD-DFTB
    if 'Casida' not in data['ExcitedState']:
        data['ExcitedState']['Casida'] = {}
    
    casida = data['ExcitedState']['Casida']
    casida['NrOfExcitations'] = n_excitations
    if 'Symmetry' not in casida:
        casida['Symmetry'] = 'Singlet'
    if 'WriteTransitions' not in casida:
        casida['WriteTransitions'] = 'Yes'
    # Diagonaliser is required for DFTB+ 24.1+
    if 'Diagonaliser' not in casida:
        casida['Diagonaliser'] = {'Stratmann': {}}
    
    # Ensure Analysis section for output
    if 'Analysis' not in data:
        data['Analysis'] = {}
    data['Analysis']['WriteEigenvectors'] = 'Yes'
    
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


def run_dftb_subprocess(
    work_dir: Path,
    dftb_binary: str = "dftb+",
    library_path: Optional[str] = None,
    num_threads: Optional[int] = None,
    timeout: int = 600
) -> bool:
    """Run DFTB+ as subprocess."""
    env = os.environ.copy()
    if library_path is not None:
        current_ld = env.get('LD_LIBRARY_PATH', '')
        env['LD_LIBRARY_PATH'] = f"{library_path}:{current_ld}" if current_ld else library_path
    if num_threads is not None:
        env['OMP_NUM_THREADS'] = str(num_threads)
    
    try:
        result = subprocess.run(
            [dftb_binary],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env
        )
        if result.returncode == 0:
            return True

        print(f"DFTB+ failed with exit code {result.returncode}")
        if result.stdout:
            stdout_tail = "\n".join(result.stdout.rstrip().splitlines()[-20:])
            if stdout_tail:
                print("DFTB+ stdout (last 20 lines):")
                print(stdout_tail)
        if result.stderr:
            stderr_tail = "\n".join(result.stderr.rstrip().splitlines()[-20:])
            if stderr_tail:
                print("DFTB+ stderr (last 20 lines):")
                print(stderr_tail)
        return False
    except subprocess.TimeoutExpired:
        print(f"DFTB+ timed out after {timeout}s")
        return False
    except Exception as e:
        print(f"DFTB+ error: {e}")
        return False


def run_tddftb_analysis(config_path: str):
    """
    Run TD-DFTB trajectory analysis.
    
    For each frame:
    1. Extract QM coordinates and optional point charges
    2. Run TD-DFTB calculation
    3. Parse excitation energies and oscillator strengths
    4. Append the requested excitations to the output file
    """
    if md is None:
        raise ModuleNotFoundError("mdtraj is required to run TD-DFTB trajectory analysis")

    from .cdftb import (
        iter_qm_coordinates,
        load_pccharges,
        load_qm_indices,
        save_pccharges,
    )

    config = load_tddftb_config(Path(config_path))
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = config.output_dir / config.work_directory
    work_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("  TD-DFTB Trajectory Analysis")
    print("="*60)
    print(f"Trajectory: {config.traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output: {config.output_dir}")
    print(f"N excitations: {config.n_excitations}")
    
    # Load QM indices
    qm_indices = load_qm_indices(config.qm_atoms_file)
    n_qm = len(qm_indices)
    print(f"QM atoms: {n_qm}")
    
    # MM indices (all non-QM atoms)
    traj_info = md.load(str(config.traj_path), top=str(config.topology_path), frame=0)
    n_atoms = traj_info.n_atoms
    mm_indices = np.array([i for i in range(n_atoms) if i not in qm_indices])
    n_mm = len(mm_indices)
    print(f"MM atoms: {n_mm}")
    
    # Load point charges template if available
    pc_template = None
    n_pc_records = None
    if config.pccharges_template and config.pccharges_template.exists():
        pc_template = load_pccharges(config.pccharges_template)
        n_pc_records = len(pc_template)
        print(f"Point charges: {n_pc_records}")
    
    # Get atom types from topology for QM atoms
    atom_types = [traj_info.topology.atom(i).element.symbol for i in qm_indices]
    print(f"QM atom types: {set(atom_types)}")
    
    # Prepare HSD template
    prepare_tddftb_hsd(
        config.hsd_template,
        work_dir,
        config.n_excitations,
        n_pc_records
    )
    
    # Storage for results
    all_excitations = []
    frame_count = 0
    initialize_excitation_output(config.excitations_file, config.n_excitations)
    
    print("\nProcessing frames...")
    
    # Determine end frame
    end_frame = None
    if config.n_frames is not None:
        end_frame = config.start_frame + config.n_frames
    
    # Iterate over trajectory
    for global_idx, time_fs, qm_coords_bohr, mm_coords_ang in iter_qm_coordinates(
        config.traj_path,
        config.topology_path,
        qm_indices,
        mm_indices,
    ):
        # Skip frames before start_frame
        if global_idx < config.start_frame:
            continue
        # Stop if we've processed enough frames
        if end_frame is not None and global_idx >= end_frame:
            break
        
        frame_count += 1
        time_current = config.t0_fs + global_idx * config.dt_fs
        
        print(f"\rFrame {global_idx} (t = {time_current:.1f} fs)...", end="", flush=True)
        
        # Save QM coordinates
        qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
        write_xyz_file(qm_coords_ang, atom_types, work_dir / "qm_coords.xyz")
        
        # Save point charges if available
        if pc_template is not None:
            # Update MM coordinates in point charges
            pc_data = pc_template.copy()
            pc_data[:, :3] = mm_coords_ang
            save_pccharges(pc_data, work_dir / "PCcharges.dat")
        
        # Run DFTB+
        success = run_dftb_subprocess(
            work_dir,
            config.dftb_binary,
            config.library_path,
            config.num_threads,
            config.timeout
        )
        
        if not success:
            print(f" [FAILED]")
            continue
        
        # Parse excitation results
        exc_dat = work_dir / "EXC.DAT"
        excitations = parse_exc_dat(exc_dat)
        
        if not excitations:
            # Try alternative parsing
            detailed_out = work_dir / "detailed.out"
            excitations = parse_detailed_out(detailed_out)
        
        if not excitations:
            print(f" [NO EXCITATIONS]")
            continue

        # Store excitations
        frame_exc_data = {
            'frame': global_idx,
            'time_fs': time_current,
            'excitations': excitations,
        }
        all_excitations.append(frame_exc_data)
        append_excitation_frame(
            config.excitations_file,
            frame_exc_data,
            config.n_excitations,
        )

        print(f" [{len(excitations)} exc]", end="")
    
    successful_frames = len(all_excitations)
    print(f"\n\nProcessed {frame_count} frames, {successful_frames} successful")

    if not all_excitations:
        print("No excitations computed!")
        return

    print(f"Saved excitations: {config.excitations_file}")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("  Summary")
    print("="*60)
    
    # Statistics on first excitation (S1)
    s1_energies = [d['excitations'][0].energy_ev for d in all_excitations if d['excitations']]
    if s1_energies:
        s1_mean = np.mean(s1_energies)
        s1_std = np.std(s1_energies)
        print(f"S1 energy: {s1_mean:.3f} ± {s1_std:.3f} eV")
    
    print("\nDone!")
