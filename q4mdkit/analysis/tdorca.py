"""
Run ORCA TD-DFT calculations along a trajectory for absorption spectra.

This module runs ORCA TD-DFT calculations for each frame in a trajectory file,
computes absorption spectra with Gaussian broadening, and averages the spectra
across frames.

Usage:
    from q4mdkit.analysis.tdorca import run_tdorca_analysis
    run_tdorca_analysis("tdorca_settings.yaml")
"""

import numpy as np
import mdtraj as md
from pathlib import Path
import shutil
import yaml
import re
import subprocess
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

# Unit conversion constants
ANG_PER_NM = 10.0  # Convert nm (MDtraj default) to Å
BOHR_PER_ANG = 1.0 / 0.529177  # Convert Å to Bohr

# Physical constants
EV_PER_HARTREE = 27.211386245988
NM_PER_EV = 1239.8419  # hc in eV*nm


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
    mm_indices: np.ndarray = None,
    chunk_size: int = 10,
):
    """Yield (global_frame_idx, time_fs, qm_coords_ang, mm_coords_ang) for each trajectory frame."""
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory not found: {traj_path}")
    if not top_path.exists():
        raise FileNotFoundError(f"Topology not found: {top_path}")

    frame_counter = 0
    for chunk in md.iterload(str(traj_path), top=str(top_path), chunk=chunk_size):
        times = chunk.time if chunk.time is not None else [None] * chunk.n_frames
        # QM coords: nm -> Å
        qm_coords_chunk = chunk.xyz[:, qm_indices, :] * ANG_PER_NM
        # MM coords: nm -> Å (if provided)
        mm_coords_chunk = None
        if mm_indices is not None:
            mm_coords_chunk = chunk.xyz[:, mm_indices, :] * ANG_PER_NM
        for local_idx in range(chunk.n_frames):
            qm_coords_ang = np.asarray(qm_coords_chunk[local_idx], dtype=np.float64)
            mm_coords_ang = None
            if mm_coords_chunk is not None:
                mm_coords_ang = np.asarray(mm_coords_chunk[local_idx], dtype=np.float64)
            time_fs = times[local_idx]
            yield frame_counter + local_idx, time_fs, qm_coords_ang, mm_coords_ang
        frame_counter += chunk.n_frames


def write_xyz_file(coords_ang: np.ndarray, atom_types: List[str], filepath: Path):
    """Write coordinates to XYZ file."""
    natoms = len(atom_types)
    with open(filepath, 'w') as f:
        f.write(f"{natoms}\n")
        f.write("ORCA TD-DFT calculation\n")
        for atom, coord in zip(atom_types, coords_ang):
            f.write(f"{atom} {coord[0]:20.10f} {coord[1]:20.10f} {coord[2]:20.10f}\n")


@dataclass
class TDORCAConfig:
    """Configuration for ORCA TD-DFT trajectory analysis."""
    # Input paths
    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    pccharges_template: Optional[Path]  # Optional: point charges for QM/MM
    orca_template: Path  # ORCA input template
    
    # Output paths
    output_dir: Path = None
    spectra_file: Path = None
    excitations_file: Path = None
    work_directory: str = "work"
    
    # Frame selection
    start_frame: int = 0
    n_frames: Optional[int] = None
    t0_fs: float = 0.0
    dt_fs: float = 4.0
    
    # ORCA settings
    orca_binary: str = "orca"
    num_procs: int = 1
    num_threads: Optional[int] = None
    memory_mb: int = 4000
    timeout: int = 3600  # 1 hour default for DFT
    
    # TD-DFT settings
    n_roots: int = 10  # Number of excited states
    broadening_ev: float = 0.1  # Gaussian broadening width in eV
    energy_range_ev: Tuple[float, float] = (0.0, 6.0)  # Energy range for spectrum
    n_energy_points: int = 1000  # Number of points in spectrum
    
    # Spectrum output
    output_individual_spectra: bool = False  # Save spectrum for each frame


def load_tdorca_config(config_path: Path) -> TDORCAConfig:
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
    if 'orca_template' not in input_cfg:
        raise ValueError("Missing required input: orca_template")
    
    traj_path = (base_dir / input_cfg['trajectory']).resolve()
    topology_path = (base_dir / input_cfg['topology']).resolve()
    qm_atoms_file = (base_dir / input_cfg['qm_atoms']).resolve()
    orca_template = (base_dir / input_cfg['orca_template']).resolve()
    
    pccharges_template = None
    if 'pccharges_template' in input_cfg:
        pccharges_template = (base_dir / input_cfg['pccharges_template']).resolve()
    
    # Parse output settings
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output_tdorca')
    spectra_file = output_cfg.get('spectra_file', 'absorption_spectrum.dat')
    excitations_file = output_cfg.get('excitations_file', 'excitations.dat')
    work_directory = output_cfg.get('work_directory', 'work')
    output_individual_spectra = output_cfg.get('individual_spectra', False)
    
    # Parse frame selection
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    n_frames = frames_cfg.get('n_frames', None)
    t0_fs = frames_cfg.get('t0_fs', 0.0)
    dt_fs = frames_cfg.get('dt_fs', 4.0)
    
    # Parse ORCA settings
    orca_cfg = data.get('orca', {})
    orca_binary = orca_cfg.get('binary', 'orca')
    num_procs = orca_cfg.get('num_procs', 1)
    num_threads = orca_cfg.get('num_threads', None)
    memory_mb = orca_cfg.get('memory_mb', 4000)
    timeout = orca_cfg.get('timeout', 3600)
    
    # Parse TD-DFT settings
    tddft_cfg = data.get('tddft', {})
    n_roots = tddft_cfg.get('n_roots', 10)
    broadening_ev = tddft_cfg.get('broadening_ev', 0.1)
    energy_range = tddft_cfg.get('energy_range_ev', [0.0, 6.0])
    n_energy_points = tddft_cfg.get('n_energy_points', 1000)
    
    return TDORCAConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        orca_template=orca_template,
        output_dir=output_dir,
        spectra_file=output_dir / spectra_file,
        excitations_file=output_dir / excitations_file,
        work_directory=work_directory,
        start_frame=start_frame,
        n_frames=n_frames,
        t0_fs=t0_fs,
        dt_fs=dt_fs,
        orca_binary=orca_binary,
        num_procs=num_procs,
        num_threads=num_threads,
        memory_mb=memory_mb,
        timeout=timeout,
        n_roots=n_roots,
        broadening_ev=broadening_ev,
        energy_range_ev=tuple(energy_range),
        n_energy_points=n_energy_points,
        output_individual_spectra=output_individual_spectra,
    )


@dataclass
class Excitation:
    """Single excitation data."""
    energy_ev: float
    oscillator_strength: float
    wavelength_nm: float = None
    
    def __post_init__(self):
        if self.wavelength_nm is None and self.energy_ev > 0:
            self.wavelength_nm = NM_PER_EV / self.energy_ev


def parse_orca_output(output_path: Path) -> List[Excitation]:
    """
    Parse ORCA output file for TD-DFT excitation data.
    
    Looks for the ABSORPTION SPECTRUM section:
    -------------------------------------------------------------------------------
                         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
    -------------------------------------------------------------------------------
    State   Energy    Wavelength  fosc         T2        TX        TY        TZ
            (cm-1)      (nm)                 (au**2)    (au)      (au)      (au)
    -------------------------------------------------------------------------------
       1   19532.6    511.96   0.000000012   0.00000   0.00000   0.00000   0.00000
       2   22945.5    435.83   0.000123456   0.00123   0.00100   0.00050   0.00020
    """
    excitations = []
    
    if not output_path.exists():
        return excitations
    
    with open(output_path, 'r') as f:
        content = f.read()
    
    # Find absorption spectrum section
    # Pattern for ORCA 5.x output
    spectrum_pattern = r'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS'
    
    if spectrum_pattern not in content:
        # Try alternative format for older ORCA versions
        return parse_orca_output_alternative(content)
    
    # Find the spectrum section
    spectrum_start = content.find(spectrum_pattern)
    if spectrum_start == -1:
        return excitations
    
    # Get the spectrum section (until next section or end)
    spectrum_section = content[spectrum_start:spectrum_start+10000]
    
    # Parse each excitation line
    # Format: State   Energy(cm-1)    Wavelength(nm)  fosc  ...
    pattern = r'^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.Ee+-]+)'
    
    for line in spectrum_section.split('\n'):
        match = re.match(pattern, line)
        if match:
            try:
                state = int(match.group(1))
                energy_cm = float(match.group(2))
                wavelength_nm = float(match.group(3))
                fosc = float(match.group(4))
                
                # Convert cm-1 to eV
                energy_ev = energy_cm / 8065.544  # 1 eV = 8065.544 cm-1
                
                excitations.append(Excitation(
                    energy_ev=energy_ev,
                    oscillator_strength=fosc,
                    wavelength_nm=wavelength_nm
                ))
            except (ValueError, IndexError):
                continue
    
    return excitations


def parse_orca_output_alternative(content: str) -> List[Excitation]:
    """
    Alternative parser for ORCA TD-DFT output.
    
    Looks for the excited states section:
    -------------------------
    TD-DFT EXCITED STATES
    -------------------------
    ...
    STATE  1:  E=   0.088989 au      2.422 eV    19532.6 cm**-1 <S**2>=   0.000
    """
    excitations = []
    
    # Pattern for excited states in ORCA output
    # STATE  1:  E=   0.088989 au      2.422 eV    19532.6 cm**-1
    pattern = r'STATE\s+(\d+):\s+E=\s*([\d.]+)\s+au\s+([\d.]+)\s+eV'
    
    states = re.findall(pattern, content)
    
    # Now find oscillator strengths
    # Look for line like: f(osc)=  0.00000012
    fosc_pattern = r'STATE\s+(\d+).*?f\(osc\)=\s*([\d.Ee+-]+)'
    fosc_dict = {}
    for match in re.finditer(fosc_pattern, content, re.DOTALL):
        state = int(match.group(1))
        fosc = float(match.group(2))
        fosc_dict[state] = fosc
    
    # Also try finding fosc in the spectrum table
    spectrum_pattern = r'^\s*(\d+)\s+[\d.]+\s+[\d.]+\s+([\d.Ee+-]+)'
    for line in content.split('\n'):
        match = re.match(spectrum_pattern, line)
        if match:
            state = int(match.group(1))
            fosc = float(match.group(2))
            fosc_dict[state] = fosc
    
    for match in states:
        try:
            state = int(match[0])
            energy_ev = float(match[2])
            fosc = fosc_dict.get(state, 0.0)
            
            excitations.append(Excitation(
                energy_ev=energy_ev,
                oscillator_strength=fosc
            ))
        except (ValueError, IndexError):
            continue
    
    return excitations


def compute_absorption_spectrum(
    excitations: List[Excitation],
    energy_range: Tuple[float, float],
    n_points: int,
    broadening: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute absorption spectrum with Gaussian broadening.
    
    Parameters
    ----------
    excitations : List[Excitation]
        List of excitation data.
    energy_range : Tuple[float, float]
        Energy range (min, max) in eV.
    n_points : int
        Number of points in spectrum.
    broadening : float
        Gaussian broadening width (sigma) in eV.
    
    Returns
    -------
    energies : np.ndarray
        Energy grid in eV.
    spectrum : np.ndarray
        Absorption intensity (oscillator strength weighted).
    """
    energies = np.linspace(energy_range[0], energy_range[1], n_points)
    spectrum = np.zeros(n_points)
    
    for exc in excitations:
        # Gaussian broadening: f(E) = f_osc * exp(-(E - E_exc)^2 / (2*sigma^2))
        gaussian = exc.oscillator_strength * np.exp(
            -(energies - exc.energy_ev)**2 / (2 * broadening**2)
        )
        spectrum += gaussian
    
    return energies, spectrum


def prepare_orca_input(
    template_path: Path,
    xyz_file: Path,
    output_path: Path,
    n_roots: int,
    num_procs: int = 1,
    memory_mb: int = 4000,
    pc_file: Optional[Path] = None
) -> Path:
    """
    Prepare ORCA input file from template.
    
    Template should contain placeholders:
    - {XYZ_FILE} or %coords ... end block will be replaced
    - {NROOTS} for number of roots
    - {NPROCS} for number of processors
    - {MEMORY} for memory in MB
    - {PC_FILE} for point charges file (optional)
    
    Or just use a complete template that we'll append geometry to.
    """
    with open(template_path, 'r') as f:
        template = f.read()
    
    # Replace placeholders
    orca_input = template
    orca_input = orca_input.replace('{NROOTS}', str(n_roots))
    orca_input = orca_input.replace('{NPROCS}', str(num_procs))
    orca_input = orca_input.replace('{MEMORY}', str(memory_mb))
    
    # Handle point charges
    if pc_file is not None and '{PC_FILE}' in orca_input:
        orca_input = orca_input.replace('{PC_FILE}', str(pc_file))
    elif pc_file is not None:
        # Add point charges section if not in template
        pc_section = f'\n%pointcharges "{pc_file}"\n'
        # Insert before geometry
        if '* xyz' in orca_input.lower() or '*xyz' in orca_input.lower():
            orca_input = pc_section + orca_input
    
    # Handle geometry
    if '{XYZ_FILE}' in orca_input:
        # Read XYZ coordinates
        with open(xyz_file, 'r') as f:
            xyz_lines = f.readlines()[2:]  # Skip header
        coords_block = ''.join(xyz_lines)
        orca_input = orca_input.replace('{XYZ_FILE}', coords_block.strip())
    elif '* xyzfile' in orca_input.lower() or '*xyzfile' in orca_input.lower():
        # Replace xyzfile reference
        orca_input = re.sub(
            r'\*\s*xyzfile\s+\d+\s+\d+\s+\S+',
            f'* xyzfile 0 1 {xyz_file}',
            orca_input,
            flags=re.IGNORECASE
        )
    else:
        # Append XYZ geometry at the end
        with open(xyz_file, 'r') as f:
            xyz_lines = f.readlines()[2:]  # Skip header
        coords_block = ''.join(xyz_lines)
        
        # Check if there's already a geometry block
        if '* xyz' not in orca_input.lower() and '*xyz' not in orca_input.lower():
            orca_input += f'\n* xyz 0 1\n{coords_block}*\n'
    
    # Write output
    with open(output_path, 'w') as f:
        f.write(orca_input)
    
    return output_path


def run_orca_subprocess(
    work_dir: Path,
    input_file: Path,
    orca_binary: str,
    num_threads: Optional[int] = None,
    timeout: int = 3600
) -> bool:
    """
    Run ORCA calculation.
    
    Returns True if calculation completed successfully.
    """
    env = os.environ.copy()
    
    if num_threads is not None:
        env['OMP_NUM_THREADS'] = str(num_threads)
    
    cmd = [orca_binary, str(input_file.name)]
    output_file = work_dir / input_file.with_suffix('.out').name
    
    try:
        # Run ORCA and capture output to file
        with open(output_file, 'w') as f_out:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                env=env,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                timeout=timeout
            )
        
        # Check for successful completion
        if output_file.exists():
            with open(output_file, 'r') as f:
                content = f.read()
                if 'ORCA TERMINATED NORMALLY' in content:
                    return True
                if '****ORCA TERMINATED NORMALLY****' in content:
                    return True
        
        return False
        
    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        return False


def run_tdorca_analysis(config_path: str):
    """
    Run ORCA TD-DFT analysis along a trajectory.
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file.
    """
    print(f"Running ORCA TD-DFT analysis with config: {config_path}")
    
    config = load_tdorca_config(Path(config_path))
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = config.output_dir / config.work_directory
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Load trajectory
    traj = md.load(str(config.traj_path), top=str(config.topology_path))
    n_frames_total = traj.n_frames
    
    # Load QM indices
    qm_indices = load_qm_indices(config.qm_atoms_file)
    n_qm = len(qm_indices)
    n_mm = traj.n_atoms - n_qm
    
    # Get atom types for QM region
    topology = traj.topology
    atom_types = [topology.atom(i).element.symbol for i in qm_indices]
    
    # Load point charges template if provided
    pc_template = None
    if config.pccharges_template is not None:
        pc_template = load_pccharges(config.pccharges_template)
    
    # Print configuration
    print("=" * 60)
    print("  ORCA TD-DFT Trajectory Analysis")
    print("=" * 60)
    print(f"Trajectory: {config.traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output: {config.output_dir}")
    print(f"N roots: {config.n_roots}")
    print(f"Broadening: {config.broadening_ev} eV")
    print(f"Energy range: {config.energy_range_ev[0]} - {config.energy_range_ev[1]} eV")
    print(f"QM atoms: {n_qm}")
    print(f"MM atoms: {n_mm}")
    print(f"QM atom types: {set(atom_types)}")
    print()
    
    # Process frames
    print("Processing frames...")
    print()
    
    all_excitations = []
    all_spectra = []
    energies = np.linspace(
        config.energy_range_ev[0],
        config.energy_range_ev[1],
        config.n_energy_points
    )
    wavelengths = np.where(energies > 0, NM_PER_EV / energies, np.inf)
    
    # Determine frames to process
    end_frame = config.start_frame + (config.n_frames if config.n_frames else n_frames_total)
    end_frame = min(end_frame, n_frames_total)
    
    frame_count = 0
    for frame_idx in range(config.start_frame, end_frame):
        global_idx = frame_idx
        time_current = config.t0_fs + frame_idx * config.dt_fs
        
        print(f"\rFrame {global_idx} (t = {time_current:.1f} fs)...", end="", flush=True)
        frame_count += 1
        
        # Get coordinates
        frame = traj[frame_idx]
        all_coords_nm = frame.xyz[0]  # Shape: (n_atoms, 3)
        all_coords_ang = all_coords_nm * 10.0  # nm to Angstrom
        
        # QM coordinates
        qm_coords_ang = all_coords_ang[qm_indices]
        
        # MM coordinates (everything except QM)
        mm_indices = [i for i in range(traj.n_atoms) if i not in qm_indices]
        mm_coords_ang = all_coords_ang[mm_indices]
        
        # Write QM geometry
        xyz_file = work_dir / "qm_coords.xyz"
        write_xyz_file(qm_coords_ang, atom_types, xyz_file)
        
        # Write point charges if available
        pc_file = None
        if pc_template is not None:
            pc_data = pc_template.copy()
            pc_data[:, :3] = mm_coords_ang
            pc_file = work_dir / "pointcharges.pc"
            # ORCA point charge format: charge x y z
            with open(pc_file, 'w') as f:
                f.write(f"{len(pc_data)}\n")
                for row in pc_data:
                    # ORCA format: q x y z
                    f.write(f"{row[3]:12.6f} {row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}\n")
        
        # Prepare ORCA input
        input_file = work_dir / "orca.inp"
        prepare_orca_input(
            config.orca_template,
            xyz_file,
            input_file,
            config.n_roots,
            config.num_procs,
            config.memory_mb,
            pc_file
        )
        
        # Run ORCA
        success = run_orca_subprocess(
            work_dir,
            input_file,
            config.orca_binary,
            config.num_threads,
            config.timeout
        )
        
        if not success:
            print(f" [FAILED]")
            continue
        
        # Parse excitation results
        output_file = work_dir / "orca.out"
        excitations = parse_orca_output(output_file)
        
        if not excitations:
            print(f" [NO EXCITATIONS]")
            continue
        
        # Store excitations
        frame_exc_data = {
            'frame': global_idx,
            'time_fs': time_current,
            'excitations': excitations
        }
        all_excitations.append(frame_exc_data)
        
        # Compute spectrum for this frame
        _, spectrum = compute_absorption_spectrum(
            excitations,
            config.energy_range_ev,
            config.n_energy_points,
            config.broadening_ev
        )
        all_spectra.append(spectrum)
        
        # Save individual spectrum if requested
        if config.output_individual_spectra:
            spec_file = config.output_dir / f"spectrum_frame{global_idx:06d}.dat"
            np.savetxt(
                spec_file,
                np.column_stack([energies, wavelengths, spectrum]),
                header="Energy(eV)  Wavelength(nm)  Intensity",
                fmt="%.6f"
            )
        
        print(f" [{len(excitations)} exc]")
    
    print(f"\nProcessed {frame_count} frames, {len(all_spectra)} successful")
    
    if not all_spectra:
        print("No spectra computed!")
        return
    
    # Compute average spectrum
    spectra_array = np.array(all_spectra)
    avg_spectrum = np.mean(spectra_array, axis=0)
    std_spectrum = np.std(spectra_array, axis=0)
    
    # Save average spectrum
    with open(config.spectra_file, 'w') as f:
        f.write(f"# ORCA TD-DFT Average Absorption Spectrum\n")
        f.write(f"# N_frames: {len(all_spectra)}\n")
        f.write(f"# Broadening: {config.broadening_ev} eV\n")
        f.write(f"# Energy(eV)  Wavelength(nm)  Intensity  StdDev\n")
        for i in range(len(energies)):
            f.write(f"{energies[i]:.6f}  {wavelengths[i]:.2f}  {avg_spectrum[i]:.8f}  {std_spectrum[i]:.8f}\n")
    
    print(f"Saved average spectrum: {config.spectra_file}")
    
    # Save excitation data
    with open(config.excitations_file, 'w') as f:
        f.write("# ORCA TD-DFT Excitation Data\n")
        f.write("# Frame  Time(fs)  State  Energy(eV)  Wavelength(nm)  OscStrength\n")
        for frame_data in all_excitations:
            for i, exc in enumerate(frame_data['excitations']):
                f.write(f"{frame_data['frame']:6d}  {frame_data['time_fs']:.2f}  ")
                f.write(f"{i+1:3d}  {exc.energy_ev:.4f}  {exc.wavelength_nm:.2f}  {exc.oscillator_strength:.6f}\n")
    
    print(f"Saved excitations: {config.excitations_file}")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("  Summary")
    print("=" * 60)
    
    # Find peak
    peak_idx = np.argmax(avg_spectrum)
    peak_energy = energies[peak_idx]
    peak_wavelength = wavelengths[peak_idx]
    print(f"Peak absorption: {peak_energy:.3f} eV ({peak_wavelength:.1f} nm)")
    
    # S1 statistics
    s1_energies = [fd['excitations'][0].energy_ev for fd in all_excitations if fd['excitations']]
    if s1_energies:
        s1_mean = np.mean(s1_energies)
        s1_std = np.std(s1_energies)
        s1_wavelength = NM_PER_EV / s1_mean
        print(f"S1 energy: {s1_mean:.3f} ± {s1_std:.3f} eV ({s1_wavelength:.1f} nm)")
    
    print("\nDone!")


def main():
    """Command line interface."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python -m q4mdkit.analysis.tdorca <config.yaml>")
        sys.exit(1)
    
    config_path = sys.argv[1]
    run_tdorca_analysis(config_path)


if __name__ == "__main__":
    main()
