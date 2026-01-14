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
from dataclasses import dataclass
from typing import List, Optional

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
    run_dftb_in_subprocess,
)


@dataclass
class DFTBConfig:
    """Configuration for standard DFTB trajectory analysis."""
    # Input paths
    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    pccharges_template: Path
    hsd_template: Path
    
    # Output paths
    output_dir: Path
    energy_file: Path
    work_directory: str
    
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
    
    return DFTBConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        hsd_template=hsd_template,
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
    read_initial_charges: bool = False,
) -> Path:
    """
    Set up HSD file for closed-shell, non-constrained DFTB calculation.
    
    Removes SpinPolarisation and ElectronicConstraints from template.
    """
    # Load and modify HSD template
    with open(hsd_template, 'r') as f:
        data = hsd.load(f)
    
    # Remove Geometry from data (will be added as text reference)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Remove SpinPolarisation for closed-shell calculation
    if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
        del data['Hamiltonian']['DFTB']['SpinPolarisation']
    
    # Remove SpinConstants if present
    if 'SpinConstants' in data['Hamiltonian']['DFTB']:
        del data['Hamiltonian']['DFTB']['SpinConstants']
    
    # Remove ElectronicConstraints for non-constrained DFTB
    if 'ElectronicConstraints' in data['Hamiltonian']['DFTB']:
        del data['Hamiltonian']['DFTB']['ElectronicConstraints']
    
    # Set ReadInitialCharges
    data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'Yes' if read_initial_charges else 'No'
    
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
    print("DFTB Trajectory Analysis (Closed-shell, Non-constrained)")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Trajectory: {config.traj_path}")
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
    print()
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up work directory path
    work_dir = config.output_dir / config.work_directory
    
    # Open output file
    energy_file = open(config.energy_file, "w")
    energy_file.write("# DFTB Trajectory Analysis (Closed-shell, Non-constrained)\n")
    energy_file.write("# Frame  Time(fs)        Energy(a.u.)\n")
    
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
            
            print(f"Frame {frame_id:05d} (t = {time_val:.3f} fs)")
            
            # Set up work directory with current coordinates
            setup_work_directory(
                work_dir, qm_coords_ang, mm_coords_ang, mm_charges,
                frame_id, time_fs, atom_types
            )
            
            # Determine if we should use initial charges
            use_initial_charges = False
            charges_file_path = work_dir / "charges.dat"
            
            if config.use_previous_charges and processed_count > 0:
                if charges_file_path.exists():
                    use_initial_charges = True
            
            # Set up HSD file (closed-shell, no constraints)
            setup_dftb_hsd(work_dir, config.hsd_template, read_initial_charges=use_initial_charges)
            
            # Run DFTB+ calculation
            energy, mcharge, error = run_dftb_in_subprocess(
                qm_coords_bohr, work_dir, write_hs=False,
                dftb_library_path=config.dftb_library_path,
                num_threads=config.num_threads,
                timeout=config.timeout
            )
            
            if error:
                print(f"  ERROR: {error}")
                energy = float("nan")
            else:
                print(f"  E = {energy:.10f} a.u.")
            
            # Save results
            energy_file.write(f"{frame_id:5d}  {time_val:12.3f}  {energy:18.10f}\n")
            energy_file.flush()
            
            processed_count += 1
    
    finally:
        energy_file.close()
    
    print("\n" + "=" * 70)
    print(f"Energies saved to {config.energy_file}")
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
