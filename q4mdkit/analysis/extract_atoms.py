"""
Extract subset of atoms from trajectory.

Extract trajectory for specific atom indices from a DCD file.
Useful for extracting QM region atoms from QM/MM simulation trajectories.
"""

import os
import argparse
import numpy as np
import mdtraj as md
from typing import Optional, Union
from pathlib import Path


def read_atom_indices(filepath: str) -> list[int]:
    """Read atom indices from a file.
    
    The file should contain space-separated atom indices.
    
    Parameters
    ----------
    filepath : str
        Path to the file containing atom indices.
    
    Returns
    -------
    list[int]
        List of atom indices (0-indexed for mdtraj).
    """
    with open(filepath, 'r') as f:
        content = f.read().strip()
    
    indices = [int(x) for x in content.split()]
    return indices


def extract_atoms(
    trajectory_path: str,
    topology_path: str,
    atom_indices: Union[str, list[int]],
    output_path: str,
    output_topology: Optional[str] = None,
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    stride: int = 1
) -> str:
    """Extract trajectory for specific atoms.
    
    Parameters
    ----------
    trajectory_path : str
        Path to the input trajectory file (DCD format).
    topology_path : str
        Path to the topology file (PDB or GRO format).
    atom_indices : str or list[int]
        Either a path to a file containing atom indices,
        or a list of atom indices directly.
        Indices are 0-indexed (mdtraj convention).
    output_path : str
        Path for the output trajectory file.
    output_topology : str, optional
        Path to save the subset topology (PDB format).
        If None, topology is saved with same name as output_path but .pdb extension.
    start_frame : int, optional
        Starting frame index (0-indexed). Default is 0.
    end_frame : int, optional
        Ending frame index (exclusive). If None, process all frames.
    stride : int, optional
        Step size for frame selection. Default is 1.
    
    Returns
    -------
    str
        Path to the output trajectory file.
    
    Example
    -------
    # Extract atoms specified in qmatoms1 file
    >>> extract_atoms("nve.dcd", "nve.pdb", "qmatoms1", "mol1.dcd")
    
    # Extract atoms by index list
    >>> extract_atoms("nve.dcd", "nve.pdb", [0, 1, 2, 3], "subset.dcd")
    """
    # Parse atom indices
    if isinstance(atom_indices, str):
        if os.path.isfile(atom_indices):
            indices = read_atom_indices(atom_indices)
            print(f"Read {len(indices)} atom indices from: {atom_indices}")
        else:
            # Try parsing as space-separated string
            indices = [int(x) for x in atom_indices.split()]
    else:
        indices = list(atom_indices)
    
    indices = np.array(indices, dtype=int)
    print(f"Extracting {len(indices)} atoms")
    print(f"  Index range: {indices.min()} - {indices.max()}")
    
    # Load trajectory
    print(f"\nLoading trajectory: {trajectory_path}")
    print(f"Topology: {topology_path}")
    traj = md.load(trajectory_path, top=topology_path)
    
    total_frames = traj.n_frames
    n_atoms_total = traj.n_atoms
    print(f"  Total frames: {total_frames}")
    print(f"  Total atoms: {n_atoms_total}")
    
    # Validate indices
    if indices.max() >= n_atoms_total:
        raise ValueError(
            f"Atom index {indices.max()} out of range. "
            f"Trajectory has {n_atoms_total} atoms (0-indexed: 0 to {n_atoms_total-1})"
        )
    if indices.min() < 0:
        raise ValueError(f"Negative atom index: {indices.min()}")
    
    # Select frame range
    if end_frame is None:
        end_frame = total_frames
    
    frame_slice = slice(start_frame, end_frame, stride)
    traj_sliced = traj[frame_slice]
    n_frames_out = traj_sliced.n_frames
    print(f"\nFrame selection: {start_frame} to {end_frame-1} (stride={stride})")
    print(f"  Output frames: {n_frames_out}")
    
    # Extract atoms
    traj_subset = traj_sliced.atom_slice(indices)
    print(f"\nExtracted trajectory:")
    print(f"  Atoms: {traj_subset.n_atoms}")
    print(f"  Frames: {traj_subset.n_frames}")
    
    # Save output trajectory
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    
    output_ext = Path(output_path).suffix.lower()
    if output_ext == '.dcd':
        traj_subset.save_dcd(output_path)
    elif output_ext == '.xyz':
        traj_subset.save_xyz(output_path)
    elif output_ext == '.pdb':
        traj_subset.save_pdb(output_path)
    else:
        traj_subset.save(output_path)
    
    print(f"\nSaved trajectory: {output_path}")
    
    # Save topology
    if output_topology is None:
        output_topology = str(Path(output_path).with_suffix('.pdb'))
    
    traj_subset[0].save_pdb(output_topology)
    print(f"Saved topology: {output_topology}")
    
    return output_path


def extract_qm_molecules(
    trajectory_path: str,
    topology_path: str,
    qmatoms1_file: str,
    qmatoms2_file: Optional[str] = None,
    output_dir: str = ".",
    output_prefix: str = "mol",
    **kwargs
) -> dict:
    """Extract QM molecule trajectories from QM/MM simulation.
    
    Convenience function to extract one or two QM molecules separately.
    
    Parameters
    ----------
    trajectory_path : str
        Path to the input trajectory file.
    topology_path : str
        Path to the topology file.
    qmatoms1_file : str
        Path to file containing atom indices for molecule 1.
    qmatoms2_file : str, optional
        Path to file containing atom indices for molecule 2.
    output_dir : str, optional
        Output directory. Default is current directory.
    output_prefix : str, optional
        Prefix for output files. Default is "mol".
    **kwargs
        Additional arguments passed to extract_atoms (start_frame, end_frame, stride).
    
    Returns
    -------
    dict
        Dictionary with paths to output files:
        {'mol1': path, 'mol2': path (if provided), 'mol1_top': path, ...}
    """
    os.makedirs(output_dir, exist_ok=True)
    result = {}
    
    # Extract molecule 1
    output1 = os.path.join(output_dir, f"{output_prefix}1.dcd")
    output1_top = os.path.join(output_dir, f"{output_prefix}1.pdb")
    
    print("=" * 60)
    print("Extracting molecule 1")
    print("=" * 60)
    extract_atoms(
        trajectory_path=trajectory_path,
        topology_path=topology_path,
        atom_indices=qmatoms1_file,
        output_path=output1,
        output_topology=output1_top,
        **kwargs
    )
    result['mol1'] = output1
    result['mol1_top'] = output1_top
    
    # Extract molecule 2 if provided
    if qmatoms2_file is not None:
        output2 = os.path.join(output_dir, f"{output_prefix}2.dcd")
        output2_top = os.path.join(output_dir, f"{output_prefix}2.pdb")
        
        print("\n" + "=" * 60)
        print("Extracting molecule 2")
        print("=" * 60)
        extract_atoms(
            trajectory_path=trajectory_path,
            topology_path=topology_path,
            atom_indices=qmatoms2_file,
            output_path=output2,
            output_topology=output2_top,
            **kwargs
        )
        result['mol2'] = output2
        result['mol2_top'] = output2_top
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)
    
    return result


def main():
    """CLI for extracting atom subsets from trajectories."""
    parser = argparse.ArgumentParser(
        description="Extract subset of atoms from trajectory file"
    )
    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Path to trajectory file (DCD)"
    )
    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Path to topology file (PDB/GRO)"
    )
    parser.add_argument(
        "-i", "--indices",
        required=True,
        help="File containing atom indices or space-separated indices"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output trajectory path"
    )
    parser.add_argument(
        "--output-topology",
        default=None,
        help="Output topology path (default: same as output with .pdb extension)"
    )
    parser.add_argument(
        "--start-frame",
        type=int,
        default=0,
        help="Starting frame index (0-indexed, default: 0)"
    )
    parser.add_argument(
        "--end-frame",
        type=int,
        default=None,
        help="Ending frame index (exclusive, default: all)"
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Frame stride (default: 1)"
    )
    
    args = parser.parse_args()
    
    extract_atoms(
        trajectory_path=args.trajectory,
        topology_path=args.topology,
        atom_indices=args.indices,
        output_path=args.output,
        output_topology=args.output_topology,
        start_frame=args.start_frame,
        end_frame=args.end_frame,
        stride=args.stride
    )


if __name__ == "__main__":
    main()
