#!/usr/bin/env python3
"""
Compute average structure from trajectory.

Usage:
    python average_structure.py --dcd mol.dcd --top mol.pdb --output mol_avg.pdb
    python average_structure.py --dcd mol_aligned.dcd --top mol.pdb  # outputs mol_aligned_average.pdb
"""

import argparse
import re
import numpy as np
import mdtraj as md
from pathlib import Path


def save_xyz_element_only(traj: md.Trajectory, filename: str, comment: str = ""):
    """
    Save trajectory as XYZ file with element symbols only (no numbers).
    
    Parameters
    ----------
    traj : md.Trajectory
        Trajectory to save (single frame)
    filename : str
        Output XYZ filename
    comment : str
        Comment line
    """
    # Get element symbols (strip numbers)
    elements = [re.sub(r'\d+', '', atom.name) for atom in traj.topology.atoms]
    
    # Convert nm to Angstrom
    coords = traj.xyz[0] * 10.0
    
    with open(filename, 'w') as f:
        f.write(f"{traj.n_atoms}\n")
        f.write(f"{comment}\n")
        for elem, (x, y, z) in zip(elements, coords):
            f.write(f"{elem:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def compute_average_structure(dcd_file: str, top_file: str, 
                               output_file: str = None) -> np.ndarray:
    """
    Compute average structure from trajectory.
    
    Parameters
    ----------
    dcd_file : str
        Input DCD trajectory file (should be pre-aligned for meaningful average)
    top_file : str
        Topology file (PDB, GRO, etc.)
    output_file : str, optional
        Output PDB file path for average structure
    
    Returns
    -------
    avg_coords : np.ndarray
        Average coordinates (n_atoms, 3) in nm
    """
    # Load trajectory
    print(f"Loading trajectory: {dcd_file}")
    print(f"Topology: {top_file}")
    traj = md.load(dcd_file, top=top_file)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    print(f"\nComputing average structure from {traj.n_frames} frames...")
    
    # Compute mean coordinates
    avg_coords = np.mean(traj.xyz, axis=0)
    
    # Determine output filename
    if output_file is None:
        input_path = Path(dcd_file)
        output_file = str(input_path.parent / f"{input_path.stem}_average.pdb")
    
    # Create trajectory with single frame containing average structure
    avg_traj = traj.slice(0)
    avg_traj.xyz[0] = avg_coords
    
    # Save as PDB
    print(f"Saving average structure: {output_file}")
    avg_traj.save_pdb(output_file)
    
    # Save as XYZ
    output_path = Path(output_file)
    xyz_file = str(output_path.parent / f"{output_path.stem}.xyz")
    print(f"Saving average structure: {xyz_file}")
    save_xyz_element_only(avg_traj, xyz_file, comment="Average structure")
    
    print("Done!")
    
    return avg_coords


def compute_average_from_traj(traj: md.Trajectory, output_file: str = None) -> np.ndarray:
    """
    Compute average structure from MDTraj Trajectory object.
    
    Parameters
    ----------
    traj : md.Trajectory
        Input trajectory (should be pre-aligned)
    output_file : str, optional
        Output PDB file path for average structure
    
    Returns
    -------
    avg_coords : np.ndarray
        Average coordinates (n_atoms, 3) in nm
    """
    print(f"\nComputing average structure from {traj.n_frames} frames...")
    
    # Compute mean coordinates
    avg_coords = np.mean(traj.xyz, axis=0)
    
    if output_file is not None:
        # Create trajectory with single frame containing average structure
        avg_traj = traj.slice(0)
        avg_traj.xyz[0] = avg_coords
        
        # Save as PDB
        print(f"Saving average structure: {output_file}")
        avg_traj.save_pdb(output_file)
        
        # Save as XYZ
        output_path = Path(output_file)
        xyz_file = str(output_path.parent / f"{output_path.stem}.xyz")
        print(f"Saving average structure: {xyz_file}")
        save_xyz_element_only(avg_traj, xyz_file, comment="Average structure")
    
    return avg_coords


def main():
    parser = argparse.ArgumentParser(
        description='Compute average structure from trajectory')
    parser.add_argument('--dcd', '-d', required=True, 
                        help='Input DCD trajectory file')
    parser.add_argument('--top', '-t', required=True, 
                        help='Topology file (PDB, GRO, etc.)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output PDB file (default: input_average.pdb)')
    
    args = parser.parse_args()
    
    compute_average_structure(
        dcd_file=args.dcd,
        top_file=args.top,
        output_file=args.output
    )


if __name__ == '__main__':
    main()
