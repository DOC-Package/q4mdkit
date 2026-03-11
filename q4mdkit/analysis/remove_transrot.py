#!/usr/bin/env python3
"""
Remove translation and rotation from trajectory.

Superpose all frames onto the first frame using MDTraj.
This removes rigid-body motion (translation + rotation).

Usage:
    python remove_transrot.py --dcd mol.dcd --top mol.pdb --output mol_aligned.dcd
    python remove_transrot.py --dcd mol.dcd --top mol.pdb  # outputs mol_aligned.dcd by default
"""

import argparse
import mdtraj as md
from pathlib import Path


def remove_translation_rotation(dcd_file: str, top_file: str, 
                                 output_file: str = None,
                                 atom_indices: list = None) -> md.Trajectory:
    """
    Remove translation and rotation from trajectory.
    
    Parameters
    ----------
    dcd_file : str
        Input DCD trajectory file
    top_file : str
        Topology file (PDB, GRO, etc.)
    output_file : str, optional
        Output DCD file path. If None, uses input_aligned.dcd
    atom_indices : list, optional
        Atom indices to use for alignment. If None, uses all atoms.
    
    Returns
    -------
    traj : md.Trajectory
        Aligned trajectory
    """
    # Load trajectory
    print(f"Loading trajectory: {dcd_file}")
    print(f"Topology: {top_file}")
    traj = md.load(dcd_file, top=top_file)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    # Use first frame as reference
    reference = traj[0]
    print(f"\nUsing frame 0 as reference")
    
    # Superpose: removes translation and rotation
    # atom_indices specifies which atoms to use for alignment calculation
    print("Aligning trajectory (removing translation and rotation)...")
    if atom_indices is not None:
        traj.superpose(reference, atom_indices=atom_indices)
        print(f"  Using {len(atom_indices)} atoms for alignment")
    else:
        traj.superpose(reference)
        print(f"  Using all {traj.n_atoms} atoms for alignment")
    
    # Determine output filename
    if output_file is None:
        input_path = Path(dcd_file)
        output_file = str(input_path.parent / f"{input_path.stem}_aligned.dcd")
    
    # Save aligned trajectory
    print(f"\nSaving aligned trajectory: {output_file}")
    traj.save_dcd(output_file)
    print("Done!")
    
    return traj


def main():
    parser = argparse.ArgumentParser(
        description='Remove translation and rotation from trajectory')
    parser.add_argument('--dcd', '-d', required=True, 
                        help='Input DCD trajectory file')
    parser.add_argument('--top', '-t', required=True, 
                        help='Topology file (PDB, GRO, etc.)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output DCD file (default: input_aligned.dcd)')
    parser.add_argument('--atoms', '-a', type=int, nargs='+', default=None,
                        help='Atom indices for alignment (0-indexed)')
    
    args = parser.parse_args()
    
    remove_translation_rotation(
        dcd_file=args.dcd,
        top_file=args.top,
        output_file=args.output,
        atom_indices=args.atoms
    )


if __name__ == '__main__':
    main()
