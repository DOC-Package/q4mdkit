#!/usr/bin/env python3
"""
Create filtered MD trajectory with selected normal modes removed.

This script removes the contribution of selected normal modes from the trajectory:
1. Loads MD trajectory and reference structure
2. Performs axis-switching (QCP alignment) for each frame
3. Projects displacements onto normal modes
4. Removes selected mode contributions
5. Reconstructs trajectory without those modes

The filtered trajectory satisfies:
    r_a^(filtered)(t) = R(t)^T @ [r_bar_a^(0) + delta_r^(filtered)(t)] + COM^(0)

where delta_r^(filtered)(t) has the selected mode contributions removed.

All coordinates are handled in Angstrom units.

Usage:
    python run_filter_trajectory.py
    
    # Or with command-line arguments:
    python run_filter_trajectory.py --modes 7 8 9 --output filtered_traj.xyz
"""

import os
import sys
import argparse
import numpy as np
from q4mdkit.analysis.normal_mode_analysis import (
    NormalModeAnalyzer,
    load_trajectory_mdtraj,
    create_mode_filtered_trajectory,
    save_filtered_trajectory_xyz,
)


def parse_mode_list(mode_str: str) -> list:
    """
    Parse mode specification string.
    
    Supports:
    - Single modes: "7"
    - Comma-separated: "7,8,9"
    - Ranges: "7-12"
    - Mixed: "7-9,15,20-22"
    
    Returns 0-indexed mode indices.
    """
    modes = []
    for part in mode_str.split(','):
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            modes.extend(range(int(start), int(end) + 1))
        else:
            modes.append(int(part))
    
    # Convert to 0-indexed (assume input is 1-indexed mode numbers like 7, 8, 9...)
    # Mode 7 corresponds to index 6 (first vibrational mode after 6 trans/rot)
    return [m - 1 for m in modes]


def main():
    parser = argparse.ArgumentParser(
        description="Create filtered trajectory with selected normal modes removed"
    )
    parser.add_argument("--dcd", "-d", default="mol1.dcd",
                        help="Input DCD trajectory file")
    parser.add_argument("--top", "-t", default="mol1.pdb",
                        help="Topology file (PDB)")
    parser.add_argument("--normal-mode-dir", "-n", default="normal_mode",
                        help="Directory containing opt.gen, hessian_eigenvectors.txt, etc.")
    parser.add_argument("--modes", "-m", default="7-12",
                        help="Modes to remove (1-indexed). E.g., '7-12' or '7,8,9,15'")
    parser.add_argument("--output", "-o", default="filtered_traj",
                        help="Output prefix (will create .xyz and .npz files)")
    parser.add_argument("--save-dcd", action="store_true",
                        help="Also save as DCD format (requires MDTraj)")
    
    args = parser.parse_args()
    
    # === File paths ===
    dcd_file = args.dcd
    top_file = args.top
    normal_mode_dir = args.normal_mode_dir
    
    gen_file = os.path.join(normal_mode_dir, "opt.gen")
    eigenvectors_file = os.path.join(normal_mode_dir, "hessian_eigenvectors.txt")
    eigenvalues_file = os.path.join(normal_mode_dir, "hessian_eigenvalues.txt")
    
    output_prefix = args.output
    
    # Parse modes to remove
    modes_to_remove = parse_mode_list(args.modes)
    
    # === Load data ===
    print("=" * 70)
    print("Filtered Trajectory Generation")
    print("=" * 70)
    print("Note: All coordinates in Angstrom units")
    print()
    
    # Check files exist
    for f in [dcd_file, top_file, gen_file, eigenvectors_file]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    # Create analyzer from files
    print(f"Loading reference structure: {gen_file}")
    print(f"Loading eigenvectors: {eigenvectors_file}")
    
    eigenvalues_path = eigenvalues_file if os.path.exists(eigenvalues_file) else None
    if eigenvalues_path:
        print(f"Loading eigenvalues: {eigenvalues_file}")
    
    analyzer = NormalModeAnalyzer.from_files(
        gen_file=gen_file,
        eigenvectors_file=eigenvectors_file,
        eigenvalues_file=eigenvalues_path,
        mode_range=None
    )
    
    print(f"\n  Reference atoms: {analyzer.n_atoms}")
    print(f"  Atom types: {set(analyzer.atom_types)}")
    print(f"  Total DOF: {analyzer.n_dof}")
    print(f"  Vibrational modes: {analyzer.n_modes}")
    
    # Load trajectory
    print(f"\nLoading trajectory: {dcd_file}")
    print(f"  Topology: {top_file}")
    
    traj_coords, n_frames = load_trajectory_mdtraj(dcd_file, top_file, unit='angstrom')
    
    print(f"  Frames: {n_frames}")
    print(f"  Atoms in trajectory: {traj_coords.shape[1]}")
    
    # Check atom count match
    if traj_coords.shape[1] != analyzer.n_atoms:
        print(f"\n*** Error: Atom count mismatch! ***")
        print(f"  Trajectory has {traj_coords.shape[1]} atoms")
        print(f"  Reference has {analyzer.n_atoms} atoms")
        sys.exit(1)
    
    # === Create filtered trajectory ===
    print("\n" + "=" * 70)
    print("Creating filtered trajectory")
    print("=" * 70)
    
    # Display modes to remove (1-indexed for user)
    modes_1indexed = [m + 1 for m in modes_to_remove]
    print(f"\nModes to remove: {modes_1indexed}")
    print(f"  (0-indexed internal: {modes_to_remove})")
    
    # Get frequencies for these modes if available
    if analyzer.frequencies_cm is not None:
        print("\nFrequencies of removed modes:")
        for m in modes_to_remove:
            if 0 <= m < len(analyzer.frequencies_cm):
                print(f"  Mode {m+1}: {analyzer.frequencies_cm[m]:.2f} cm^-1")
    
    print("\nPerforming filtering...")
    
    result = analyzer.create_filtered_trajectory(traj_coords, modes_to_remove)
    
    filtered_traj = result['filtered_trajectory']
    rmsd_orig = result['rmsd_original']
    rmsd_filt = result['rmsd_filtered']
    removed_mode_coords = result['removed_mode_coords']
    
    print(f"\nFiltering complete!")
    print(f"\n  RMSD from reference (original trajectory):")
    print(f"    Mean: {np.mean(rmsd_orig):.6f} Å")
    print(f"    Std:  {np.std(rmsd_orig):.6f} Å")
    print(f"    Max:  {np.max(rmsd_orig):.6f} Å")
    
    print(f"\n  RMSD from reference (filtered trajectory):")
    print(f"    Mean: {np.mean(rmsd_filt):.6f} Å")
    print(f"    Std:  {np.std(rmsd_filt):.6f} Å")
    print(f"    Max:  {np.max(rmsd_filt):.6f} Å")
    
    print(f"\n  RMSD reduction: {(1 - np.mean(rmsd_filt)/np.mean(rmsd_orig))*100:.1f}%")
    
    # Statistics of removed mode coordinates
    print(f"\n  Removed mode coordinate statistics:")
    print(f"    {'Mode':>6} {'Mean':>15} {'Std':>15} {'Max|d|':>15}")
    print("    " + "-" * 55)
    for i, m in enumerate(modes_to_remove):
        coords = removed_mode_coords[:, i]
        print(f"    {m+1:6d} {np.mean(coords):15.6e} {np.std(coords):15.6e} {np.max(np.abs(coords)):15.6e}")
    
    # === Save results ===
    print("\n" + "=" * 70)
    print("Saving results")
    print("=" * 70)
    
    # Save filtered trajectory as XYZ
    xyz_file = f"{output_prefix}.xyz"
    save_filtered_trajectory_xyz(
        xyz_file,
        filtered_traj,
        analyzer.atom_types,
        comment_prefix=f"filtered (modes {','.join(map(str, modes_1indexed))} removed)"
    )
    
    # Save as DCD if requested
    if args.save_dcd:
        try:
            import mdtraj as md
            dcd_out = f"{output_prefix}.dcd"
            
            # Create trajectory object and save
            # Load original topology
            orig_traj = md.load(dcd_file, top=top_file, frame=0)
            topology = orig_traj.topology
            
            # Convert Angstrom to nm for MDTraj
            filtered_nm = filtered_traj / 10.0
            
            # Create new trajectory
            new_traj = md.Trajectory(filtered_nm, topology)
            new_traj.save_dcd(dcd_out)
            print(f"  -> {dcd_out}")
        except ImportError:
            print("  Warning: MDTraj not available, skipping DCD output")
        except Exception as e:
            print(f"  Warning: Failed to save DCD: {e}")
    
    # Save analysis data as npz
    npz_file = f"{output_prefix}.npz"
    np.savez(npz_file,
             filtered_trajectory=filtered_traj,
             original_trajectory=traj_coords,
             modes_removed=np.array(modes_to_remove),
             modes_removed_1indexed=np.array(modes_1indexed),
             removed_mode_coords=removed_mode_coords,
             rmsd_original=rmsd_orig,
             rmsd_filtered=rmsd_filt,
             atom_types=analyzer.atom_types,
             ref_coords=analyzer.ref_coords)
    print(f"  -> {npz_file}")
    
    # Save RMSD comparison
    rmsd_file = f"{output_prefix}_rmsd.txt"
    with open(rmsd_file, 'w') as f:
        f.write(f"# RMSD comparison: original vs filtered trajectory\n")
        f.write(f"# Modes removed: {modes_1indexed}\n")
        f.write(f"# frame  rmsd_original(A)  rmsd_filtered(A)  reduction(%)\n")
        for i in range(n_frames):
            reduction = (1 - rmsd_filt[i] / rmsd_orig[i]) * 100 if rmsd_orig[i] > 0 else 0
            f.write(f"{i:8d} {rmsd_orig[i]:18.8e} {rmsd_filt[i]:18.8e} {reduction:12.2f}\n")
    print(f"  -> {rmsd_file}")
    
    # Save removed mode coordinates
    modes_file = f"{output_prefix}_removed_modes.txt"
    with open(modes_file, 'w') as f:
        f.write(f"# Removed mode coordinates\n")
        f.write(f"# Modes: {modes_1indexed}\n")
        f.write("# frame " + " ".join([f"d_{m+1:03d}" for m in modes_to_remove]) + "\n")
        for i in range(n_frames):
            line = f"{i:8d}"
            for d in removed_mode_coords[i]:
                line += f" {d:15.8e}"
            f.write(line + "\n")
    print(f"  -> {modes_file}")
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)
    print(f"\nFiltered trajectory saved to: {xyz_file}")
    print(f"Use this trajectory for further analysis without modes {modes_1indexed}")


if __name__ == "__main__":
    main()
