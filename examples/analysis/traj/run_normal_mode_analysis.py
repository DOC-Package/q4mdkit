#!/usr/bin/env python3
"""
Run normal mode coordinate analysis on mol1.dcd trajectory.

This script extracts normal mode coordinates from MD trajectory by:
1. Loading MD trajectory (mol1.dcd) 
2. Loading reference structure (opt.gen)
3. Loading Hessian eigenvectors and eigenvalues
4. Performing axis-switching (QCP alignment) for each frame
5. Projecting onto reference normal modes

All coordinates are handled in Angstrom units.

Usage:
    python run_normal_mode_analysis.py
"""

import os
import sys
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from q4mdkit.analysis.normal_mode_analysis import (
    NormalModeAnalyzer,
    load_trajectory_mdtraj,
    read_gen_file,
    calculate_mode_correlation,
    plot_mode_correlation,
    plot_strong_correlations,
    print_correlation_summary,
    rotation_angle,
)


def main():
    # === File paths ===
    # Trajectory
    dcd_file = "mol1.dcd"
    top_file = "mol1.pdb"
    
    # Reference structure and Hessian data
    normal_mode_dir = "normal_mode"
    gen_file = os.path.join(normal_mode_dir, "opt.gen")
    eigenvectors_file = os.path.join(normal_mode_dir, "hessian_eigenvectors.txt")
    eigenvalues_file = os.path.join(normal_mode_dir, "hessian_eigenvalues.txt")
    
    # Output
    output_prefix = "mode_coords"
    
    # === Load data ===
    print("=" * 60)
    print("Normal Mode Coordinate Analysis")
    print("=" * 60)
    print("Note: All coordinates in Angstrom units")
    
    # Check files exist
    for f in [dcd_file, top_file, gen_file, eigenvectors_file, eigenvalues_file]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    # Create analyzer from files (Angstrom units)
    print(f"\nLoading reference structure: {gen_file}")
    print(f"Loading eigenvectors: {eigenvectors_file}")
    print(f"Loading eigenvalues: {eigenvalues_file}")
    
    analyzer = NormalModeAnalyzer.from_files(
        gen_file=gen_file,
        eigenvectors_file=eigenvectors_file,
        eigenvalues_file=eigenvalues_file,
        mode_range=None  # Use default (skip first 6 modes)
    )
    
    print(f"\n  Reference atoms: {analyzer.n_atoms}")
    print(f"  Atom types: {set(analyzer.atom_types)}")
    print(f"  Total DOF: {analyzer.n_dof}")
    print(f"  Vibrational modes: {analyzer.n_modes} (indices {analyzer.mode_range[0]} to {analyzer.mode_range[1]})")
    
    # Load trajectory (Angstrom units)
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
    
    # === Perform analysis ===
    print("\n" + "=" * 60)
    print("Running analysis with QCP alignment...")
    print("=" * 60)
    
    results = analyzer.analyze_trajectory(traj_coords, align=True)
    
    mode_coords = results['mode_coords']
    rmsd = results['rmsd']
    rot_angles = results['rotation_angles']
    
    print(f"\nResults:")
    print(f"  Mode coordinates shape: {mode_coords.shape}")
    print(f"  RMSD from reference (Å):")
    print(f"    Mean: {np.mean(rmsd):.6f}")
    print(f"    Std:  {np.std(rmsd):.6f}")
    print(f"    Max:  {np.max(rmsd):.6f}")
    print(f"  Rotation angle (degrees):")
    print(f"    Mean: {np.mean(rot_angles):.4f}")
    print(f"    Std:  {np.std(rot_angles):.4f}")
    print(f"    Max:  {np.max(rot_angles):.4f}")
    
    # === Statistics ===
    print("\n" + "=" * 60)
    print("Mode statistics")
    print("=" * 60)
    
    stats = analyzer.calculate_statistics(mode_coords, temperature=300.0)
    frequencies = analyzer.get_mode_frequencies()
    
    print(f"\nFirst 10 vibrational modes (mode 7-16):")
    print(f"{'Mode':>6} {'Freq(cm-1)':>12} {'Mean':>15} {'Std':>15} {'Var':>15}")
    print("-" * 65)
    
    for i in range(min(10, analyzer.n_modes)):
        mode_idx = analyzer.mode_indices[i] + 1  # 1-indexed
        freq = frequencies[i] if frequencies is not None else 0.0
        print(f"{mode_idx:6d} {freq:12.2f} {stats['mean'][i]:15.6e} "
              f"{stats['std'][i]:15.6e} {stats['variance'][i]:15.6e}")
    
    # === Save results ===
    print("\n" + "=" * 60)
    print("Saving results")
    print("=" * 60)
    
    # Save as npz
    npz_file = f"{output_prefix}.npz"
    np.savez(npz_file,
             mode_coords=mode_coords,
             rmsd=rmsd,
             rotation_angles=rot_angles,
             mode_indices=analyzer.mode_indices,
             frequencies_cm=frequencies,
             mean=stats['mean'],
             std=stats['std'],
             variance=stats['variance'])
    print(f"  -> {npz_file}")
    
    # Save mode coordinates as text
    txt_file = f"{output_prefix}.txt"
    with open(txt_file, 'w') as f:
        f.write(f"# Normal mode coordinates from {dcd_file}\n")
        f.write(f"# Reference: {gen_file}\n")
        f.write(f"# Frames: {n_frames}, Modes: {analyzer.n_modes}\n")
        f.write(f"# Mode indices: {analyzer.mode_range[0]+1} to {analyzer.mode_range[1]} (1-indexed)\n")
        f.write(f"# Columns: frame_index, d_7, d_8, ..., d_3N\n")
        f.write("#\n")
        for i in range(n_frames):
            line = f"{i:8d}"
            for d in mode_coords[i]:
                line += f" {d:15.8e}"
            f.write(line + "\n")
    print(f"  -> {txt_file}")
    
    # Save RMSD and rotation angles
    rmsd_file = f"{output_prefix}_rmsd.txt"
    with open(rmsd_file, 'w') as f:
        f.write(f"# RMSD and rotation angles from reference structure\n")
        f.write(f"# frame  rmsd_angstrom  rotation_deg\n")
        for i, (r, angle) in enumerate(zip(rmsd, rot_angles)):
            f.write(f"{i:8d} {r:15.8e} {angle:12.6f}\n")
    print(f"  -> {rmsd_file}")
    
    # Save summary statistics
    stats_file = f"{output_prefix}_stats.txt"
    with open(stats_file, 'w') as f:
        f.write(f"# Normal mode statistics\n")
        f.write(f"# Reference: {gen_file}\n")
        f.write(f"# Temperature: 300 K\n")
        f.write("#\n")
        f.write(f"# {'Mode':>5} {'Freq(cm-1)':>12} {'Mean':>15} {'Std':>15} {'Variance':>15}\n")
        for i in range(analyzer.n_modes):
            mode_idx = analyzer.mode_indices[i] + 1
            freq = frequencies[i] if frequencies is not None else 0.0
            f.write(f"  {mode_idx:5d} {freq:12.4f} {stats['mean'][i]:15.8e} "
                   f"{stats['std'][i]:15.8e} {stats['variance'][i]:15.8e}\n")
    print(f"  -> {stats_file}")
    
    # === Correlation analysis ===
    print("\n" + "=" * 60)
    print("Mode Correlation Analysis")
    print("=" * 60)
    
    # Mode indices for labeling (1-indexed)
    mode_labels = analyzer.mode_indices + 1
    
    # Print correlation summary
    print_correlation_summary(mode_coords, mode_labels, top_n=15)
    
    # Plot correlation matrix heatmap
    corr_matrix_file = f"{output_prefix}_correlation_matrix.png"
    plot_mode_correlation(
        mode_coords, 
        mode_labels,
        output_file=corr_matrix_file,
        figsize=(12, 10),
        title="Normal Mode Correlation Matrix"
    )
    
    # Plot scatter plots of strongly correlated pairs
    scatter_file = f"{output_prefix}_strong_correlations.png"
    plot_strong_correlations(
        mode_coords,
        mode_labels,
        threshold=0.7,
        max_pairs=12,
        output_file=scatter_file
    )
    
    # Save correlation data
    corr_result = calculate_mode_correlation(mode_coords, mode_labels)
    corr_file = f"{output_prefix}_correlation.npz"
    np.savez(corr_file,
             correlation_matrix=corr_result['correlation_matrix'],
             mode_indices=corr_result['mode_indices'],
             max_abs=corr_result['statistics']['max_abs'],
             mean_abs=corr_result['statistics']['mean_abs'],
             median_abs=corr_result['statistics']['median_abs'])
    print(f"  -> {corr_file}")
    
    # Save strong correlations as text
    strong_corr_file = f"{output_prefix}_strong_correlations.txt"
    with open(strong_corr_file, 'w') as f:
        f.write("# Strongly correlated mode pairs (|r| > 0.5)\n")
        f.write("# mode_i  mode_j  correlation\n")
        for mode_i, mode_j, r in corr_result['strong_correlations']:
            f.write(f"{mode_i:8d} {mode_j:8d} {r:12.6f}\n")
    print(f"  -> {strong_corr_file}")
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
