#!/usr/bin/env python3
"""
Compute difference of Mulliken charges between cation and neutral states.

Reads mulliken_charges_cation.dat and mulliken_charges_neutral.dat,
computes the difference (cation - neutral), and outputs the result.

This represents the excess charge distribution in the cation state.
"""

import numpy as np
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================
CATION_FILE = 'mulliken_charges_cation.dat'
NEUTRAL_FILE = 'mulliken_charges_neutral.dat'
OUTPUT_FILE = 'mulliken_charges_diff.dat'


def read_mulliken_charges(filepath: str) -> tuple:
    """
    Read Mulliken charges from file.
    
    Expected format:
    # Mulliken charges for each frame
    # Atom types: C C C ...
    # Columns: Frame, Time(fs), Q_1, Q_2, ..., Q_N
    # Units: elementary charge (e)
    5     12.300  0.12345678  0.23456789  ...
    
    Returns
    -------
    frames : np.ndarray
        Frame indices
    times : np.ndarray
        Time values (fs)
    charges : np.ndarray
        Charges array (n_frames, n_atoms)
    atom_types : list or None
        Atom types if found in header
    """
    frames = []
    times = []
    charges = []
    atom_types = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line_strip = line.strip()
            
            # Parse header for atom types
            if line_strip.startswith('# Atom types:'):
                atom_types = line_strip.replace('# Atom types:', '').split()
                continue
            
            # Skip other comments
            if not line_strip or line_strip.startswith('#'):
                continue
            
            parts = line_strip.split()
            if len(parts) >= 3:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                charges.append([float(x) for x in parts[2:]])
    
    frames = np.array(frames)
    times = np.array(times)
    charges = np.array(charges)
    
    print(f"  Loaded {filepath}: {len(frames)} frames, {charges.shape[1]} atoms")
    
    return frames, times, charges, atom_types


def main():
    print("=" * 70)
    print("Mulliken Charge Difference: Cation - Neutral")
    print("=" * 70)
    
    # Check files exist
    if not Path(CATION_FILE).exists():
        raise FileNotFoundError(f"Cation file not found: {CATION_FILE}")
    if not Path(NEUTRAL_FILE).exists():
        raise FileNotFoundError(f"Neutral file not found: {NEUTRAL_FILE}")
    
    # Read data
    print("\nReading charge files...")
    frames_c, times_c, charges_c, atom_types = read_mulliken_charges(CATION_FILE)
    frames_n, times_n, charges_n, _ = read_mulliken_charges(NEUTRAL_FILE)
    
    # Validate
    if len(frames_c) != len(frames_n):
        raise ValueError(f"Frame count mismatch: cation={len(frames_c)}, neutral={len(frames_n)}")
    
    if charges_c.shape[1] != charges_n.shape[1]:
        raise ValueError(f"Atom count mismatch: cation={charges_c.shape[1]}, neutral={charges_n.shape[1]}")
    
    n_frames = len(frames_c)
    n_atoms = charges_c.shape[1]
    
    # Compute difference
    print("\nComputing difference (cation - neutral)...")
    charges_diff = charges_c - charges_n
    
    # Statistics
    diff_mean = np.mean(charges_diff, axis=0)  # Mean over frames for each atom
    diff_std = np.std(charges_diff, axis=0)
    total_excess = np.sum(charges_diff, axis=1)  # Sum over atoms for each frame
    
    print(f"\n=== Statistics ===")
    print(f"  Frames: {n_frames}")
    print(f"  Atoms: {n_atoms}")
    print(f"  Total excess charge (mean): {np.mean(total_excess):.6f} e")
    print(f"  Total excess charge (std):  {np.std(total_excess):.6f} e")
    
    # Per-atom statistics
    print(f"\n=== Per-atom charge difference (mean ± std) ===")
    for i in range(n_atoms):
        atom_label = atom_types[i] if atom_types else f"Atom{i}"
        print(f"  {i:3d} ({atom_label:2s}): {diff_mean[i]:+.6f} ± {diff_std[i]:.6f} e")
    
    # Write output
    print(f"\nWriting: {OUTPUT_FILE}")
    with open(OUTPUT_FILE, 'w') as f:
        f.write("# Mulliken charge difference: cation - neutral\n")
        if atom_types:
            f.write(f"# Atom types: {' '.join(atom_types)}\n")
        f.write(f"# Columns: Frame, Time(fs), dQ_1, dQ_2, ..., dQ_N, Total\n")
        f.write(f"# Units: elementary charge (e)\n")
        
        for i in range(n_frames):
            frame = frames_c[i]
            time = times_c[i]
            charges_str = "  ".join(f"{q:12.8f}" for q in charges_diff[i])
            total = total_excess[i]
            f.write(f"{frame:5d}  {time:12.3f}  {charges_str}  {total:12.8f}\n")
    
    # Write summary statistics
    summary_file = OUTPUT_FILE.replace('.dat', '_summary.dat')
    print(f"Writing: {summary_file}")
    with open(summary_file, 'w') as f:
        f.write("# Mulliken charge difference summary (time-averaged)\n")
        f.write("# Atom  Type  Mean_dQ(e)      Std_dQ(e)\n")
        for i in range(n_atoms):
            atom_label = atom_types[i] if atom_types else "?"
            f.write(f"{i:5d}  {atom_label:4s}  {diff_mean[i]:+12.8f}  {diff_std[i]:12.8f}\n")
        f.write(f"# Total excess charge: {np.mean(total_excess):.6f} ± {np.std(total_excess):.6f} e\n")
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
