#!/usr/bin/env python
"""
Calculate energy difference between cation and neutral states.

Reads:
- elstat_energy_cation.dat: Cation energies (multiple molecule columns)
- elstat_energy_neutral.dat: Neutral energies (multiple molecule columns)

Outputs:
- energy_diff.dat: Energy difference (Cation - Neutral) in eV for each molecule
"""

import numpy as np
import argparse

def read_energies(filename):
    """Read energies from file (multiple columns)
    
    Returns:
        frames: array of frame numbers
        times: array of time values
        energies: 2D array (n_frames, n_molecules)
        headers: list of header lines
        col_names: list of column names for each molecule
    """
    frames = []
    times = []
    energies = []
    headers = []
    col_names = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.rstrip())
                # Parse column names from header (looking for "E_mol*" or similar)
                if 'E_mol' in line or 'Frame' in line:
                    parts = line.replace('#', '').split()
                    # Extract molecule column names (skip Frame and Time)
                    col_names = [p for p in parts if p.startswith('E_mol') or p.startswith('dE_mol')]
                continue
            if line.strip() == '':
                continue
            
            cols = line.split()
            # Need at least 3 columns (frame, time, 1 energy)
            if len(cols) < 3:
                continue
            try:
                frame = int(cols[0])
                time = float(cols[1])
                energy_vals = [float(c) for c in cols[2:]]
                frames.append(frame)
                times.append(time)
                energies.append(energy_vals)
            except (ValueError, IndexError):
                continue
    
    return np.array(frames), np.array(times), np.array(energies), headers, col_names

def main():
    parser = argparse.ArgumentParser(description='Calculate energy difference between cation and neutral states')
    parser.add_argument('--cation', '-c', default='elstat_energy_cation.dat',
                        help='Cation energy file (default: elstat_energy_cation.dat)')
    parser.add_argument('--neutral', '-n', default='elstat_energy_neutral.dat',
                        help='Neutral energy file (default: elstat_energy_neutral.dat)')
    parser.add_argument('--output', '-o', default='energy_diff.dat',
                        help='Output file (default: energy_diff.dat)')
    args = parser.parse_args()
    
    # Read input files
    cation_frames, cation_times, cation_energies, cation_headers, col_names = read_energies(args.cation)
    neutral_frames, neutral_times, neutral_energies, neutral_headers, _ = read_energies(args.neutral)
    
    # Check that the number of data points match
    n_cation = len(cation_energies)
    n_neutral = len(neutral_energies)
    
    if n_cation != n_neutral:
        print(f"Warning: Number of data points differ: cation={n_cation}, neutral={n_neutral}")
        n_points = min(n_cation, n_neutral)
    else:
        n_points = n_cation
    
    # Check that the number of columns match
    n_cols_cation = cation_energies.shape[1] if len(cation_energies.shape) > 1 else 1
    n_cols_neutral = neutral_energies.shape[1] if len(neutral_energies.shape) > 1 else 1
    
    if n_cols_cation != n_cols_neutral:
        print(f"Error: Number of energy columns differ: cation={n_cols_cation}, neutral={n_cols_neutral}")
        return
    
    n_cols = n_cols_cation
    
    # Verify time alignment
    time_diff = np.abs(cation_times[:n_points] - neutral_times[:n_points])
    if np.max(time_diff) > 1e-6:
        print("Warning: Time values do not match between files!")
        print(f"Max time difference: {np.max(time_diff)} fs")
    
    # Calculate energy difference (Cation - Neutral) for each molecule column
    dE = cation_energies[:n_points] - neutral_energies[:n_points]
    
    # Generate column names if not found in header
    if not col_names or len(col_names) != n_cols:
        col_names = [f'dE_col{i+1}' for i in range(n_cols)]
    else:
        # Convert E_mol* to dE_mol*
        col_names = [name.replace('E_mol', 'dE_mol') for name in col_names]
    
    # Write output file
    output_file = args.output
    with open(output_file, 'w') as f:
        f.write("# Energy Difference (Cation - Neutral)\n")
        # Copy molecule IDs line if present
        for h in cation_headers:
            if 'Molecule IDs' in h:
                f.write(h + '\n')
                break
        
        # Write column header
        header_line = "# Frame  Time(fs)  "
        header_line += "  ".join([f"{name:>14}" for name in col_names])
        f.write(header_line + "\n")
        
        for i in range(n_points):
            frame = cation_frames[i]
            time = cation_times[i]
            line = f"{frame:5d}  {time:12.3f}"
            for j in range(n_cols):
                line += f"  {dE[i, j]:14.10f}"
            f.write(line + "\n")
    
    print(f"Energy difference calculated for {n_points} frames, {n_cols} molecules")
    print(f"Output written to: {output_file}")
    
    # Print statistics for each molecule
    print("\nStatistics (eV):")
    print(f"{'Column':<15} {'Mean':>12} {'Std':>12}")
    print("-" * 41)
    for j in range(n_cols):
        mean_val = np.mean(dE[:n_points, j])
        std_val = np.std(dE[:n_points, j])
        print(f"{col_names[j]:<15} {mean_val:>12.6f} {std_val:>12.6f}")

if __name__ == '__main__':
    main()
