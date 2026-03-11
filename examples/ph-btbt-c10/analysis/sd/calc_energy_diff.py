#!/usr/bin/env python
"""
Calculate energy difference between cation and neutral states.

Reads:
- energies_cation.dat: Cation energies (6th column)
- energies_neutral.dat: Neutral energies (6th column)

Outputs:
- energy_diff.dat: Energy difference (Cation - Neutral) in a.u. and eV
"""

import numpy as np

# Conversion factor from Hartree to eV
HARTREE_TO_EV = 27.211386245988

def read_energies(filename):
    """Read energies from file (6th column)"""
    frames = []
    times = []
    energies = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.split()
            # Skip incomplete rows (need at least 6 columns)
            if len(cols) < 6:
                continue
            try:
                frame = int(cols[0])
                time = float(cols[1])
                energy = float(cols[5])  # 6th column
                frames.append(frame)
                times.append(time)
                energies.append(energy)
            except (ValueError, IndexError):
                continue
    return np.array(frames), np.array(times), np.array(energies)

def main():
    # Read input files
    cation_frames, cation_times, cation_energies = read_energies('energies_cation.dat')
    neutral_frames, neutral_times, neutral_energies = read_energies('energies_neutral.dat')
    
    # Check that the number of data points match
    n_cation = len(cation_energies)
    n_neutral = len(neutral_energies)
    
    if n_cation != n_neutral:
        print(f"Warning: Number of data points differ: cation={n_cation}, neutral={n_neutral}")
        n_points = min(n_cation, n_neutral)
    else:
        n_points = n_cation
    
    # Verify time alignment
    time_diff = np.abs(cation_times[:n_points] - neutral_times[:n_points])
    if np.max(time_diff) > 1e-6:
        print("Warning: Time values do not match between files!")
        print(f"Max time difference: {np.max(time_diff)} fs")
    
    # Calculate energy difference (Cation - Neutral)
    dE_au = cation_energies[:n_points] - neutral_energies[:n_points]
    dE_eV = dE_au * HARTREE_TO_EV
    
    # Write output file
    output_file = 'energy_diff_elstat.dat'
    with open(output_file, 'w') as f:
        f.write("# Energy Difference (Cation - Neutral)\n")
        f.write("# Frame  Time(fs)        dE(a.u.)           dE(eV)\n")
        for i in range(n_points):
            frame = cation_frames[i]
            time = cation_times[i]
            f.write(f"{frame:5d}  {time:12.3f}  {dE_au[i]:18.10f}  {dE_eV[i]:12.6f}\n")
    
    print(f"Energy difference calculated for {n_points} frames")
    print(f"Output written to: {output_file}")
    print(f"Mean dE: {np.mean(dE_eV):.6f} eV")
    print(f"Std dE:  {np.std(dE_eV):.6f} eV")

if __name__ == '__main__':
    main()
