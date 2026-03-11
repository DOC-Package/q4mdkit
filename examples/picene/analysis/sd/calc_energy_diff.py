#!/usr/bin/env python
"""
Calculate energy difference between cation and neutral states.

Reads:
- qm_energy_sampled.dat: Cation energies
- energies.dat: Neutral energies

Outputs:
- energy_diff.dat: Energy difference (Cation - Neutral) in a.u. and eV
"""

import numpy as np

# Conversion factor from Hartree to eV
HARTREE_TO_EV = 27.211386245988

def read_neutral_energies(filename):
    """Read neutral energies from qm_energy_sampled.dat"""
    data = np.loadtxt(filename, comments='#')
    # columns: call_index, time[fs], energy[a.u.]
    times = data[:, 1]
    energies = data[:, 2]
    return times, energies

def read_cation_energies(filename):
    """Read cation energies from energies.dat"""
    data = np.loadtxt(filename, comments='#')
    # columns: Frame, Time(fs), Energy(a.u.)
    frames = data[:, 0].astype(int)
    times = data[:, 1]
    energies = data[:, 2]
    return frames, times, energies

def main():
    # Read input files
    #neutral_times, neutral_energies = read_neutral_energies('qm_energy_sampled.dat')
    #cation_frames, cation_times, cation_energies = read_cation_energies('energies.dat')
    neutral_times, neutral_energies = read_neutral_energies('energies_cation_wopc.dat')
    cation_frames, cation_times, cation_energies = read_cation_energies('energies_neutral_wopc.dat')
    
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
    
    # Write output file in the same format as energy_diff_all_pentacene.dat
    output_file = 'energy_diff_wopc.dat'
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
