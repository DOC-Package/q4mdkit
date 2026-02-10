#!/usr/bin/env python
"""Calculate energy difference between cation and neutral states."""

import numpy as np

def calc_energy_diff(cation_file, neutral_file, output_file):
    """Calculate E_cation - E_neutral and save to file."""
    
    # Load data (skip 2 header lines)
    cation_data = np.loadtxt(cation_file, skiprows=2)
    neutral_data = np.loadtxt(neutral_file, skiprows=2)
    
    # Extract columns: Frame, Time(fs), Energy(a.u.)
    frames = cation_data[:, 0].astype(int)
    times = cation_data[:, 1]
    energy_c = cation_data[:, 2]
    energy_n = neutral_data[:, 2]
    
    # Calculate difference (cation - neutral)
    energy_diff = energy_c - energy_n
    
    # Convert to eV (1 a.u. = 27.211386 eV)
    hartree_to_ev = 27.211386
    energy_diff_ev = energy_diff * hartree_to_ev
    
    # Save to file
    with open(output_file, 'w') as f:
        f.write("# Energy Difference (Cation - Neutral)\n")
        f.write("# Frame  Time(fs)        dE(a.u.)           dE(eV)\n")
        for i in range(len(frames)):
            f.write(f"{frames[i]:5d}  {times[i]:12.3f}  {energy_diff[i]:18.10f}  {energy_diff_ev[i]:12.6f}\n")
    
    print(f"Output: {output_file}")
    print(f"  Number of frames: {len(frames)}")
    print(f"  Mean dE: {np.mean(energy_diff_ev):.6f} eV")
    print(f"  Std dE:  {np.std(energy_diff_ev):.6f} eV")

if __name__ == "__main__":
    # 790 series
    calc_energy_diff("energies_c790.dat", "energies_n790.dat", "energy_diff_790.dat")
    #print()
    
    # 1500 series
    #calc_energy_diff("energies_c1500.dat", "energies_n1500.dat", "energy_diff_1500.dat")
