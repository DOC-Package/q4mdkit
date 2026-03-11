#!/usr/bin/env python
"""
Subtract energy_diff_elstat.dat from energy_diff.dat
"""

import numpy as np

def read_energy_diff(filename):
    """Read energy difference file"""
    frames = []
    times = []
    dE_au = []
    dE_eV = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            try:
                frames.append(int(cols[0]))
                times.append(float(cols[1]))
                dE_au.append(float(cols[2]))
                dE_eV.append(float(cols[3]))
            except (ValueError, IndexError):
                continue
    return np.array(frames), np.array(times), np.array(dE_au), np.array(dE_eV)

def main():
    # Read input files
    frames1, times1, dE_au1, dE_eV1 = read_energy_diff('energy_diff.dat')
    frames2, times2, dE_au2, dE_eV2 = read_energy_diff('energy_diff_elstat.dat')
    
    n_points = min(len(frames1), len(frames2))
    
    # Subtract: energy_diff - energy_diff_elstat
    dE_au_diff = dE_au1[:n_points] - dE_au2[:n_points]
    dE_eV_diff = dE_eV1[:n_points] - dE_eV2[:n_points]
    
    # Write output
    output_file = 'energy_diff_minus_elstat.dat'
    with open(output_file, 'w') as f:
        f.write("# Energy Difference: (energy_diff - energy_diff_elstat)\n")
        f.write("# Frame  Time(fs)        dE(a.u.)           dE(eV)\n")
        for i in range(n_points):
            f.write(f"{frames1[i]:5d}  {times1[i]:12.3f}  {dE_au_diff[i]:18.10f}  {dE_eV_diff[i]:12.6f}\n")
    
    print(f"Output written to: {output_file}")
    print(f"Points: {n_points}")
    print(f"Mean dE: {np.nanmean(dE_eV_diff):.6f} eV")
    print(f"Std dE:  {np.nanstd(dE_eV_diff):.6f} eV")

if __name__ == '__main__':
    main()
