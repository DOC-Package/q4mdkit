#!/usr/bin/env python
"""
Subtract mean-subtracted energy_diff_elstat from mean-subtracted energy_diff
(energy_diff - mean(energy_diff)) - (energy_diff_elstat - mean(energy_diff_elstat))
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
    
    # Get valid data for mean calculation
    dE_eV1_valid = dE_eV1[:n_points]
    dE_eV2_valid = dE_eV2[:n_points]
    dE_au1_valid = dE_au1[:n_points]
    dE_au2_valid = dE_au2[:n_points]
    
    # Calculate means
    mean1_eV = np.nanmean(dE_eV1_valid)
    mean2_eV = np.nanmean(dE_eV2_valid)
    mean1_au = np.nanmean(dE_au1_valid)
    mean2_au = np.nanmean(dE_au2_valid)
    
    print(f"Mean energy_diff: {mean1_eV:.6f} eV")
    print(f"Mean energy_diff_elstat: {mean2_eV:.6f} eV")
    
    # Subtract means then calculate difference
    # (dE1 - mean1) - (dE2 - mean2)
    dE_au_diff = (dE_au1_valid - mean1_au) - (dE_au2_valid - mean2_au)
    dE_eV_diff = (dE_eV1_valid - mean1_eV) - (dE_eV2_valid - mean2_eV)
    
    # Write output
    output_file = 'energy_diff_mean_subtracted.dat'
    with open(output_file, 'w') as f:
        f.write("# (energy_diff - mean) - (energy_diff_elstat - mean)\n")
        f.write(f"# Mean energy_diff: {mean1_eV:.6f} eV, Mean energy_diff_elstat: {mean2_eV:.6f} eV\n")
        f.write("# Frame  Time(fs)        dE(a.u.)           dE(eV)\n")
        for i in range(n_points):
            f.write(f"{frames1[i]:5d}  {times1[i]:12.3f}  {dE_au_diff[i]:18.10f}  {dE_eV_diff[i]:12.6f}\n")
    
    print(f"\nOutput written to: {output_file}")
    print(f"Points: {n_points}")
    print(f"Mean of result: {np.nanmean(dE_eV_diff):.6f} eV")
    print(f"Std of result:  {np.nanstd(dE_eV_diff):.6f} eV")

if __name__ == '__main__':
    main()
