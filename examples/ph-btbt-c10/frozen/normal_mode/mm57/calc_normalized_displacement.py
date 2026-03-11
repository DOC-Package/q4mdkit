#!/usr/bin/env python3
"""
Calculate the normalized displacement from mode coordinates.

The normalized displacement is defined as:
    q̃_k(t) = Δq_k(t) / sqrt(<Δq_k²>_t)
where:
    Δq_k(t) = q_k(t) - <q_k>_t

This script reads mode_coords.txt and outputs mode_coords_normalized.txt
"""

import numpy as np

def main():
    # Read mode coordinates
    input_file = "mode_coords.txt"
    output_file = "mode_coords_normalized.txt"
    
    print(f"Reading {input_file}...")
    
    # Load data, skipping header lines (lines starting with #)
    data = np.loadtxt(input_file, comments='#')
    
    # First column is frame index, remaining columns are mode coordinates
    frame_indices = data[:, 0].astype(int)
    mode_coords = data[:, 1:]  # q_k(t) for all modes
    
    n_frames, n_modes = mode_coords.shape
    print(f"Number of frames: {n_frames}")
    print(f"Number of modes: {n_modes}")
    
    # Calculate the mean for each mode: <q_k>_t (time average)
    mean_q = np.mean(mode_coords, axis=0)  # shape: (n_modes,)
    
    # Calculate the deviation: Δq_k(t) = q_k(t) - <q_k>_t
    delta_q = mode_coords - mean_q  # shape: (n_frames, n_modes)
    
    # Calculate the variance: <Δq_k²>_t (time average of squared deviation)
    variance_q = np.mean(delta_q**2, axis=0)  # shape: (n_modes,)
    std_q = np.sqrt(variance_q)  # standard deviation
    
    # Calculate the normalized displacement: q̃_k(t) = Δq_k(t) / sqrt(<Δq_k²>_t)
    # Note: Be careful with division by zero for modes with very small variance
    with np.errstate(divide='ignore', invalid='ignore'):
        normalized_disp = delta_q / std_q  # shape: (n_frames, n_modes)
        # Replace inf/nan with 0 for modes with zero variance
        normalized_disp = np.nan_to_num(normalized_disp, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Print statistics
    print("\nStatistics:")
    print(f"Mean of q_k: min={mean_q.min():.6e}, max={mean_q.max():.6e}")
    print(f"Variance <Δq_k²>: min={variance_q.min():.6e}, max={variance_q.max():.6e}")
    print(f"Normalized displacement: min={normalized_disp.min():.6e}, max={normalized_disp.max():.6e}")
    
    # Save to file
    print(f"\nSaving to {output_file}...")
    
    # Create header
    header = f"""Normalized displacement from {input_file}
Definition: q̃_k(t) = Δq_k(t) / sqrt(<Δq_k²>_t), where Δq_k(t) = q_k(t) - <q_k>_t
Frames: {n_frames}, Modes: {n_modes}
Mode indices: 7 to 108 (1-indexed)
Columns: frame_index, q̃_7, q̃_8, ..., q̃_3N"""
    
    # Combine frame indices with normalized displacement
    output_data = np.column_stack([frame_indices, normalized_disp])
    
    # Save with scientific notation
    np.savetxt(output_file, output_data, 
               header=header, 
               fmt=['%8d'] + ['%15.8e'] * n_modes,
               comments='# ')
    
    print("Done!")
    
    # Also save statistics to a separate file
    stats_file = "mode_coords_normalized_stats.txt"
    print(f"\nSaving statistics to {stats_file}...")
    
    with open(stats_file, 'w') as f:
        f.write("# Statistics for normalized displacement calculation\n")
        f.write(f"# Input file: {input_file}\n")
        f.write(f"# Number of frames: {n_frames}\n")
        f.write(f"# Number of modes: {n_modes}\n")
        f.write("#\n")
        f.write("# Mode_index  Mean_q_k  Variance_<Δq_k²>  Std_dev\n")
        
        for i in range(n_modes):
            mode_idx = i + 7  # Mode indices start at 7
            f.write(f"{mode_idx:4d}  {mean_q[i]:15.8e}  {variance_q[i]:15.8e}  {np.sqrt(variance_q[i]):15.8e}\n")
    
    print("Done!")

if __name__ == "__main__":
    main()
