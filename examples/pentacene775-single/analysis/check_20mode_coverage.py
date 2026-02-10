#!/usr/bin/env python3
"""
Check how well 20 selected modes can cover the target residual frequencies.
"""
import numpy as np

# Load mode information
with open('g_k_coefficients_below_2000cm.dat') as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith('#') and 'R²' not in l]
    data = [l.split() for l in lines if len(l.split()) >= 3]

mode_indices = np.array([int(d[0]) for d in data])
mode_freqs = np.array([float(d[1]) for d in data])

# Filter for 400-1800 cm⁻¹ range
fit_mask = (mode_freqs >= 400) & (mode_freqs <= 1800)
mode_indices_all = mode_indices[fit_mask]
mode_freqs_all = mode_freqs[fit_mask]

print("="*70)
print("Coverage Analysis: 20 Modes Selection")
print("="*70)

# Target residual peak frequencies
target_freqs = [524, 800, 1092, 1231, 1562]
tolerance = 30  # cm⁻¹

print(f"\nTarget residual peaks: {target_freqs} cm⁻¹")
print(f"Tolerance: ±{tolerance} cm⁻¹")

# Strategy 1: Use mode_pairs_for_residuals.dat to find most frequent modes
print("\n" + "="*70)
print("Strategy 1: Select based on frequency in existing pairs")
print("="*70)

# Read the existing mode pairs data to find most important modes
mode_counts = {}
try:
    with open('mode_pairs_for_residuals.dat', 'r') as f:
        for line in f:
            if line.startswith('#') or 'Mode' in line or '=' in line or 'Top' in line:
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    k = int(parts[0])
                    l = int(parts[1])
                    mode_counts[k] = mode_counts.get(k, 0) + 1
                    mode_counts[l] = mode_counts.get(l, 0) + 1
                except:
                    continue
    
    if mode_counts:
        # Get top 20 modes by frequency
        sorted_modes = sorted(mode_counts.items(), key=lambda x: x[1], reverse=True)
        top_20_mode_indices = [m[0] for m in sorted_modes[:20]]
        
        print(f"\nTop 20 most frequent modes in existing pairs:")
        for i, (mode_idx, count) in enumerate(sorted_modes[:20], 1):
            # Find frequency
            freq_idx = np.where(mode_indices_all == mode_idx)[0]
            if len(freq_idx) > 0:
                freq = mode_freqs_all[freq_idx[0]]
                print(f"  {i:2d}. Mode {mode_idx:2d}: {freq:6.1f} cm⁻¹ (appears {count:2d} times)")
    else:
        raise ValueError("No pairs found")
except Exception as e:
    print(f"Could not read mode_pairs_for_residuals.dat: {e}")
    print("Using fallback: strategically selected 20 modes")
    # Fallback: select modes to cover wide frequency range
    # Target: 524, 800, 1092, 1231, 1562 cm⁻¹
    # Need modes around: 250-800 cm⁻¹ for differences, and higher for sums
    # Select every 3rd or 4th mode from full range
    step = len(mode_indices_all) // 20
    if step == 0:
        step = 1
    top_20_mode_indices = mode_indices_all[::step][:20].tolist()

# Get frequencies for selected 20 modes
selected_mode_freqs = []
for mode_idx in top_20_mode_indices:
    freq_idx = np.where(mode_indices_all == mode_idx)[0]
    if len(freq_idx) > 0:
        selected_mode_freqs.append(mode_freqs_all[freq_idx[0]])

selected_mode_freqs = np.array(selected_mode_freqs)

print(f"\n" + "="*70)
print(f"Coverage Analysis with {len(selected_mode_freqs)} modes")
print("="*70)

# Check all pairs
n_modes = len(selected_mode_freqs)
total_pairs = n_modes * (n_modes - 1) // 2

print(f"\nTotal possible pairs: {total_pairs}")
print(f"Modes range: {selected_mode_freqs.min():.1f} - {selected_mode_freqs.max():.1f} cm⁻¹")

# Find pairs for each target frequency
coverage = {}
for target in target_freqs:
    sum_pairs = []
    diff_pairs = []
    
    for i in range(n_modes):
        for j in range(i+1, n_modes):
            freq_i = selected_mode_freqs[i]
            freq_j = selected_mode_freqs[j]
            
            # Check sum
            if abs((freq_i + freq_j) - target) <= tolerance:
                sum_pairs.append((top_20_mode_indices[i], top_20_mode_indices[j], 
                                freq_i, freq_j, freq_i + freq_j))
            
            # Check difference
            diff = abs(freq_i - freq_j)
            if abs(diff - target) <= tolerance:
                diff_pairs.append((top_20_mode_indices[i], top_20_mode_indices[j],
                                 freq_i, freq_j, diff))
    
    coverage[target] = {'sum': sum_pairs, 'diff': diff_pairs}
    
    print(f"\n{target} cm⁻¹:")
    print(f"  Sum combinations (ω_k + ω_l): {len(sum_pairs)}")
    print(f"  Diff combinations (|ω_k - ω_l|): {len(diff_pairs)}")
    
    # Show a few examples
    if sum_pairs:
        print(f"  Examples (sum):")
        for k, l, fk, fl, total in sum_pairs[:3]:
            print(f"    Mode {k} ({fk:.1f}) + Mode {l} ({fl:.1f}) = {total:.1f} cm⁻¹")
    
    if diff_pairs:
        print(f"  Examples (diff):")
        for k, l, fk, fl, total in diff_pairs[:3]:
            print(f"    Mode {k} ({fk:.1f}) - Mode {l} ({fl:.1f}) = {total:.1f} cm⁻¹")

# Summary
print("\n" + "="*70)
print("Summary")
print("="*70)

total_coverage = sum(len(coverage[t]['sum']) + len(coverage[t]['diff']) for t in target_freqs)
print(f"\nTotal pairs covering target frequencies: {total_coverage}")
print(f"Out of {total_pairs} possible pairs: {total_coverage/total_pairs*100:.1f}%")

print(f"\nBreakdown by target:")
for target in target_freqs:
    n_sum = len(coverage[target]['sum'])
    n_diff = len(coverage[target]['diff'])
    print(f"  {target:4d} cm⁻¹: {n_sum:2d} sum + {n_diff:2d} diff = {n_sum+n_diff:2d} pairs")

# Compare with 53 modes (420-1300 cm⁻¹)
fit_mask_53 = (mode_freqs_all >= 420) & (mode_freqs_all <= 1300)
n_modes_53 = fit_mask_53.sum()
total_pairs_53 = n_modes_53 * (n_modes_53 - 1) // 2

print(f"\n" + "="*70)
print("Comparison")
print("="*70)
print(f"\n20 modes (selected):  {total_pairs} pairs, {total_coverage} covering targets ({total_coverage/total_pairs*100:.1f}%)")
print(f"53 modes (420-1300):  {total_pairs_53} pairs, 232 covering targets (16.9%)")
print(f"\nEfficiency (coverage/pairs):")
print(f"  20 modes: {total_coverage/total_pairs:.3f}")
print(f"  53 modes: {232/total_pairs_53:.3f}")

# Calculate 3rd and 4th order terms to compute
n_3rd = 2 * total_pairs  # Two types of 3rd order
n_4th = total_pairs

print(f"\n" + "="*70)
print("Computational Cost")
print("="*70)
print(f"\n3rd order correlations: {n_3rd} terms")
print(f"4th order correlations: {n_4th} terms")
print(f"Total: {n_3rd + n_4th} correlation functions to compute")
print(f"\nParameters to fit: {total_pairs} (d_kl)")
