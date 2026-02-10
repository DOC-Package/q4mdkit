#!/usr/bin/env python3
"""
各三次相関項（Type1〜Type4）がどのようなピーク位置を持つか解析する

3rd order correlations:
  Type 1: FT[<q_k(t)q_l(t)q_k(0)>] → ω_k + ω_l, |ω_k - ω_l| にピーク
  Type 2: FT[<q_k(t)q_l(t)q_l(0)>] → ω_k + ω_l, |ω_k - ω_l| にピーク  
  Type 3: FT[<q_k(t)q_k(0)q_l(0)>] → ω_k にピーク
  Type 4: FT[<q_l(t)q_l(0)q_k(0)>] → ω_l にピーク
"""
import numpy as np
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("Analysis: Peak Positions of Each 3rd Order Correlation Type")
print("="*70)

# ============================================================================
# Load data
# ============================================================================
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff = energy_diff_data[:, 4]

# Load mode information
with open('g_k_coefficients_below_2000cm.dat') as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith('#') and 'R²' not in l]
    data = [l.split() for l in lines if len(l.split()) >= 3]

mode_indices = np.array([int(d[0]) for d in data])
mode_freqs = np.array([float(d[1]) for d in data])

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords.txt')
q_k_all = q_k_data[:, 1:]
mode_coords = q_k_all[:, mode_indices - 7]

N_frames = min(len(energy_diff), len(mode_coords))
energy_diff = energy_diff[:N_frames]
mode_coords = mode_coords[:N_frames, :]
N_modes = len(mode_freqs)

print(f"  Frames: {N_frames}")
print(f"  Modes: {N_modes}")

# ============================================================================
# Helper functions
# ============================================================================
def compute_spectral_from_corr(corr, freq):
    """Compute spectral density from pre-computed correlation"""
    t = np.arange(len(corr)) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = corr * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

# Frequency grid - wider range to see all peaks
freq_grid = np.linspace(0, 2000, 2001)
print(f"\nFrequency grid: {len(freq_grid)} points (0-2000 cm⁻¹)")

# ============================================================================
# Select representative mode pairs to analyze
# ============================================================================
print("\n" + "="*70)
print("Selecting representative mode pairs for analysis")
print("="*70)

# Find modes close to specific frequencies for clear demonstration
def find_mode_near_freq(target, mode_freqs, tolerance=30):
    """Find mode index closest to target frequency"""
    diffs = np.abs(mode_freqs - target)
    if np.min(diffs) < tolerance:
        return np.argmin(diffs)
    return None

# Select pairs that give interesting sum/difference frequencies
test_pairs = []

# Pair 1: Two low-frequency modes (sum around 500-600)
idx1 = find_mode_near_freq(250, mode_freqs, 50)
idx2 = find_mode_near_freq(270, mode_freqs, 50)
if idx1 is not None and idx2 is not None and idx1 != idx2:
    test_pairs.append((idx1, idx2, 'Low+Low'))

# Pair 2: Low + medium frequency modes (sum around 800)  
idx1 = find_mode_near_freq(300, mode_freqs, 50)
idx2 = find_mode_near_freq(500, mode_freqs, 50)
if idx1 is not None and idx2 is not None:
    test_pairs.append((idx1, idx2, 'Low+Mid'))

# Pair 3: Two modes with sum around 1092
idx1 = find_mode_near_freq(500, mode_freqs, 50)
idx2 = find_mode_near_freq(600, mode_freqs, 50)
if idx1 is not None and idx2 is not None:
    test_pairs.append((idx1, idx2, 'Mid+Mid→1092'))

# Pair 4: Two modes with difference around 524
idx1 = find_mode_near_freq(1000, mode_freqs, 50)
idx2 = find_mode_near_freq(476, mode_freqs, 50)  
if idx1 is not None and idx2 is not None:
    test_pairs.append((idx1, idx2, 'Diff→524'))

# Pair 5: Two similar frequency modes (small difference)
idx1 = find_mode_near_freq(500, mode_freqs, 30)
idx2 = find_mode_near_freq(530, mode_freqs, 50)
if idx1 is not None and idx2 is not None and idx1 != idx2:
    test_pairs.append((idx1, idx2, 'Similar freq'))

# Pair 6: High frequency modes
idx1 = find_mode_near_freq(1200, mode_freqs, 50)
idx2 = find_mode_near_freq(1400, mode_freqs, 50)
if idx1 is not None and idx2 is not None:
    test_pairs.append((idx1, idx2, 'High+High'))

print(f"\nSelected {len(test_pairs)} test pairs:")
for i, (idx1, idx2, name) in enumerate(test_pairs):
    freq1, freq2 = mode_freqs[idx1], mode_freqs[idx2]
    sum_f = freq1 + freq2
    diff_f = abs(freq1 - freq2)
    print(f"  {i+1}. {name}: Mode {mode_indices[idx1]} ({freq1:.1f}) × Mode {mode_indices[idx2]} ({freq2:.1f})")
    print(f"      Sum = {sum_f:.1f} cm⁻¹, Diff = {diff_f:.1f} cm⁻¹")

# ============================================================================
# Compute 3rd order correlations for each pair
# ============================================================================
print("\n" + "="*70)
print("Computing 3rd order correlation spectral densities for each type")
print("="*70)

t_arr = np.arange(max_lag)
window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag))

results = []

for pair_num, (i, j, pair_name) in enumerate(test_pairs):
    print(f"\n--- Pair {pair_num+1}: {pair_name} ---")
    
    freq_i = mode_freqs[i]
    freq_j = mode_freqs[j]
    q_k = mode_coords[:, i]
    q_l = mode_coords[:, j]
    
    # Center
    q_k_c = q_k - np.mean(q_k)
    q_l_c = q_l - np.mean(q_l)
    
    # Product at same time
    product_kl = q_k_c * q_l_c
    
    print(f"  Mode k: {mode_indices[i]} ({freq_i:.1f} cm⁻¹)")
    print(f"  Mode l: {mode_indices[j]} ({freq_j:.1f} cm⁻¹)")
    print(f"  Expected peaks:")
    print(f"    Sum frequency: {freq_i + freq_j:.1f} cm⁻¹")
    print(f"    Diff frequency: {abs(freq_i - freq_j):.1f} cm⁻¹")
    
    # Type 1: <q_k(t)q_l(t)q_k(0)>
    corr_type1 = np.zeros(max_lag)
    for lag in range(max_lag):
        corr_type1[lag] = np.mean(product_kl[lag:] * q_k_c[:N_frames-lag])
    corr_type1 *= window
    J_type1 = compute_spectral_from_corr(corr_type1, freq_grid)
    
    # Type 2: <q_k(t)q_l(t)q_l(0)>
    corr_type2 = np.zeros(max_lag)
    for lag in range(max_lag):
        corr_type2[lag] = np.mean(product_kl[lag:] * q_l_c[:N_frames-lag])
    corr_type2 *= window
    J_type2 = compute_spectral_from_corr(corr_type2, freq_grid)
    
    # Type 3: <q_k(t)q_k(0)q_l(0)>
    corr_type3 = np.zeros(max_lag)
    for lag in range(max_lag):
        corr_type3[lag] = np.mean(q_k_c[lag:] * product_kl[:N_frames-lag])
    corr_type3 *= window
    J_type3 = compute_spectral_from_corr(corr_type3, freq_grid)
    
    # Type 4: <q_l(t)q_l(0)q_k(0)>
    corr_type4 = np.zeros(max_lag)
    for lag in range(max_lag):
        corr_type4[lag] = np.mean(q_l_c[lag:] * product_kl[:N_frames-lag])
    corr_type4 *= window
    J_type4 = compute_spectral_from_corr(corr_type4, freq_grid)
    
    # Find peaks
    def find_peaks(J, freq_grid, threshold_ratio=0.1):
        """Find peak positions"""
        from scipy.signal import find_peaks as sp_find_peaks
        abs_J = np.abs(J)
        threshold = threshold_ratio * np.max(abs_J) if np.max(abs_J) > 0 else 0
        peaks, props = sp_find_peaks(abs_J, height=threshold)
        return freq_grid[peaks], J[peaks]
    
    # Detect and report peaks
    print(f"\n  Peak analysis:")
    
    peak_data = {}
    for type_name, J in [('Type1', J_type1), ('Type2', J_type2), ('Type3', J_type3), ('Type4', J_type4)]:
        peak_freqs, peak_vals = find_peaks(J, freq_grid, 0.15)
        peak_data[type_name] = (peak_freqs, peak_vals)
        if len(peak_freqs) > 0:
            sorted_idx = np.argsort(np.abs(peak_vals))[::-1]
            top_peaks = peak_freqs[sorted_idx][:5]
            print(f"    {type_name}: peaks at {', '.join([f'{p:.0f}' for p in top_peaks])} cm⁻¹")
        else:
            print(f"    {type_name}: no significant peaks detected")
    
    results.append({
        'name': pair_name,
        'i': i, 'j': j,
        'mode_i': mode_indices[i], 'mode_j': mode_indices[j],
        'freq_i': freq_i, 'freq_j': freq_j,
        'J_type1': J_type1, 'J_type2': J_type2,
        'J_type3': J_type3, 'J_type4': J_type4,
        'peak_data': peak_data
    })

# ============================================================================
# Generate detailed plots
# ============================================================================
print("\n" + "="*70)
print("Generating plots")
print("="*70)

# Create comprehensive plot for each pair
for pair_num, res in enumerate(results):
    freq_i = res['freq_i']
    freq_j = res['freq_j']
    sum_freq = freq_i + freq_j
    diff_freq = abs(freq_i - freq_j)
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    
    # Title
    fig.suptitle(f'3rd Order Correlation Types: Mode {res["mode_i"]} ({freq_i:.1f} cm⁻¹) × Mode {res["mode_j"]} ({freq_j:.1f} cm⁻¹)\n'
                 f'Sum = {sum_freq:.1f} cm⁻¹, Diff = {diff_freq:.1f} cm⁻¹',
                 fontsize=14, fontweight='bold')
    
    # Plot each type
    types = [
        ('Type 1: ⟨q_k(t)q_l(t)q_k(0)⟩', res['J_type1'], 'blue'),
        ('Type 2: ⟨q_k(t)q_l(t)q_l(0)⟩', res['J_type2'], 'green'),
        ('Type 3: ⟨q_k(t)q_k(0)q_l(0)⟩', res['J_type3'], 'red'),
        ('Type 4: ⟨q_l(t)q_l(0)q_k(0)⟩', res['J_type4'], 'purple'),
    ]
    
    ax_positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    
    for (title, J, color), (row, col) in zip(types, ax_positions):
        ax = axes[row, col]
        ax.plot(freq_grid, J, color=color, lw=1.5)
        ax.axhline(0, color='k', ls='--', lw=0.5)
        
        # Mark expected peak positions
        ax.axvline(freq_i, color='gray', ls=':', alpha=0.5, label=f'ω_k={freq_i:.0f}')
        ax.axvline(freq_j, color='gray', ls='--', alpha=0.5, label=f'ω_l={freq_j:.0f}')
        ax.axvline(sum_freq, color='orange', ls='-', alpha=0.7, label=f'Sum={sum_freq:.0f}')
        if diff_freq > 10:  # Only show if diff is significant
            ax.axvline(diff_freq, color='cyan', ls='-', alpha=0.7, label=f'Diff={diff_freq:.0f}')
        
        ax.set_xlabel('Frequency (cm⁻¹)')
        ax.set_ylabel('J(ν)')
        ax.set_title(title, fontweight='bold')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, min(2000, sum_freq + 200)])
    
    # Combined plot
    ax5 = axes[0, 2]
    ax5.plot(freq_grid, res['J_type1'], 'b-', lw=1, alpha=0.7, label='Type 1')
    ax5.plot(freq_grid, res['J_type2'], 'g-', lw=1, alpha=0.7, label='Type 2')
    ax5.plot(freq_grid, res['J_type3'], 'r-', lw=1, alpha=0.7, label='Type 3')
    ax5.plot(freq_grid, res['J_type4'], 'purple', lw=1, alpha=0.7, label='Type 4')
    ax5.axhline(0, color='k', ls='--', lw=0.5)
    ax5.axvline(sum_freq, color='orange', ls='-', alpha=0.7, lw=2)
    ax5.axvline(diff_freq, color='cyan', ls='-', alpha=0.7, lw=2)
    ax5.set_xlabel('Frequency (cm⁻¹)')
    ax5.set_ylabel('J(ν)')
    ax5.set_title('All Types Comparison', fontweight='bold')
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim([0, min(2000, sum_freq + 200)])
    
    # Combination used in v4 (κ_k*(Type1+Type4) + κ_l*(Type2+Type3))
    ax6 = axes[1, 2]
    J_combined_ki = res['J_type1'] + res['J_type4']  # multiplied by κ_k
    J_combined_kj = res['J_type2'] + res['J_type3']  # multiplied by κ_l
    
    ax6.plot(freq_grid, J_combined_ki, 'b-', lw=1.5, alpha=0.8, label='Type1+Type4 (×κ_k)')
    ax6.plot(freq_grid, J_combined_kj, 'r-', lw=1.5, alpha=0.8, label='Type2+Type3 (×κ_l)')
    ax6.axhline(0, color='k', ls='--', lw=0.5)
    ax6.axvline(sum_freq, color='orange', ls='-', alpha=0.7, lw=2, label=f'Sum={sum_freq:.0f}')
    ax6.axvline(diff_freq, color='cyan', ls='-', alpha=0.7, lw=2, label=f'Diff={diff_freq:.0f}')
    ax6.axvline(freq_i, color='gray', ls=':', alpha=0.5, label=f'ω_k={freq_i:.0f}')
    ax6.axvline(freq_j, color='gray', ls='--', alpha=0.5, label=f'ω_l={freq_j:.0f}')
    ax6.set_xlabel('Frequency (cm⁻¹)')
    ax6.set_ylabel('J(ν)')
    ax6.set_title('v4 Combinations', fontweight='bold')
    ax6.legend(fontsize=8, loc='best')
    ax6.grid(True, alpha=0.3)
    ax6.set_xlim([0, min(2000, sum_freq + 200)])
    
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    filename = f'3rd_order_peaks_pair{pair_num+1}_{res["name"].replace("+", "_").replace("→", "_").replace(" ", "_")}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Saved: {filename}")
    plt.close()

# ============================================================================
# Summary plot comparing all pairs
# ============================================================================
print("\nGenerating summary plot...")

n_pairs = len(results)
fig, axes = plt.subplots(n_pairs, 4, figsize=(18, 3*n_pairs))

for pair_num, res in enumerate(results):
    freq_i = res['freq_i']
    freq_j = res['freq_j']
    sum_freq = freq_i + freq_j
    diff_freq = abs(freq_i - freq_j)
    
    for type_idx, (J, color, title) in enumerate([
        (res['J_type1'], 'blue', 'Type 1'),
        (res['J_type2'], 'green', 'Type 2'),
        (res['J_type3'], 'red', 'Type 3'),
        (res['J_type4'], 'purple', 'Type 4')
    ]):
        ax = axes[pair_num, type_idx] if n_pairs > 1 else axes[type_idx]
        ax.plot(freq_grid, J, color=color, lw=1)
        ax.axhline(0, color='k', ls='--', lw=0.3)
        ax.axvline(freq_i, color='gray', ls=':', alpha=0.5)
        ax.axvline(freq_j, color='gray', ls='--', alpha=0.5)
        ax.axvline(sum_freq, color='orange', ls='-', alpha=0.7)
        ax.axvline(diff_freq, color='cyan', ls='-', alpha=0.7)
        ax.set_xlim([0, min(2000, sum_freq + 200)])
        
        if pair_num == 0:
            ax.set_title(title, fontsize=11, fontweight='bold')
        if type_idx == 0:
            ax.set_ylabel(f'{res["name"]}\nω_k={freq_i:.0f}, ω_l={freq_j:.0f}', fontsize=9)
        if pair_num == n_pairs - 1:
            ax.set_xlabel('Frequency (cm⁻¹)')

fig.suptitle('3rd Order Correlation Types: Peak Positions\n'
             'Gray lines: ω_k, ω_l  |  Orange: Sum (ω_k+ω_l)  |  Cyan: Diff (|ω_k-ω_l|)',
             fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('3rd_order_peaks_summary.png', dpi=150, bbox_inches='tight')
print("Saved: 3rd_order_peaks_summary.png")
plt.close()

# ============================================================================
# Theoretical summary
# ============================================================================
print("\n" + "="*70)
print("THEORETICAL SUMMARY: Peak Positions of 3rd Order Correlation Types")
print("="*70)

print("""
For a pair of harmonic modes q_k (frequency ω_k) and q_l (frequency ω_l):

Type 1: ⟨q_k(t)q_l(t)q_k(0)⟩
  - Time evolution: q_k(t)q_l(t) oscillates at ω_k+ω_l and |ω_k-ω_l|
  - Correlation with q_k(0) selects components that match ω_k
  → Expected peaks: ω_k+ω_l, |ω_k-ω_l| (combination/difference frequencies)
  → Additional peak possible at ω_k (direct correlation)

Type 2: ⟨q_k(t)q_l(t)q_l(0)⟩  
  - Time evolution: q_k(t)q_l(t) oscillates at ω_k+ω_l and |ω_k-ω_l|
  - Correlation with q_l(0) selects components that match ω_l
  → Expected peaks: ω_k+ω_l, |ω_k-ω_l| (combination/difference frequencies)
  → Additional peak possible at ω_l (direct correlation)

Type 3: ⟨q_k(t)q_k(0)q_l(0)⟩
  - q_k(t) correlates with q_k(0) → gives peak at ω_k
  - q_l(0) is static (t=0), acts as amplitude modulation
  → Primary peak at ω_k (mode k's fundamental frequency)
  → Modified by cross-correlation with q_l(0)

Type 4: ⟨q_l(t)q_l(0)q_k(0)⟩
  - q_l(t) correlates with q_l(0) → gives peak at ω_l
  - q_k(0) is static (t=0), acts as amplitude modulation
  → Primary peak at ω_l (mode l's fundamental frequency)
  → Modified by cross-correlation with q_k(0)

In the v4 fitting formulation:
  J^(3)_kl = d_kl × [κ_k(Type1+Type4) + κ_l(Type2+Type3)]

  - κ_k(Type1+Type4): Peaks at combination freqs (ω_k+ω_l, |ω_k-ω_l|) + ω_l
  - κ_l(Type2+Type3): Peaks at combination freqs (ω_k+ω_l, |ω_k-ω_l|) + ω_k

Key insight: 3rd order terms can generate spectral density at frequencies
different from the individual mode frequencies, allowing them to fill in
gaps not covered by the linear (quadratic coupling) terms.
""")

print("\n" + "="*70)
print("Analysis complete!")
print("="*70)
