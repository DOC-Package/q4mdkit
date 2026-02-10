#!/usr/bin/env python3
"""
Third-order correlation function analysis:
Compute C3_type1(t) = <q_k(t)q_l(t)q_j(0)> and C3_type2(t) = <q_j(t)q_k(0)q_l(0)>
and analyze their Fourier transforms.
Compare with 4th order correlations for triplets with high-frequency peaks.
"""
import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# ============================================================================
# Setup logging to both terminal and file
# ============================================================================
class Logger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w')
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = Logger('third_order_correlation.log')

# Constants
c_cm_fs = 2.99792458e-5  # cm/fs
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("Third-Order Correlation Function Analysis")
print("="*70)

# Load data
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff = energy_diff_data[:, 4]

# Load mode information
mode_indices_all = []
mode_freqs_all = []
with open('mode_coords_new_stats.txt') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            mode_indices_all.append(int(parts[0]))
            mode_freqs_all.append(float(parts[1]))

mode_indices = np.array(mode_indices_all)
mode_freqs = np.array(mode_freqs_all)

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords_new.txt')
q_k_all = q_k_data[:, 1:]
mode_coords = q_k_all[:, mode_indices - 7]

N_frames = min(len(energy_diff), len(mode_coords))
energy_diff = energy_diff[:N_frames]
mode_coords = mode_coords[:N_frames, :]
N_modes = len(mode_freqs)

print(f"  Frames: {N_frames}")
print(f"  Total modes: {N_modes}")

# Filter modes: 100-1600 cm⁻¹
freq_min, freq_max = 100, 1600
mask = (mode_freqs >= freq_min) & (mode_freqs <= freq_max)
idx_selected = np.where(mask)[0]
n_selected = len(idx_selected)
print(f"  Modes in {freq_min}-{freq_max} cm⁻¹: {n_selected}")

# ============================================================================
# Find top 100 mode triplets by |<q_k(0)q_l(0)q_j(0)>| (three-body correlation)
# ============================================================================
print("\n" + "="*70)
print("Finding top 100 mode triplets by three-body correlation")
print("="*70)

# Compute three-body correlation for selected modes
print("\nComputing three-body correlations <q_k q_l q_j>...")
q_selected = mode_coords[:, idx_selected]
q_mean = np.mean(q_selected, axis=0)
q_centered = q_selected - q_mean

# Find all triplets (k, l, j) with k <= l <= j
triplet_values = []
total_triplets = n_selected * (n_selected + 1) * (n_selected + 2) // 6
print(f"  Total triplets to compute: {total_triplets}")

count = 0
for i in range(n_selected):
    for j in range(i, n_selected):
        for k in range(j, n_selected):
            # <q_i q_j q_k>
            corr3 = np.mean(q_centered[:, i] * q_centered[:, j] * q_centered[:, k])
            triplet_values.append((i, j, k, np.abs(corr3), corr3))
            count += 1
    if (i + 1) % 10 == 0:
        print(f"    Progress: {i+1}/{n_selected} modes")

print(f"  Computed {count} triplets")

# Sort by absolute correlation value
triplet_values.sort(key=lambda x: x[3], reverse=True)
top_500_triplets = triplet_values[:500]

print(f"\nTop 20 mode triplets by |<q_k q_l q_j>|:")
for rank, (i, j, k, absval, val) in enumerate(top_500_triplets[:20], 1):
    idx_i = idx_selected[i]
    idx_j = idx_selected[j]
    idx_k = idx_selected[k]
    freq_i = mode_freqs[idx_i]
    freq_j = mode_freqs[idx_j]
    freq_k = mode_freqs[idx_k]
    mode_i = mode_indices[idx_i]
    mode_j = mode_indices[idx_j]
    mode_k = mode_indices[idx_k]
    print(f"  {rank:2d}. Mode {mode_i:2d}({freq_i:5.0f}) × {mode_j:2d}({freq_j:5.0f}) × {mode_k:2d}({freq_k:5.0f}): <q³> = {val:+.4e}")

# ============================================================================
# Compute third-order correlation functions
# ============================================================================
print("\n" + "="*70)
print("Computing third-order correlation functions")
print("="*70)

def compute_third_order_corr_type1(q_k, q_l, q_j):
    """
    Compute C3_type1(t) = <q_k(t) * q_l(t) * q_j(0)>
    This probes how the product q_k*q_l at time t correlates with q_j at time 0.
    """
    q_j_mean = np.mean(q_j)
    q_j_centered = q_j - q_j_mean
    
    # q_k(t) * q_l(t) product
    q_kl = q_k * q_l
    q_kl_mean = np.mean(q_kl)
    q_kl_centered = q_kl - q_kl_mean
    
    N = len(q_k)
    C3 = np.zeros(max_lag)
    
    for lag in range(max_lag):
        # <q_k(t)q_l(t)q_j(0)> = <q_kl(lag:) * q_j(:N-lag)>
        C3[lag] = np.mean(q_kl_centered[lag:] * q_j_centered[:N-lag])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return C3 #* window

def compute_third_order_corr_type2(q_j, q_k, q_l):
    """
    Compute C3_type2(t) = <q_j(t) * q_k(0) * q_l(0)>
    This probes how q_j at time t correlates with the product q_k*q_l at time 0.
    """
    q_kl = q_k * q_l
    q_kl_mean = np.mean(q_kl)
    q_kl_centered = q_kl - q_kl_mean
    
    q_j_mean = np.mean(q_j)
    q_j_centered = q_j - q_j_mean
    
    N = len(q_j)
    C3 = np.zeros(max_lag)
    
    for lag in range(max_lag):
        # <q_j(t)q_k(0)q_l(0)> = <q_j(lag:) * q_kl(:N-lag)>
        C3[lag] = np.mean(q_j_centered[lag:] * q_kl_centered[:N-lag])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return C3 #* window

def compute_fourth_order_corr(q_i, q_j, q_k, q_l):
    """
    Compute C4(t) = <q_i(t) * q_j(t) * q_k(0) * q_l(0)>
    This probes how the product q_i*q_j at time t correlates with q_k*q_l at time 0.
    """
    q_ij = q_i * q_j
    q_ij_mean = np.mean(q_ij)
    q_ij_centered = q_ij - q_ij_mean
    
    q_kl = q_k * q_l
    q_kl_mean = np.mean(q_kl)
    q_kl_centered = q_kl - q_kl_mean
    
    N = len(q_i)
    C4 = np.zeros(max_lag)
    
    for lag in range(max_lag):
        C4[lag] = np.mean(q_ij_centered[lag:] * q_kl_centered[:N-lag])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return C4 #* window

def compute_fourier_transform(C, freq_grid):
    """Compute Fourier transform of correlation function"""
    t = np.arange(max_lag) * dt
    
    J = np.zeros(len(freq_grid))
    for i, nu in enumerate(freq_grid):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    return J

def find_spectral_peaks(J, freq, threshold=0.1, min_separation=20):
    """Find peaks in spectral density"""
    max_J = np.max(np.abs(J))
    if max_J <= 0:
        return np.array([]), np.array([])
    
    # Find peaks in positive part
    peaks_pos_idx, _ = find_peaks(
        J, 
        height=threshold * max_J,
        prominence=threshold * max_J * 0.3,
        distance=int(min_separation / (freq[1] - freq[0]))
    )
    
    # Find peaks in negative part (inverted)
    peaks_neg_idx, _ = find_peaks(
        -J, 
        height=threshold * max_J,
        prominence=threshold * max_J * 0.3,
        distance=int(min_separation / (freq[1] - freq[0]))
    )
    
    # Combine and sort
    all_peaks_idx = np.concatenate([peaks_pos_idx, peaks_neg_idx])
    if len(all_peaks_idx) == 0:
        return np.array([]), np.array([])
    
    peak_freqs = freq[all_peaks_idx]
    peak_heights = J[all_peaks_idx]
    
    # Sort by absolute height
    sorted_idx = np.argsort(np.abs(peak_heights))[::-1]
    return peak_freqs[sorted_idx], peak_heights[sorted_idx]

# Frequency grid
freq_grid = np.linspace(0, 3000, 3001)  # Extended to capture sum frequencies
print(f"\nFrequency grid: 0-3000 cm⁻¹ ({len(freq_grid)} points)")

# ============================================================================
# Analyze third-order correlations for top 500 triplets
# ============================================================================
print("\nComputing third-order correlation functions for top 500 triplets...")

results_type1 = []  # <q_i(t)q_j(t)q_k(0)>
results_type2 = []  # <q_k(t)q_i(0)q_j(0)>

for triplet_idx, (i_local, j_local, k_local, absval, corr_val) in enumerate(top_500_triplets):
    if (triplet_idx + 1) % 50 == 0:
        print(f"  Processing triplet {triplet_idx + 1}/500...")
    
    idx_i = idx_selected[i_local]
    idx_j = idx_selected[j_local]
    idx_k = idx_selected[k_local]
    
    q_i = mode_coords[:, idx_i]
    q_j = mode_coords[:, idx_j]
    q_k = mode_coords[:, idx_k]
    
    freq_i = mode_freqs[idx_i]
    freq_j = mode_freqs[idx_j]
    freq_k = mode_freqs[idx_k]
    mode_i = mode_indices[idx_i]
    mode_j = mode_indices[idx_j]
    mode_k = mode_indices[idx_k]
    
    # Type 1: <q_i(t)q_j(t)q_k(0)> - product of two modes at t, third at 0
    C3_1 = compute_third_order_corr_type1(q_i, q_j, q_k)
    J3_1 = compute_fourier_transform(C3_1, freq_grid)
    peaks_1, heights_1 = find_spectral_peaks(J3_1, freq_grid)
    
    results_type1.append({
        'i': idx_i, 'j': idx_j, 'k': idx_k,
        'mode_i': mode_i, 'mode_j': mode_j, 'mode_k': mode_k,
        'freq_i': freq_i, 'freq_j': freq_j, 'freq_k': freq_k,
        'corr3': corr_val,
        'C3': C3_1,
        'J3': J3_1,
        'peaks': peaks_1,
        'heights': heights_1,
        'type': 'type1'
    })
    
    # Type 2: <q_k(t)q_i(0)q_j(0)> - one mode at t, product of two at 0
    C3_2 = compute_third_order_corr_type2(q_k, q_i, q_j)
    J3_2 = compute_fourier_transform(C3_2, freq_grid)
    peaks_2, heights_2 = find_spectral_peaks(J3_2, freq_grid)
    
    results_type2.append({
        'i': idx_i, 'j': idx_j, 'k': idx_k,
        'mode_i': mode_i, 'mode_j': mode_j, 'mode_k': mode_k,
        'freq_i': freq_i, 'freq_j': freq_j, 'freq_k': freq_k,
        'corr3': corr_val,
        'C3': C3_2,
        'J3': J3_2,
        'peaks': peaks_2,
        'heights': heights_2,
        'type': 'type2'
    })

print(f"\nTotal Type 1 results: {len(results_type1)}")
print(f"Total Type 2 results: {len(results_type2)}")

# ============================================================================
# Analyze peak frequencies
# ============================================================================
print("\n" + "="*70)
print("Analysis of Peak Frequencies in Third-Order Spectra")
print("="*70)

# Collect all peaks and their expected frequencies
all_peaks_type1 = []
all_peaks_type2 = []

print("\n--- Type 1: <q_i(t)q_j(t)q_k(0)> ---")
print("Expected peaks: ν_i ± ν_j (from q_i*q_j oscillation)")
print("\nTop 20 results by spectral integral:")

# Sort by spectral integral
results_type1_sorted = sorted(results_type1, 
                               key=lambda x: np.trapezoid(np.abs(x['J3']), freq_grid), 
                               reverse=True)

for rank, res in enumerate(results_type1_sorted[:20], 1):
    freq_i = res['freq_i']
    freq_j = res['freq_j']
    freq_k = res['freq_k']
    
    # Expected peaks for Type 1: ν_i ± ν_j
    expected_sum = freq_i + freq_j
    expected_diff = abs(freq_i - freq_j)
    
    print(f"\n{rank:2d}. Mode {res['mode_i']:2d}({freq_i:.0f}) × {res['mode_j']:2d}({freq_j:.0f}) × {res['mode_k']:2d}({freq_k:.0f})")
    print(f"    <q³> = {res['corr3']:+.4e}")
    print(f"    Expected: {expected_diff:.0f} (diff), {expected_sum:.0f} (sum) cm⁻¹")
    
    if len(res['peaks']) > 0:
        print(f"    Detected peaks:")
        for pk, ht in zip(res['peaks'][:5], res['heights'][:5]):
            # Check if close to expected
            marker = ""
            if abs(pk - expected_diff) < 50:
                marker = " ← DIFF(i-j)"
            elif abs(pk - expected_sum) < 50:
                marker = " ← SUM(i+j)"
            elif abs(pk - freq_i) < 30:
                marker = " ← ν_i"
            elif abs(pk - freq_j) < 30:
                marker = " ← ν_j"
            elif abs(pk - freq_k) < 30:
                marker = " ← ν_k"
            print(f"      {pk:7.1f} cm⁻¹ (height {ht:+.3e}){marker}")
            all_peaks_type1.append(pk)
    else:
        print(f"    No significant peaks")

print("\n\n--- Type 2: <q_k(t)q_i(0)q_j(0)> ---")
print("Expected peaks: ν_k (from q_k oscillation)")
print("\nTop 20 results by spectral integral:")

results_type2_sorted = sorted(results_type2, 
                               key=lambda x: np.trapezoid(np.abs(x['J3']), freq_grid), 
                               reverse=True)

for rank, res in enumerate(results_type2_sorted[:20], 1):
    freq_i = res['freq_i']
    freq_j = res['freq_j']
    freq_k = res['freq_k']
    
    # Expected peak for Type 2: ν_k
    expected_peak = freq_k
    
    print(f"\n{rank:2d}. Mode {res['mode_i']:2d}({freq_i:.0f}) × {res['mode_j']:2d}({freq_j:.0f}) × {res['mode_k']:2d}({freq_k:.0f})")
    print(f"    <q³> = {res['corr3']:+.4e}")
    print(f"    Expected: {expected_peak:.0f} cm⁻¹ (ν_k)")
    
    if len(res['peaks']) > 0:
        print(f"    Detected peaks:")
        for pk, ht in zip(res['peaks'][:5], res['heights'][:5]):
            marker = ""
            if abs(pk - freq_k) < 30:
                marker = " ← ν_k"
            elif abs(pk - freq_i) < 30:
                marker = " ← ν_i"
            elif abs(pk - freq_j) < 30:
                marker = " ← ν_j"
            print(f"      {pk:7.1f} cm⁻¹ (height {ht:+.3e}){marker}")
            all_peaks_type2.append(pk)
    else:
        print(f"    No significant peaks")

# ============================================================================
# Statistical summary of peak frequencies
# ============================================================================
print("\n" + "="*70)
print("Statistical Summary of Peak Frequencies")
print("="*70)

if len(all_peaks_type1) > 0:
    print(f"\nType 1 <q_i(t)q_j(t)q_k(0)>:")
    print(f"  Total peaks detected: {len(all_peaks_type1)}")
    print(f"  Range: {min(all_peaks_type1):.1f} - {max(all_peaks_type1):.1f} cm⁻¹")
    print(f"  Mean: {np.mean(all_peaks_type1):.1f} cm⁻¹")
    print(f"  Median: {np.median(all_peaks_type1):.1f} cm⁻¹")

if len(all_peaks_type2) > 0:
    print(f"\nType 2 <q_k(t)q_i(0)q_j(0)>:")
    print(f"  Total peaks detected: {len(all_peaks_type2)}")
    print(f"  Range: {min(all_peaks_type2):.1f} - {max(all_peaks_type2):.1f} cm⁻¹")
    print(f"  Mean: {np.mean(all_peaks_type2):.1f} cm⁻¹")
    print(f"  Median: {np.median(all_peaks_type2):.1f} cm⁻¹")

# ============================================================================
# Find triplets with peaks >= 500 cm⁻¹
# ============================================================================
print("\n" + "="*70)
print("Triplets with peaks >= 500 cm⁻¹")
print("="*70)

high_freq_threshold = 500.0

# Type 1 results with high frequency peaks
print("\n--- Type 1: <q_i(t)q_j(t)q_k(0)> with peaks >= 500 cm⁻¹ ---")
type1_high_freq = []
for res in results_type1:
    if len(res['peaks']) > 0:
        high_peaks = res['peaks'][res['peaks'] >= high_freq_threshold]
        if len(high_peaks) > 0:
            max_peak = np.max(high_peaks)
            # Get corresponding height
            idx_max = np.where(res['peaks'] == max_peak)[0][0]
            max_height = res['heights'][idx_max]
            type1_high_freq.append({
                'res': res,
                'high_peaks': high_peaks,
                'max_peak': max_peak,
                'max_height': max_height
            })

# Sort by max peak frequency
type1_high_freq.sort(key=lambda x: x['max_peak'], reverse=True)

print(f"\nFound {len(type1_high_freq)} triplets with peaks >= {high_freq_threshold:.0f} cm⁻¹")
print("\nAll triplets (sorted by highest peak frequency):")
for rank, item in enumerate(type1_high_freq, 1):
    res = item['res']
    freq_i, freq_j, freq_k = res['freq_i'], res['freq_j'], res['freq_k']
    mode_i, mode_j, mode_k = res['mode_i'], res['mode_j'], res['mode_k']
    
    print(f"\n{rank:3d}. Mode {mode_i:2d}({freq_i:.0f}) × {mode_j:2d}({freq_j:.0f}) × {mode_k:2d}({freq_k:.0f})")
    print(f"     <q³> = {res['corr3']:+.4e}")
    print(f"     Expected sum: ν_i+ν_j = {freq_i+freq_j:.0f} cm⁻¹")
    print(f"     High-freq peaks (>= 500 cm⁻¹):")
    for pk in item['high_peaks']:
        idx_pk = np.where(res['peaks'] == pk)[0][0]
        ht = res['heights'][idx_pk]
        # Check if close to sum frequency
        marker = ""
        if abs(pk - (freq_i + freq_j)) < 50:
            marker = " ← SUM(i+j)"
        elif abs(pk - (freq_j + freq_k)) < 50:
            marker = " ← SUM(j+k)"
        elif abs(pk - (freq_i + freq_k)) < 50:
            marker = " ← SUM(i+k)"
        print(f"       {pk:7.1f} cm⁻¹ (height {ht:+.3e}){marker}")

# Type 2 results with high frequency peaks
print("\n\n--- Type 2: <q_k(t)q_i(0)q_j(0)> with peaks >= 500 cm⁻¹ ---")
type2_high_freq = []
for res in results_type2:
    if len(res['peaks']) > 0:
        high_peaks = res['peaks'][res['peaks'] >= high_freq_threshold]
        if len(high_peaks) > 0:
            max_peak = np.max(high_peaks)
            idx_max = np.where(res['peaks'] == max_peak)[0][0]
            max_height = res['heights'][idx_max]
            type2_high_freq.append({
                'res': res,
                'high_peaks': high_peaks,
                'max_peak': max_peak,
                'max_height': max_height
            })

type2_high_freq.sort(key=lambda x: x['max_peak'], reverse=True)

print(f"\nFound {len(type2_high_freq)} triplets with peaks >= {high_freq_threshold:.0f} cm⁻¹")
print("\nAll triplets (sorted by highest peak frequency):")
for rank, item in enumerate(type2_high_freq, 1):
    res = item['res']
    freq_i, freq_j, freq_k = res['freq_i'], res['freq_j'], res['freq_k']
    mode_i, mode_j, mode_k = res['mode_i'], res['mode_j'], res['mode_k']
    
    print(f"\n{rank:3d}. Mode {mode_i:2d}({freq_i:.0f}) × {mode_j:2d}({freq_j:.0f}) × {mode_k:2d}({freq_k:.0f})")
    print(f"     <q³> = {res['corr3']:+.4e}")
    print(f"     Expected: ν_k = {freq_k:.0f} cm⁻¹")
    print(f"     High-freq peaks (>= 500 cm⁻¹):")
    for pk in item['high_peaks']:
        idx_pk = np.where(res['peaks'] == pk)[0][0]
        ht = res['heights'][idx_pk]
        marker = ""
        if abs(pk - freq_k) < 50:
            marker = " ← ν_k"
        elif abs(pk - freq_i) < 50:
            marker = " ← ν_i"
        elif abs(pk - freq_j) < 50:
            marker = " ← ν_j"
        print(f"       {pk:7.1f} cm⁻¹ (height {ht:+.3e}){marker}")

# ============================================================================
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Generating visualizations...")
print("="*70)

# Figure 1: Top 6 Type 1 spectra
fig1, axes1 = plt.subplots(2, 3, figsize=(15, 10))
axes1 = axes1.flatten()

for idx, res in enumerate(results_type1_sorted[:6]):
    ax = axes1[idx]
    ax.plot(freq_grid, res['J3'], 'b-', linewidth=1)
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    
    # Mark expected peaks
    freq_i, freq_j = res['freq_i'], res['freq_j']
    ax.axvline(freq_i + freq_j, color='r', linestyle='--', alpha=0.7, label=f'Sum={freq_i+freq_j:.0f}')
    ax.axvline(abs(freq_i - freq_j), color='g', linestyle='--', alpha=0.7, label=f'Diff={abs(freq_i-freq_j):.0f}')
    
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('FT[C3]')
    ax.set_title(f"Mode {res['mode_i']}({freq_i:.0f})×{res['mode_j']}({freq_j:.0f})×{res['mode_k']}({res['freq_k']:.0f})")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 2500)
    ax.grid(True, alpha=0.3)

plt.suptitle('Type 1: <q_i(t)q_j(t)q_k(0)> Fourier Transform', fontsize=14)
plt.tight_layout()
plt.savefig('third_order_type1_spectra.png', dpi=150, bbox_inches='tight')
print("Saved: third_order_type1_spectra.png")
plt.close()

# Figure 2: Top 6 Type 2 spectra
fig2, axes2 = plt.subplots(2, 3, figsize=(15, 10))
axes2 = axes2.flatten()

for idx, res in enumerate(results_type2_sorted[:6]):
    ax = axes2[idx]
    ax.plot(freq_grid, res['J3'], 'b-', linewidth=1)
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    
    # Mark expected peak (ν_k)
    freq_k = res['freq_k']
    ax.axvline(freq_k, color='r', linestyle='--', alpha=0.7, label=f'ν_k={freq_k:.0f}')
    
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('FT[C3]')
    ax.set_title(f"Mode {res['mode_i']}({res['freq_i']:.0f})×{res['mode_j']}({res['freq_j']:.0f})×{res['mode_k']}({freq_k:.0f})")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 2500)
    ax.grid(True, alpha=0.3)

plt.suptitle('Type 2: <q_k(t)q_i(0)q_j(0)> Fourier Transform', fontsize=14)
plt.tight_layout()
plt.savefig('third_order_type2_spectra.png', dpi=150, bbox_inches='tight')
print("Saved: third_order_type2_spectra.png")
plt.close()

# Figure 3: Histogram of peak frequencies
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))

if len(all_peaks_type1) > 0:
    axes3[0].hist(all_peaks_type1, bins=50, range=(0, 3000), color='blue', alpha=0.7, edgecolor='black')
    axes3[0].set_xlabel('Peak Frequency (cm⁻¹)')
    axes3[0].set_ylabel('Count')
    axes3[0].set_title(f'Type 1 Peak Distribution (n={len(all_peaks_type1)})')
    axes3[0].grid(True, alpha=0.3)

if len(all_peaks_type2) > 0:
    axes3[1].hist(all_peaks_type2, bins=50, range=(0, 3000), color='green', alpha=0.7, edgecolor='black')
    axes3[1].set_xlabel('Peak Frequency (cm⁻¹)')
    axes3[1].set_ylabel('Count')
    axes3[1].set_title(f'Type 2 Peak Distribution (n={len(all_peaks_type2)})')
    axes3[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('third_order_peak_histogram.png', dpi=150, bbox_inches='tight')
print("Saved: third_order_peak_histogram.png")
plt.close()

# Figure 4: Correlation functions (time domain)
fig4, axes4 = plt.subplots(2, 3, figsize=(15, 10))

# Top 3 Type 1
for idx, res in enumerate(results_type1_sorted[:3]):
    ax = axes4[0, idx]
    t_axis = np.arange(max_lag) * dt / 1000  # ps
    ax.plot(t_axis, res['C3'], 'b-', linewidth=1)
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('C3(t)')
    ax.set_title(f"Type1: {res['mode_i']}×{res['mode_j']}×{res['mode_k']}")
    ax.grid(True, alpha=0.3)

# Top 3 Type 2
for idx, res in enumerate(results_type2_sorted[:3]):
    ax = axes4[1, idx]
    t_axis = np.arange(max_lag) * dt / 1000  # ps
    ax.plot(t_axis, res['C3'], 'g-', linewidth=1)
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('C3(t)')
    ax.set_title(f"Type2: {res['mode_i']}×{res['mode_j']}×{res['mode_k']}")
    ax.grid(True, alpha=0.3)

plt.suptitle('Third-Order Correlation Functions (Time Domain)', fontsize=14)
plt.tight_layout()
plt.savefig('third_order_correlation_time.png', dpi=150, bbox_inches='tight')
print("Saved: third_order_correlation_time.png")
plt.close()

# ============================================================================
# Compare 3rd and 4th order correlations for triplets with peaks >= 500 cm⁻¹
# ============================================================================
print("\n" + "="*70)
print("Comparing 3rd and 4th Order Correlations")
print("="*70)

# Select triplets with high-frequency peaks from Type 1
comparison_data = []
n_compare = min(50, len(type1_high_freq))  # Compare top 50 or all available

print(f"\nComputing 4th order correlations for {n_compare} triplets with peaks >= 500 cm⁻¹...")

for idx, item in enumerate(type1_high_freq[:n_compare]):
    if (idx + 1) % 10 == 0:
        print(f"  Processing {idx + 1}/{n_compare}...")
    
    res = item['res']
    idx_i, idx_j, idx_k = res['i'], res['j'], res['k']
    
    q_i = mode_coords[:, idx_i]
    q_j = mode_coords[:, idx_j]
    q_k = mode_coords[:, idx_k]
    
    # 3rd order: <q_i(t)q_j(t)q_k(0)>
    C3 = res['C3']
    J3 = res['J3']
    
    # 4th order variations:
    # C4a: <q_i(t)q_j(t)q_k(0)q_k(0)> - same k at t=0
    C4a = compute_fourth_order_corr(q_i, q_j, q_k, q_k)
    J4a = compute_fourier_transform(C4a, freq_grid)
    
    # C4b: <q_i(t)q_j(t)q_i(0)q_j(0)> - same pair at t and 0
    C4b = compute_fourth_order_corr(q_i, q_j, q_i, q_j)
    J4b = compute_fourier_transform(C4b, freq_grid)
    
    # C4c: <q_i(t)q_k(t)q_j(0)q_k(0)> - mixed
    C4c = compute_fourth_order_corr(q_i, q_k, q_j, q_k)
    J4c = compute_fourier_transform(C4c, freq_grid)
    
    # Compute spectral integrals (absolute values for comparison)
    int_J3 = np.trapezoid(np.abs(J3), freq_grid)
    int_J4a = np.trapezoid(np.abs(J4a), freq_grid)
    int_J4b = np.trapezoid(np.abs(J4b), freq_grid)
    int_J4c = np.trapezoid(np.abs(J4c), freq_grid)
    
    # C(0) values (equal-time correlations)
    C3_0 = C3[0]
    C4a_0 = C4a[0]
    C4b_0 = C4b[0]
    C4c_0 = C4c[0]
    
    comparison_data.append({
        'res': res,
        'max_peak': item['max_peak'],
        'C3': C3, 'J3': J3, 'C3_0': C3_0, 'int_J3': int_J3,
        'C4a': C4a, 'J4a': J4a, 'C4a_0': C4a_0, 'int_J4a': int_J4a,
        'C4b': C4b, 'J4b': J4b, 'C4b_0': C4b_0, 'int_J4b': int_J4b,
        'C4c': C4c, 'J4c': J4c, 'C4c_0': C4c_0, 'int_J4c': int_J4c,
    })

# Print comparison summary
print("\n" + "-"*70)
print("Comparison Summary: 3rd vs 4th Order Correlations")
print("-"*70)
print(f"{'Rank':>4} {'Modes':>20} {'Peak':>8} {'|C3(0)|':>12} {'|C4a(0)|':>12} {'|C4b(0)|':>12} {'Ratio':>10}")
print("-"*70)

for idx, data in enumerate(comparison_data[:30], 1):
    res = data['res']
    modes_str = f"{res['mode_i']}×{res['mode_j']}×{res['mode_k']}"
    
    c3_0 = np.abs(data['C3_0'])
    c4a_0 = np.abs(data['C4a_0'])
    c4b_0 = np.abs(data['C4b_0'])
    
    ratio = c4a_0 / c3_0 if c3_0 > 0 else 0
    
    print(f"{idx:4d} {modes_str:>20} {data['max_peak']:8.1f} {c3_0:12.4e} {c4a_0:12.4e} {c4b_0:12.4e} {ratio:10.2f}")

# ============================================================================
# Visualization: 3rd vs 4th order comparison
# ============================================================================
print("\n" + "="*70)
print("Generating comparison visualizations...")
print("="*70)

# Figure 5: Bar chart comparison of C(0) values
fig5, axes5 = plt.subplots(2, 2, figsize=(14, 12))

# Prepare data for bar charts
n_bars = min(20, len(comparison_data))
labels = [f"{d['res']['mode_i']}×{d['res']['mode_j']}×{d['res']['mode_k']}" for d in comparison_data[:n_bars]]
c3_0_vals = [np.abs(d['C3_0']) for d in comparison_data[:n_bars]]
c4a_0_vals = [np.abs(d['C4a_0']) for d in comparison_data[:n_bars]]
c4b_0_vals = [np.abs(d['C4b_0']) for d in comparison_data[:n_bars]]
c4c_0_vals = [np.abs(d['C4c_0']) for d in comparison_data[:n_bars]]

x = np.arange(n_bars)
width = 0.2

# Plot 1: C(0) comparison
ax1 = axes5[0, 0]
ax1.bar(x - 1.5*width, c3_0_vals, width, label='C3(0)', color='blue', alpha=0.7)
ax1.bar(x - 0.5*width, c4a_0_vals, width, label='C4a(0): ij-kk', color='red', alpha=0.7)
ax1.bar(x + 0.5*width, c4b_0_vals, width, label='C4b(0): ij-ij', color='green', alpha=0.7)
ax1.bar(x + 1.5*width, c4c_0_vals, width, label='C4c(0): ik-jk', color='orange', alpha=0.7)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax1.set_ylabel('|C(0)|')
ax1.set_title('Equal-time Correlation Comparison: 3rd vs 4th Order')
ax1.legend(fontsize=8)
ax1.set_yscale('log')
ax1.grid(True, alpha=0.3, axis='y')

# Plot 2: Spectral integral comparison
int_j3_vals = [d['int_J3'] for d in comparison_data[:n_bars]]
int_j4a_vals = [d['int_J4a'] for d in comparison_data[:n_bars]]
int_j4b_vals = [d['int_J4b'] for d in comparison_data[:n_bars]]
int_j4c_vals = [d['int_J4c'] for d in comparison_data[:n_bars]]

ax2 = axes5[0, 1]
ax2.bar(x - 1.5*width, int_j3_vals, width, label='∫|J3|', color='blue', alpha=0.7)
ax2.bar(x - 0.5*width, int_j4a_vals, width, label='∫|J4a|', color='red', alpha=0.7)
ax2.bar(x + 0.5*width, int_j4b_vals, width, label='∫|J4b|', color='green', alpha=0.7)
ax2.bar(x + 1.5*width, int_j4c_vals, width, label='∫|J4c|', color='orange', alpha=0.7)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax2.set_ylabel('Spectral Integral')
ax2.set_title('Spectral Integral Comparison: 3rd vs 4th Order')
ax2.legend(fontsize=8)
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3, axis='y')

# Plot 3: Ratio C4a(0)/C3(0) vs peak frequency
ax3 = axes5[1, 0]
peak_freqs = [d['max_peak'] for d in comparison_data]
ratios_4a = [np.abs(d['C4a_0'])/np.abs(d['C3_0']) if np.abs(d['C3_0']) > 0 else 0 for d in comparison_data]
ratios_4b = [np.abs(d['C4b_0'])/np.abs(d['C3_0']) if np.abs(d['C3_0']) > 0 else 0 for d in comparison_data]

ax3.scatter(peak_freqs, ratios_4a, alpha=0.6, label='C4a/C3', color='red', s=30)
ax3.scatter(peak_freqs, ratios_4b, alpha=0.6, label='C4b/C3', color='green', s=30)
ax3.axhline(1, color='k', linestyle='--', linewidth=1, label='Ratio=1')
ax3.set_xlabel('Peak Frequency (cm⁻¹)')
ax3.set_ylabel('C4(0)/C3(0) Ratio')
ax3.set_title('4th/3rd Order Ratio vs Peak Frequency')
ax3.legend()
ax3.set_yscale('log')
ax3.grid(True, alpha=0.3)

# Plot 4: Scatter plot C3(0) vs C4a(0)
ax4 = axes5[1, 1]
ax4.scatter(c3_0_vals, c4a_0_vals, alpha=0.6, label='C4a vs C3', color='red', s=50)
ax4.scatter(c3_0_vals, c4b_0_vals, alpha=0.6, label='C4b vs C3', color='green', s=50)

# Add diagonal line
min_val = min(min(c3_0_vals), min(c4a_0_vals), min(c4b_0_vals))
max_val = max(max(c3_0_vals), max(c4a_0_vals), max(c4b_0_vals))
ax4.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1, label='y=x')
ax4.set_xlabel('|C3(0)|')
ax4.set_ylabel('|C4(0)|')
ax4.set_title('3rd vs 4th Order Correlation Magnitudes')
ax4.legend()
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('third_fourth_order_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: third_fourth_order_comparison.png")
plt.close()

# Figure 6: Spectral comparison for top 6 triplets
fig6, axes6 = plt.subplots(3, 2, figsize=(14, 15))
axes6 = axes6.flatten()

for idx, data in enumerate(comparison_data[:6]):
    ax = axes6[idx]
    res = data['res']
    
    # Normalize for visualization
    max_val = max(np.max(np.abs(data['J3'])), np.max(np.abs(data['J4a'])), 
                  np.max(np.abs(data['J4b'])), 1e-10)
    
    ax.plot(freq_grid, data['J3']/max_val, 'b-', linewidth=1.5, label='J3', alpha=0.8)
    ax.plot(freq_grid, data['J4a']/max_val, 'r--', linewidth=1.5, label='J4a (ij-kk)', alpha=0.8)
    ax.plot(freq_grid, data['J4b']/max_val, 'g:', linewidth=1.5, label='J4b (ij-ij)', alpha=0.8)
    
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    ax.axvline(500, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Normalized J(ω)')
    ax.set_title(f"Mode {res['mode_i']}×{res['mode_j']}×{res['mode_k']} (peak: {data['max_peak']:.0f} cm⁻¹)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1500)
    ax.grid(True, alpha=0.3)

plt.suptitle('Spectral Comparison: 3rd vs 4th Order Correlations', fontsize=14)
plt.tight_layout()
plt.savefig('spectral_3rd_4th_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_3rd_4th_comparison.png")
plt.close()

# Figure 7: Time-domain comparison for top 6 triplets
fig7, axes7 = plt.subplots(3, 2, figsize=(14, 15))
axes7 = axes7.flatten()
t_axis = np.arange(max_lag) * dt / 1000  # ps

for idx, data in enumerate(comparison_data[:6]):
    ax = axes7[idx]
    res = data['res']
    
    # Normalize for visualization
    max_val = max(np.max(np.abs(data['C3'])), np.max(np.abs(data['C4a'])), 
                  np.max(np.abs(data['C4b'])), 1e-10)
    
    ax.plot(t_axis, data['C3']/max_val, 'b-', linewidth=1.5, label='C3', alpha=0.8)
    ax.plot(t_axis, data['C4a']/max_val, 'r--', linewidth=1.5, label='C4a (ij-kk)', alpha=0.8)
    ax.plot(t_axis, data['C4b']/max_val, 'g:', linewidth=1.5, label='C4b (ij-ij)', alpha=0.8)
    
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Normalized C(t)')
    ax.set_title(f"Mode {res['mode_i']}×{res['mode_j']}×{res['mode_k']}")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 5)
    ax.grid(True, alpha=0.3)

plt.suptitle('Time-domain Comparison: 3rd vs 4th Order Correlations', fontsize=14)
plt.tight_layout()
plt.savefig('time_domain_3rd_4th_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: time_domain_3rd_4th_comparison.png")
plt.close()

print("\n" + "="*70)
print("Analysis complete!")
print("Log saved to: third_order_correlation.log")
print("="*70)
