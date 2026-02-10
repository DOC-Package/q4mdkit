#!/usr/bin/env python3
"""
Two-stage spectral domain fitting:
1. Linear terms with NNLS (c_k ≥ 0)
2. Cross terms on residual (d_ij unrestricted)

This combines physical constraints on linear terms with flexibility for cross terms.
"""
import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks
from scipy.optimize import nnls
from sklearn.linear_model import Ridge
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("Two-Stage Spectral Fitting: Linear (NNLS) + Cross Terms")
print("="*70)

# Load data
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff = energy_diff_data[:, 4]

# Load mode information
with open('g_k_coefficients_below_2000cm.dat') as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith('#') and 'R²' not in l]
    data = [l.split() for l in lines if len(l.split()) >= 3]

mode_indices = np.array([int(d[0]) for d in data])
mode_freqs = np.array([float(d[1]) for d in data])

# Filter modes for Stage 1: 450 ≤ freq ≤ 1800 cm⁻¹
stage1_mask = (mode_freqs >= 450) & (mode_freqs <= 1800)
print(f"\nTotal modes: {len(mode_indices)}")
print(f"Stage 1 modes (450-1800 cm⁻¹): {np.sum(stage1_mask)}")

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

def compute_autocorr(x):
    """Compute autocorrelation with Hanning window"""
    x_mean = np.mean(x)
    x_centered = x - x_mean
    N = len(x)
    
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(x_centered[:N-lag] * x_centered[lag:])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return acf * window

def compute_crosscorr(x, y):
    """Compute cross-correlation with Hanning window"""
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_centered = x - x_mean
    y_centered = y - y_mean
    N = len(x)
    
    ccf = np.zeros(max_lag)
    for lag in range(max_lag):
        ccf[lag] = np.mean(x_centered[:N-lag] * y_centered[lag:])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return ccf * window

def compute_spectral_density(x, freq):
    """Compute spectral density using Simpson integration"""
    C = compute_autocorr(x)
    t = np.arange(max_lag) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, t)
    
    # Multiply by prefactor
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

def compute_cross_spectral_density(x, y, freq):
    """Compute cross-spectral density using Simpson integration"""
    C = compute_crosscorr(x, y)
    t = np.arange(max_lag) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, t)
    
    # Multiply by prefactor
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

def find_spectral_peaks(J, freq, threshold=0.1, min_separation=20):
    """
    Find peaks in spectral density.
    
    Parameters:
    - J: spectral density array
    - freq: frequency grid
    - threshold: minimum peak height relative to max(J)
    - min_separation: minimum separation between peaks (cm⁻¹)
    
    Returns:
    - peak_freqs: frequencies of detected peaks
    - peak_heights: heights of peaks
    """
    # Find peaks with minimum height and prominence
    max_J = np.max(J)
    if max_J <= 0:
        return np.array([]), np.array([])
    
    peaks_idx, properties = find_peaks(
        J, 
        height=threshold * max_J,
        prominence=threshold * max_J * 0.5,
        distance=int(min_separation / (freq[1] - freq[0]))
    )
    
    if len(peaks_idx) == 0:
        return np.array([]), np.array([])
    
    peak_freqs = freq[peaks_idx]
    peak_heights = J[peaks_idx]
    
    # Sort by height (descending)
    sorted_idx = np.argsort(peak_heights)[::-1]
    
    return peak_freqs[sorted_idx], peak_heights[sorted_idx]

# Frequency grid (420-1800 cm⁻¹)
freq_grid = np.linspace(420, 1800, 1381)
print(f"\nFrequency grid: {len(freq_grid)} points (420-1800 cm⁻¹)")

print("\nComputing J_total...")
J_total = compute_spectral_density(energy_diff, freq_grid)
J_total_integral = np.trapezoid(J_total, freq_grid)
print(f"  ∫J_total = {J_total_integral:.6e} cm⁻¹²")

print("\nComputing J_k for all modes...")
J_k_matrix = np.zeros((len(freq_grid), N_modes))
J_k_integrals = np.zeros(N_modes)
for k in range(N_modes):
    if (k+1) % 10 == 0 or k == 0:
        print(f"  Mode {k+1}/{N_modes}")
    J_k_matrix[:, k] = compute_spectral_density(mode_coords[:, k], freq_grid)
    J_k_integrals[k] = np.trapezoid(J_k_matrix[:, k], freq_grid)

# ============================================================================
# Stage 1: Linear terms with NNLS (only 450-1800 cm⁻¹ modes)
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, c_k ≥ 0)")
print(f"Using {np.sum(stage1_mask)} modes with 450 ≤ freq ≤ 1800 cm⁻¹")
print("="*70)

# Use only filtered modes for Stage 1
J_k_stage1 = J_k_matrix[:, stage1_mask]
J_k_integrals_stage1 = J_k_integrals[stage1_mask]
n_modes_stage1 = np.sum(stage1_mask)

# Normalize
J_k_normalized = J_k_stage1 / J_k_integrals_stage1[np.newaxis, :]
J_total_normalized = J_total / J_total_integral

# NNLS fit
print("\nFitting linear terms with NNLS...")
c_k_nnls_stage1, residual_norm = nnls(J_k_normalized, J_total_normalized)
c_k_physical_stage1 = c_k_nnls_stage1 * J_total_integral / J_k_integrals_stage1

# Expand to full mode array
c_k_physical = np.zeros(N_modes)
c_k_physical[stage1_mask] = c_k_physical_stage1

# Reconstruct
J_linear = J_k_matrix @ c_k_physical
J_linear_integral = np.trapezoid(J_linear, freq_grid)

# Residual
J_residual_1 = J_total - J_linear
J_residual_1_integral = np.trapezoid(J_residual_1, freq_grid)

# R²
ss_res_1 = np.sum(J_residual_1**2)
ss_tot = np.sum((J_total - np.mean(J_total))**2)
R2_linear = 1 - ss_res_1 / ss_tot

n_active = np.sum(c_k_physical > 1e-10)

print(f"\nLinear model results:")
print(f"  R² = {R2_linear:.4f}")
print(f"  ∫J_linear = {J_linear_integral:.6e} cm⁻¹² ({100*J_linear_integral/J_total_integral:.1f}%)")
print(f"  ∫J_residual = {J_residual_1_integral:.6e} cm⁻¹² ({100*J_residual_1_integral/J_total_integral:.1f}%)")
print(f"  Active modes: {n_active}/{N_modes}")

# ============================================================================
# Stage 2: Cross terms on residual (specific frequency range)
# ============================================================================
print("\n" + "="*70)
print("Stage 2: Cross Terms on Residual (400-800 cm⁻¹)")
print("="*70)

# Select modes for cross terms
freq_min, freq_max = 0, 800
mask_cross = (mode_freqs >= freq_min) & (mode_freqs <= freq_max)
idx_cross = np.where(mask_cross)[0]
n_cross_modes = len(idx_cross)

print(f"\nModes in {freq_min}-{freq_max} cm⁻¹: {n_cross_modes}")
print(f"Number of quadratic terms: {n_cross_modes * (n_cross_modes + 1) // 2} (diagonal: {n_cross_modes}, off-diagonal: {n_cross_modes * (n_cross_modes - 1) // 2})")

if n_cross_modes >= 2:
    # Compute cross term spectral densities
    print("\nComputing cross term spectral densities...")
    cross_pairs = []
    J_cross_list = []
    
    for i_idx, i in enumerate(idx_cross):
        # Diagonal term: q_i^2
        q_diag = mode_coords[:, i] ** 2
        J_diag = compute_spectral_density(q_diag, freq_grid)
        J_cross_list.append(J_diag)
        cross_pairs.append((i, i))
        
        # Off-diagonal terms: q_i * q_j (i < j)
        for j in idx_cross[i_idx+1:]:
            q_cross = mode_coords[:, i] * mode_coords[:, j]
            J_cross = compute_spectral_density(q_cross, freq_grid)
            J_cross_list.append(J_cross)
            cross_pairs.append((i, j))
        
        if (i_idx + 1) % 5 == 0:
            print(f"  Progress: {i_idx+1}/{n_cross_modes} modes")
    
    n_cross = len(cross_pairs)
    J_cross_matrix = np.column_stack(J_cross_list)
    
    print(f"\nTotal cross terms computed: {n_cross}")
    
    # ========================================================================
    # Add high-frequency mode combinations for 540 cm⁻¹ peak
    # ========================================================================
    print("\n" + "="*70)
    print("Adding high-frequency mode pairs for 540 cm⁻¹ peak")
    print("="*70)
    
    target_freq = 540.0  # cm⁻¹
    tolerance = 100.0    # cm⁻¹
    
    # Find mode pairs where |ν_i - ν_j| ≈ 540 cm⁻¹
    additional_pairs = []
    for i in range(N_modes):
        for j in range(i+1, N_modes):
            freq_i = mode_freqs[i]
            freq_j = mode_freqs[j]
            diff_freq = abs(freq_i - freq_j)
            
            # Check if difference is near 540 cm⁻¹
            if abs(diff_freq - target_freq) < tolerance:
                # Skip if already included in idx_cross range
                if i in idx_cross and j in idx_cross:
                    continue
                
                additional_pairs.append((i, j, diff_freq))
    
    print(f"\nFound {len(additional_pairs)} additional mode pairs with |Δν| ≈ 540 cm⁻¹")
    
    if len(additional_pairs) > 0:
        print("\nTop 10 additional pairs by frequency difference:")
        # Sort by how close to 540 cm⁻¹
        additional_pairs_sorted = sorted(additional_pairs, key=lambda x: abs(x[2] - target_freq))
        for rank, (i, j, diff) in enumerate(additional_pairs_sorted[:10], 1):
            print(f"  {rank:2d}. Mode {mode_indices[i]:2d} ({mode_freqs[i]:6.1f}) × Mode {mode_indices[j]:2d} ({mode_freqs[j]:6.1f}): Δν = {diff:.1f} cm⁻¹")
        
        print(f"\nComputing spectral densities for additional pairs...")
        for idx, (i, j, diff) in enumerate(additional_pairs):
            q_cross = mode_coords[:, i] * mode_coords[:, j]
            J_cross = compute_spectral_density(q_cross, freq_grid)
            J_cross_list.append(J_cross)
            cross_pairs.append((i, j))
            
            if (idx + 1) % 10 == 0:
                print(f"  Progress: {idx+1}/{len(additional_pairs)} pairs")
        
        # Update matrix
        n_cross_updated = len(cross_pairs)
        J_cross_matrix = np.column_stack(J_cross_list)
        print(f"\nTotal cross terms (including 540 cm⁻¹ pairs): {n_cross_updated}")
        n_cross = n_cross_updated
    
    # ========================================================================
    # Analyze peak positions in cross terms
    # ========================================================================
    print("\n" + "="*70)
    print("Peak Detection in Cross Term Spectral Densities")
    print("="*70)
    
    # Analyze peaks for top 10 cross terms by integral
    cross_integrals = np.array([np.trapezoid(J, freq_grid) for J in J_cross_list])
    top_cross_by_integral = np.argsort(np.abs(cross_integrals))[-10:][::-1]
    
    peak_data = []
    print("\nTop 10 cross terms by spectral integral:")
    for rank, idx in enumerate(top_cross_by_integral, 1):
        i, j = cross_pairs[idx]
        freq_i = mode_freqs[i]
        freq_j = mode_freqs[j]
        
        # Expected peak positions
        sum_freq = freq_i + freq_j
        diff_freq = abs(freq_i - freq_j)
        
        # Detect actual peaks
        J_cross_single = J_cross_list[idx]
        peak_freqs, peak_heights = find_spectral_peaks(J_cross_single, freq_grid, threshold=0.15)
        
        # Find peaks near expected positions
        sum_peak_idx = None
        diff_peak_idx = None
        
        if len(peak_freqs) > 0:
            # Look for sum peak (should be strongest)
            sum_diffs = np.abs(peak_freqs - sum_freq)
            if np.min(sum_diffs) < 50:  # Within 50 cm⁻¹
                sum_peak_idx = np.argmin(sum_diffs)
            
            # Look for difference peak
            diff_diffs = np.abs(peak_freqs - diff_freq)
            if np.min(diff_diffs) < 50 and diff_freq > 50:  # Within 50 cm⁻¹ and not too low
                diff_peak_idx = np.argmin(diff_diffs)
        
        print(f"\n{rank:2d}. Mode {mode_indices[i]:2d} ({freq_i:6.1f}) × Mode {mode_indices[j]:2d} ({freq_j:6.1f})")
        print(f"    Expected: Sum={sum_freq:7.1f}, Diff={diff_freq:6.1f} cm⁻¹")
        print(f"    Integral: {cross_integrals[idx]:+.6e} cm⁻¹²")
        
        if len(peak_freqs) > 0:
            print(f"    Detected peaks ({len(peak_freqs)} total):")
            for k, (pf, ph) in enumerate(zip(peak_freqs[:5], peak_heights[:5]), 1):
                marker = ""
                if sum_peak_idx is not None and k-1 == sum_peak_idx:
                    marker = " ← SUM"
                elif diff_peak_idx is not None and k-1 == diff_peak_idx:
                    marker = " ← DIFF"
                print(f"      {k}. {pf:7.1f} cm⁻¹ (height {ph:.3e}){marker}")
        else:
            print(f"    No significant peaks detected")
        
        peak_data.append({
            'idx': idx,
            'i': i, 'j': j,
            'freq_i': freq_i, 'freq_j': freq_j,
            'sum_expected': sum_freq,
            'diff_expected': diff_freq,
            'peaks': peak_freqs,
            'heights': peak_heights,
            'integral': cross_integrals[idx]
        })
    
    # Fit cross terms to residual (NNLS: d_ij ≥ 0)
    print("\n" + "="*70)
    print("Fitting Cross Terms to Residual (NNLS: d_ij ≥ 0)")
    print("="*70)
    
    # NNLS fitting (non-negative constraint)
    from scipy.optimize import nnls
    d_ij, residual_norm = nnls(J_cross_matrix, J_residual_1)
    
    J_cross_model = J_cross_matrix @ d_ij
    J_cross_integral = np.trapezoid(J_cross_model, freq_grid)
    
    # Total model
    J_total_model = J_linear + J_cross_model
    J_residual_2 = J_total - J_total_model
    
    ss_res_2 = np.sum(J_residual_2**2)
    R2_total = 1 - ss_res_2 / ss_tot
    
    J_total_model_integral = np.trapezoid(J_total_model, freq_grid)
    J_residual_2_integral = np.trapezoid(J_residual_2, freq_grid)
    
    n_active = np.sum(d_ij > 1e-10)
    
    print(f"  R² (total) = {R2_total:.4f} (linear: {R2_linear:.4f}, improvement: +{R2_total-R2_linear:.4f})")
    print(f"  Total capture = {100 * J_total_model_integral / J_total_integral:.1f}%")
    print(f"    Linear contribution: {100 * J_linear_integral / J_total_integral:.1f}%")
    print(f"    Cross contribution: +{100 * J_cross_integral / J_total_integral:.1f}%")
    print(f"  Active cross terms: {n_active}/{len(d_ij)} (all d_ij ≥ 0)")
    
    # Store single result
    results_cross = [{
        'alpha': 0.0,
        'R2': R2_total,
        'd_ij': d_ij,
        'J_cross': J_cross_model,
        'J_total': J_total_model,
        'capture': 100 * J_total_model_integral / J_total_integral,
        'cross_contrib': 100 * J_cross_integral / J_total_integral,
        'n_pos': n_active,
        'n_neg': 0
    }]
    
    # ========================================================================
    # Individual stage visualizations
    # ========================================================================
    print("\nGenerating individual stage visualizations...")
    
    # Stage 1 visualization
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # Top panel: J_total vs J_linear
    ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid, J_linear, 'b-', lw=1.5, label=f'Stage 1: Linear NNLS (R²={R2_linear:.4f})', alpha=0.7)
    ax1.fill_between(freq_grid, 0, J_linear, alpha=0.2, color='blue')
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Stage 1: Linear Terms Only ({100*J_linear_integral/J_total_integral:.1f}% capture)', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Bottom panel: Residual
    ax2.plot(freq_grid, J_residual_1, 'r-', lw=1.5, label=f'Residual ({100*J_residual_1_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax2.set_ylabel('Residual (cm⁻¹²)', fontsize=12)
    ax2.set_title('Stage 1 Residual', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('spectral_stage1_linear.png', dpi=150, bbox_inches='tight')
    print(f"Saved: spectral_stage1_linear.png")
    plt.close()
    
    # Stage 2 visualization
    fig2, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    # Top panel: J_total vs J_total_model
    ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid, J_linear, 'b--', lw=1, label=f'Stage 1: Linear only', alpha=0.5)
    ax1.plot(freq_grid, J_total_model, 'r-', lw=1.5, label=f'Stage 2: Linear+Quadratic (R²={R2_total:.4f})', alpha=0.8)
    ax1.fill_between(freq_grid, 0, J_total_model, alpha=0.2, color='red')
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Stage 2: Linear + Quadratic Terms ({100*J_total_model_integral/J_total_integral:.1f}% capture)', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Middle panel: Quadratic contribution
    ax2.plot(freq_grid, J_cross_model, 'purple', lw=1.5, label=f'Quadratic contribution ({100*J_cross_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax2.set_title('Stage 2: Quadratic Term Contribution', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Bottom panel: Final residual
    ax3.plot(freq_grid, J_residual_1, 'b-', lw=1, label=f'After Stage 1 ({100*J_residual_1_integral/J_total_integral:.1f}%)', alpha=0.5)
    ax3.plot(freq_grid, J_residual_2, 'r-', lw=1.5, label=f'After Stage 2 ({100*J_residual_2_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax3.axhline(0, color='k', ls='--', lw=0.5)
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax3.set_ylabel('Residual (cm⁻¹²)', fontsize=12)
    ax3.set_title('Residual Comparison', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('spectral_stage2_combined.png', dpi=150, bbox_inches='tight')
    print(f"Saved: spectral_stage2_combined.png")
    plt.close()
    
    # Select best model (only one with NNLS)
    best = results_cross[0]
    
    print(f"\n" + "="*70)
    print("Best Combined Model")
    print("="*70)
    print(f"\nRegularization: NNLS (d_ij ≥ 0)")
    print(f"  Total R² = {best['R2']:.4f} (linear: {R2_linear:.4f}, improvement: +{best['R2']-R2_linear:.4f})")
    print(f"  Total capture = {best['capture']:.1f}%")
    print(f"    Linear contribution: {100*J_linear_integral/J_total_integral:.1f}%")
    print(f"    Cross contribution: {best['cross_contrib']:+.1f}%")
    print(f"  Active cross terms: {best['n_pos']}/{len(d_ij)} (all d_ij ≥ 0)")
    
    # Top cross terms
    top_cross_idx = np.argsort(np.abs(best['d_ij']))[-10:][::-1]
    print(f"\nTop 10 cross terms by |d_ij|:")
    for rank, idx in enumerate(top_cross_idx, 1):
        i, j = cross_pairs[idx]
        sign = '+' if best['d_ij'][idx] > 0 else '-'
        print(f"  {rank:2d}. Mode {mode_indices[i]:2d} ({mode_freqs[i]:6.1f}) × "
              f"Mode {mode_indices[j]:2d} ({mode_freqs[j]:6.1f}): "
              f"{sign} {abs(best['d_ij'][idx]):.6e}")
    
    # ========================================================================
    # Visualization
    # ========================================================================
    print("\nGenerating visualization...")
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Panel 1: Spectral density decomposition
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid, J_linear, 'b--', lw=1.5, label=f'Stage 1: Linear ({100*J_linear_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax1.plot(freq_grid, best['J_total'], 'r-', lw=1.5, label=f'Stage 2: +Product+Diag ({best["capture"]:.1f}%)', alpha=0.8)
    ax1.fill_between(freq_grid, 0, best['J_total'], alpha=0.2, color='red')
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Two-Stage Spectral Density Fitting (R²={best["R2"]:.4f})', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Residual comparison
    ax2 = fig.add_subplot(gs[1, 0])
    J_residual_2 = J_total - best['J_total']
    ax2.plot(freq_grid, J_residual_1, 'b-', lw=1, alpha=0.7, label='After Stage 1')
    ax2.plot(freq_grid, J_residual_2, 'r-', lw=1, alpha=0.7, label='After Stage 2')
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax2.set_ylabel('Residual (cm⁻¹²)', fontsize=11)
    ax2.set_title('Residual Reduction', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Cross term contribution
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(freq_grid, best['J_cross'], 'r-', lw=1.5, label=f'Quadratic ({best["cross_contrib"]:+.1f}%)')
    ax3.axhline(0, color='k', ls='--', lw=0.5)
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax3.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax3.set_title('Quadratic Term Contribution', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Cross coefficient distribution
    ax4 = fig.add_subplot(gs[1, 2])
    d_sorted = np.sort(np.abs(best['d_ij']))[::-1]
    ax4.semilogy(d_sorted, 'o-', markersize=3, lw=1, color='purple')
    ax4.set_xlabel('Cross term index (sorted)', fontsize=11)
    ax4.set_ylabel('|d_ij|', fontsize=11)
    ax4.set_title(f'Cross Coefficient Distribution', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3, which='both')
    
    # Panel 5: Top linear modes
    ax5 = fig.add_subplot(gs[2, 0])
    top_linear = np.argsort(c_k_physical)[-10:][::-1]
    ax5.barh(range(10), c_k_physical[top_linear], color='blue', alpha=0.7)
    ax5.set_yticks(range(10))
    ax5.set_yticklabels([f"Mode {mode_indices[i]} ({mode_freqs[i]:.0f} cm⁻¹)" for i in top_linear], fontsize=9)
    ax5.set_xlabel('Linear coeff c_k (≥0)', fontsize=11)
    ax5.set_title('Top 10 Linear Modes', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='x')
    
    # Panel 6: Top cross terms (by absolute value)
    ax6 = fig.add_subplot(gs[2, 1])
    colors = ['green' if d > 0 else 'red' for d in best['d_ij'][top_cross_idx]]
    ax6.barh(range(10), best['d_ij'][top_cross_idx], color=colors, alpha=0.7)
    ax6.set_yticks(range(10))
    labels = []
    for idx in top_cross_idx:
        i, j = cross_pairs[idx]
        labels.append(f"{mode_indices[i]}×{mode_indices[j]} ({mode_freqs[i]:.0f},{mode_freqs[j]:.0f})")
    ax6.set_yticklabels(labels, fontsize=8)
    ax6.set_xlabel('Cross coeff d_ij', fontsize=11)
    ax6.set_title('Top 10 Cross Terms', fontsize=11, fontweight='bold')
    ax6.axvline(0, color='k', ls='-', lw=0.5)
    ax6.grid(True, alpha=0.3, axis='x')
    
    # Panel 7: R² vs alpha
    ax7 = fig.add_subplot(gs[2, 2])
    alphas_plot = [r['alpha'] for r in results_cross]
    R2_plot = [r['R2'] for r in results_cross]
    ax7.plot(alphas_plot, R2_plot, 'o-', color='purple', lw=2, markersize=6)
    ax7.axhline(R2_linear, color='blue', ls='--', lw=1, label='Linear only')
    ax7.set_xscale('log')
    ax7.set_xlabel('Cross regularization α', fontsize=11)
    ax7.set_ylabel('Total R²', fontsize=11)
    ax7.set_title('Cross Term Regularization', fontsize=11, fontweight='bold')
    ax7.legend(fontsize=9)
    ax7.grid(True, alpha=0.3)
    
    plt.savefig('spectral_fit_with_cross2.png', dpi=150, bbox_inches='tight')
    print(f"Saved: spectral_fit_with_cross2.png")
    
    # ========================================================================
    # Save results
    # ========================================================================
    print("\nSaving coefficients...")
    
    # Linear coefficients
    with open('spectral_coefficients_linear_nnls.dat', 'w') as f:
        f.write(f"# Linear coefficients from NNLS (Stage 1)\n")
        f.write(f"# All c_k ≥ 0 (physical constraint)\n")
        f.write(f"# R² (linear only) = {R2_linear:.6f}\n")
        f.write(f"# Spectral capture = {100*J_linear_integral/J_total_integral:.2f}%\n")
        f.write(f"# Active modes: {n_active}/{N_modes}\n")
        f.write(f"# Mode  Frequency(cm⁻¹)  c_k  Active\n")
        for i in range(N_modes):
            active = 'Yes' if c_k_physical[i] > 1e-10 else 'No'
            f.write(f"{mode_indices[i]:4d}  {mode_freqs[i]:10.2f}  {c_k_physical[i]:12.6e}  {active}\n")
    
    print(f"Saved: spectral_coefficients_linear_nnls.dat")
    
    # Cross coefficients
    with open('spectral_coefficients_cross.dat', 'w') as f:
        f.write(f"# Cross term coefficients (Stage 2)\n")
        f.write(f"# Fitted to residual after linear NNLS\n")
        f.write(f"# Regularization: α={best['alpha']}\n")
        f.write(f"# R² (total) = {best['R2']:.6f}\n")
        f.write(f"# Cross contribution = {best['cross_contrib']:+.2f}%\n")
        f.write(f"# Positive: {best['n_pos']}, Negative: {best['n_neg']}\n")
        f.write(f"# Mode_i  Freq_i  Mode_j  Freq_j  d_ij  |d_ij|  Sign\n")
        
        # Sort by absolute value
        sorted_idx = np.argsort(np.abs(best['d_ij']))[::-1]
        for idx in sorted_idx:
            i, j = cross_pairs[idx]
            sign = '+' if best['d_ij'][idx] > 0 else '-'
            f.write(f"{mode_indices[i]:4d}  {mode_freqs[i]:7.2f}  "
                   f"{mode_indices[j]:4d}  {mode_freqs[j]:7.2f}  "
                   f"{best['d_ij'][idx]:12.6e}  {abs(best['d_ij'][idx]):12.6e}  {sign}\n")
    
    print(f"Saved: spectral_coefficients_cross.dat")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nTwo-Stage Fitting Summary:")
    print(f"  Stage 1 (Linear NNLS):    R²={R2_linear:.4f}, capture={100*J_linear_integral/J_total_integral:.1f}%")
    print(f"  Stage 2 (Quadratic NNLS): ΔR²=+{best['R2']-R2_linear:.4f}, capture={best['cross_contrib']:+.1f}%")
    print(f"  Final: R²={best['R2']:.4f}, {best['capture']:.1f}% total spectral capture")
    print(f"\n  Active terms:")
    print(f"    Linear modes: {n_active}/{N_modes}")
    print(f"    Quadratic terms (incl. diagonal): {best['n_pos']}/{n_cross} (all ≥ 0)")
    
    # ========================================================================
    # Peak Decomposition: 520-530 cm⁻¹ and 750 cm⁻¹ regions
    # ========================================================================
    print("\n" + "="*70)
    print("Peak Decomposition Analysis")
    print("="*70)
    
    # Define target regions
    peak1_center = 525.0  # cm⁻¹
    peak1_width = 30.0    # ±30 cm⁻¹ (495-555 cm⁻¹)
    peak2_center = 790.0  # cm⁻¹
    peak2_width = 30.0    # ±30 cm⁻¹ (760-820 cm⁻¹)
    
    # Collect contributions for peak 1 (520-530 cm⁻¹)
    peak1_linear = []
    peak1_cross = []
    
    # Linear contributions to peak 1
    for i in range(N_modes):
        if abs(mode_freqs[i] - peak1_center) < peak1_width and c_k_physical[i] > 1e-10:
            J_mode = c_k_physical[i] * J_k_matrix[:, i]
            peak1_linear.append((i, mode_freqs[i], c_k_physical[i], J_mode))
    
    # Cross term contributions to peak 1
    for idx, (i, j) in enumerate(cross_pairs):
        if best['d_ij'][idx] > 1e-10:
            freq_i = mode_freqs[i]
            freq_j = mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            
            # Check if sum or difference frequency is in peak 1 region
            if abs(sum_freq - peak1_center) < peak1_width or abs(diff_freq - peak1_center) < peak1_width:
                J_cross_single = best['d_ij'][idx] * J_cross_list[idx]
                peak1_cross.append((i, j, freq_i, freq_j, sum_freq, diff_freq, best['d_ij'][idx], J_cross_single))
    
    # Collect contributions for peak 2 (750 cm⁻¹)
    peak2_linear = []
    peak2_cross = []
    
    # Linear contributions to peak 2
    for i in range(N_modes):
        if abs(mode_freqs[i] - peak2_center) < peak2_width and c_k_physical[i] > 1e-10:
            J_mode = c_k_physical[i] * J_k_matrix[:, i]
            peak2_linear.append((i, mode_freqs[i], c_k_physical[i], J_mode))
    
    # Cross term contributions to peak 2
    for idx, (i, j) in enumerate(cross_pairs):
        if best['d_ij'][idx] > 1e-10:
            freq_i = mode_freqs[i]
            freq_j = mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            
            # Check if sum or difference frequency is in peak 2 region
            if abs(sum_freq - peak2_center) < peak2_width or abs(diff_freq - peak2_center) < peak2_width:
                J_cross_single = best['d_ij'][idx] * J_cross_list[idx]
                peak2_cross.append((i, j, freq_i, freq_j, sum_freq, diff_freq, best['d_ij'][idx], J_cross_single))
    
    print(f"\nPeak 1 ({peak1_center:.0f} cm⁻¹ region):")
    print(f"  Linear contributions: {len(peak1_linear)}")
    print(f"  Cross term contributions: {len(peak1_cross)}")
    
    print(f"\nPeak 2 ({peak2_center:.0f} cm⁻¹ region):")
    print(f"  Linear contributions: {len(peak2_linear)}")
    print(f"  Cross term contributions: {len(peak2_cross)}")
    
    # Visualize Peak 1 decomposition
    print(f"\nGenerating peak 1 ({peak1_center:.0f} cm⁻¹) decomposition plot...")
    
    fig1 = plt.figure(figsize=(14, 10))
    gs1 = fig1.add_gridspec(3, 2, hspace=0.35, wspace=0.3)
    
    # Panel 1: Overall comparison
    ax1 = fig1.add_subplot(gs1[0, :])
    mask_region = (freq_grid >= peak1_center - peak1_width) & (freq_grid <= peak1_center + peak1_width)
    ax1.plot(freq_grid[mask_region], J_total[mask_region], 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid[mask_region], best['J_total'][mask_region], 'r--', lw=1.5, label='Full model', alpha=0.8)
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Spectral Density around {peak1_center:.0f} cm⁻¹', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Linear contributions
    ax2 = fig1.add_subplot(gs1[1, 0])
    J_linear_sum = np.zeros_like(freq_grid)
    for i, freq, coeff, J_mode in peak1_linear:
        ax2.plot(freq_grid[mask_region], J_mode[mask_region], '-', alpha=0.6, lw=1, 
                label=f'Mode {mode_indices[i]} ({freq:.1f} cm⁻¹)')
        J_linear_sum += J_mode
    ax2.plot(freq_grid[mask_region], J_linear_sum[mask_region], 'b-', lw=2, alpha=0.8, label='Sum')
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax2.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax2.set_title(f'Linear Contributions ({len(peak1_linear)} modes)', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Cross term contributions
    ax3 = fig1.add_subplot(gs1[1, 1])
    J_cross_sum = np.zeros_like(freq_grid)
    colors_cross = plt.cm.tab10(np.linspace(0, 1, min(10, len(peak1_cross))))
    for idx, (i, j, fi, fj, sum_f, diff_f, coeff, J_cross_single) in enumerate(peak1_cross[:10]):
        peak_type = "Σ" if abs(sum_f - peak1_center) < abs(diff_f - peak1_center) else "Δ"
        ax3.plot(freq_grid[mask_region], J_cross_single[mask_region], '-', alpha=0.6, lw=1,
                color=colors_cross[idx % 10],
                label=f'{mode_indices[i]}×{mode_indices[j]} ({peak_type}={sum_f if peak_type=="Σ" else diff_f:.0f})')
        J_cross_sum += J_cross_single
    ax3.plot(freq_grid[mask_region], J_cross_sum[mask_region], 'r-', lw=2, alpha=0.8, label='Sum')
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax3.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax3.set_title(f'Quadratic Contributions (top {min(10, len(peak1_cross))})', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=7, loc='best', ncol=2)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Stacked contributions
    ax4 = fig1.add_subplot(gs1[2, :])
    ax4.plot(freq_grid[mask_region], J_total[mask_region], 'k-', lw=2, label='J_total (data)', alpha=0.8, zorder=3)
    ax4.fill_between(freq_grid[mask_region], 0, J_linear_sum[mask_region], 
                     alpha=0.5, color='blue', label='Linear', zorder=1)
    ax4.fill_between(freq_grid[mask_region], J_linear_sum[mask_region], 
                     J_linear_sum[mask_region] + J_cross_sum[mask_region],
                     alpha=0.5, color='red', label='Quadratic', zorder=2)
    ax4.plot(freq_grid[mask_region], (J_linear_sum + J_cross_sum)[mask_region], 
            'g--', lw=1.5, label='Total model', alpha=0.8, zorder=3)
    ax4.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax4.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax4.set_title(f'Decomposition at {peak1_center:.0f} cm⁻¹', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.savefig(f'peak_decomposition_{int(peak1_center)}cm.png', dpi=150, bbox_inches='tight')
    print(f"Saved: peak_decomposition_{int(peak1_center)}cm.png")
    
    # Visualize Peak 2 decomposition
    print(f"Generating peak 2 ({peak2_center:.0f} cm⁻¹) decomposition plot...")
    
    fig2 = plt.figure(figsize=(14, 10))
    gs2 = fig2.add_gridspec(3, 2, hspace=0.35, wspace=0.3)
    
    # Panel 1: Overall comparison
    ax1 = fig2.add_subplot(gs2[0, :])
    mask_region2 = (freq_grid >= peak2_center - peak2_width) & (freq_grid <= peak2_center + peak2_width)
    ax1.plot(freq_grid[mask_region2], J_total[mask_region2], 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid[mask_region2], best['J_total'][mask_region2], 'r--', lw=1.5, label='Full model', alpha=0.8)
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Spectral Density around {peak2_center:.0f} cm⁻¹', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Linear contributions
    ax2 = fig2.add_subplot(gs2[1, 0])
    J_linear_sum2 = np.zeros_like(freq_grid)
    for i, freq, coeff, J_mode in peak2_linear:
        ax2.plot(freq_grid[mask_region2], J_mode[mask_region2], '-', alpha=0.6, lw=1,
                label=f'Mode {mode_indices[i]} ({freq:.1f} cm⁻¹)')
        J_linear_sum2 += J_mode
    ax2.plot(freq_grid[mask_region2], J_linear_sum2[mask_region2], 'b-', lw=2, alpha=0.8, label='Sum')
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax2.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax2.set_title(f'Linear Contributions ({len(peak2_linear)} modes)', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Cross term contributions
    ax3 = fig2.add_subplot(gs2[1, 1])
    J_cross_sum2 = np.zeros_like(freq_grid)
    colors_cross2 = plt.cm.tab10(np.linspace(0, 1, min(10, len(peak2_cross))))
    for idx, (i, j, fi, fj, sum_f, diff_f, coeff, J_cross_single) in enumerate(peak2_cross[:10]):
        peak_type = "Σ" if abs(sum_f - peak2_center) < abs(diff_f - peak2_center) else "Δ"
        ax3.plot(freq_grid[mask_region2], J_cross_single[mask_region2], '-', alpha=0.6, lw=1,
                color=colors_cross2[idx % 10],
                label=f'{mode_indices[i]}×{mode_indices[j]} ({peak_type}={sum_f if peak_type=="Σ" else diff_f:.0f})')
        J_cross_sum2 += J_cross_single
    ax3.plot(freq_grid[mask_region2], J_cross_sum2[mask_region2], 'r-', lw=2, alpha=0.8, label='Sum')
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax3.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax3.set_title(f'Quadratic Contributions (top {min(10, len(peak2_cross))})', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=7, loc='best', ncol=2)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Stacked contributions
    ax4 = fig2.add_subplot(gs2[2, :])
    ax4.plot(freq_grid[mask_region2], J_total[mask_region2], 'k-', lw=2, label='J_total (data)', alpha=0.8, zorder=3)
    ax4.fill_between(freq_grid[mask_region2], 0, J_linear_sum2[mask_region2],
                     alpha=0.5, color='blue', label='Linear', zorder=1)
    ax4.fill_between(freq_grid[mask_region2], J_linear_sum2[mask_region2],
                     J_linear_sum2[mask_region2] + J_cross_sum2[mask_region2],
                     alpha=0.5, color='red', label='Quadratic', zorder=2)
    ax4.plot(freq_grid[mask_region2], (J_linear_sum2 + J_cross_sum2)[mask_region2],
            'g--', lw=1.5, label='Total model', alpha=0.8, zorder=3)
    ax4.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax4.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax4.set_title(f'Decomposition at {peak2_center:.0f} cm⁻¹', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.savefig(f'peak_decomposition_{int(peak2_center)}cm.png', dpi=150, bbox_inches='tight')
    print(f"Saved: peak_decomposition_{int(peak2_center)}cm.png")
    
    # Print detailed contributions
    print("\n" + "="*70)
    print(f"Detailed contributions to {peak1_center:.0f} cm⁻¹ peak:")
    print("="*70)
    
    if peak1_linear:
        print(f"\nLinear terms ({len(peak1_linear)}):")
        for i, freq, coeff, J_mode in sorted(peak1_linear, key=lambda x: x[2], reverse=True):
            integral = np.trapezoid(J_mode, freq_grid)
            print(f"  Mode {mode_indices[i]:2d} ({freq:6.1f} cm⁻¹): c_k = {coeff:.6e}, ∫J = {integral:.6e}")
    
    if peak1_cross:
        print(f"\nQuadratic terms (top 10 of {len(peak1_cross)}):")
        peak1_cross_sorted = sorted(peak1_cross, key=lambda x: x[6], reverse=True)
        for i, j, fi, fj, sum_f, diff_f, coeff, J_cross_single in peak1_cross_sorted[:10]:
            integral = np.trapezoid(J_cross_single, freq_grid)
            peak_type = "Sum" if abs(sum_f - peak1_center) < abs(diff_f - peak1_center) else "Diff"
            target_f = sum_f if peak_type == "Sum" else diff_f
            print(f"  Mode {mode_indices[i]:2d}×{mode_indices[j]:2d} ({fi:6.1f}×{fj:6.1f}): "
                  f"{peak_type}={target_f:6.1f}, d_ij = {coeff:.6e}, ∫J = {integral:.6e}")
    
    print("\n" + "="*70)
    print(f"Detailed contributions to {peak2_center:.0f} cm⁻¹ peak:")
    print("="*70)
    
    if peak2_linear:
        print(f"\nLinear terms ({len(peak2_linear)}):")
        for i, freq, coeff, J_mode in sorted(peak2_linear, key=lambda x: x[2], reverse=True):
            integral = np.trapezoid(J_mode, freq_grid)
            print(f"  Mode {mode_indices[i]:2d} ({freq:6.1f} cm⁻¹): c_k = {coeff:.6e}, ∫J = {integral:.6e}")
    
    if peak2_cross:
        print(f"\nQuadratic terms (top 10 of {len(peak2_cross)}):")
        peak2_cross_sorted = sorted(peak2_cross, key=lambda x: x[6], reverse=True)
        for i, j, fi, fj, sum_f, diff_f, coeff, J_cross_single in peak2_cross_sorted[:10]:
            integral = np.trapezoid(J_cross_single, freq_grid)
            peak_type = "Sum" if abs(sum_f - peak2_center) < abs(diff_f - peak2_center) else "Diff"
            target_f = sum_f if peak_type == "Sum" else diff_f
            print(f"  Mode {mode_indices[i]:2d}×{mode_indices[j]:2d} ({fi:6.1f}×{fj:6.1f}): "
                  f"{peak_type}={target_f:6.1f}, d_ij = {coeff:.6e}, ∫J = {integral:.6e}")

else:
    print("\nNot enough modes in frequency range for cross terms!")

