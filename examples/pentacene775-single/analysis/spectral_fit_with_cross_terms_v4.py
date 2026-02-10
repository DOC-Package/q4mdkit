#!/usr/bin/env python3
"""
Two-stage spectral domain fitting:
1. Linear terms with NNLS (c_k ≥ 0)
2. 3rd order correlations on residual (κ_k d_kl terms)

3rd order correlations:
  J^(3)_kl(ω) combines:
    Type 1: FT[<q_k(t)q_l(t)q_k(0)>] with coeff κ_k d_kl
    Type 2: FT[<q_k(t)q_l(t)q_l(0)>] with coeff κ_l d_kl
    Type 3: FT[<q_k(t)q_k(0)q_l(0)>] with coeff κ_l d_kl
    Type 4: FT[<q_l(t)q_l(0)q_k(0)>] with coeff κ_k d_kl

where κ_k are the linear coefficients from Stage 1.
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import nnls, least_squares
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
print("Two-Stage Spectral Fitting v4: Linear (NNLS) + 3rd Order Correlations")
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

# ============================================================================
# Helper functions
# ============================================================================
def compute_autocorr(x):
    """Compute autocorrelation with Hanning window"""
    x_mean = np.mean(x)
    x_centered = x - x_mean
    N = len(x)
    
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(x_centered[:N-lag] * x_centered[lag:])
    
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return acf * window

def compute_spectral_density(x, freq):
    """Compute spectral density using Simpson integration"""
    C = compute_autocorr(x)
    t = np.arange(max_lag) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

def compute_spectral_from_corr(corr, freq):
    """Compute spectral density from pre-computed correlation"""
    t = np.arange(len(corr)) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = corr * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

# Frequency grid
freq_grid = np.linspace(420, 1800, 1381)
print(f"\nFrequency grid: {len(freq_grid)} points (420-1800 cm⁻¹)")

# ============================================================================
# Load cached linear spectral densities
# ============================================================================
cache_file_linear = 'spectral_cache_linear.npz'

if os.path.exists(cache_file_linear):
    print(f"\nLoading cached linear spectral densities from {cache_file_linear}...")
    cache_data = np.load(cache_file_linear)
    J_total = cache_data['J_total']
    J_k_matrix = cache_data['J_k_matrix']
    J_k_integrals = cache_data['J_k_integrals']
    cached_freq_grid = cache_data['freq_grid']
    
    if len(cached_freq_grid) == len(freq_grid) and np.allclose(cached_freq_grid, freq_grid):
        print(f"  Loaded J_total and J_k_matrix ({J_k_matrix.shape[1]} modes)")
    else:
        print(f"  Cache frequency grid mismatch, recomputing...")
        os.remove(cache_file_linear)
        cache_data = None
else:
    cache_data = None

if cache_data is None or not os.path.exists(cache_file_linear):
    print("\nComputing J_total...")
    J_total = compute_spectral_density(energy_diff, freq_grid)
    
    print("\nComputing J_k for all modes...")
    J_k_matrix = np.zeros((len(freq_grid), N_modes))
    J_k_integrals = np.zeros(N_modes)
    for k in range(N_modes):
        if (k+1) % 10 == 0 or k == 0:
            print(f"  Mode {k+1}/{N_modes}")
        J_k_matrix[:, k] = compute_spectral_density(mode_coords[:, k], freq_grid)
        J_k_integrals[k] = np.trapezoid(J_k_matrix[:, k], freq_grid)
    
    np.savez_compressed(cache_file_linear,
                        J_total=J_total,
                        J_k_matrix=J_k_matrix,
                        J_k_integrals=J_k_integrals,
                        freq_grid=freq_grid)
    print(f"\n  Saved linear spectral densities to {cache_file_linear}")

J_total_integral = np.trapezoid(J_total, freq_grid)
print(f"  ∫J_total = {J_total_integral:.6e} cm⁻¹²")

# ============================================================================
# Stage 1: Linear terms with NNLS
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, c_k ≥ 0)")
print(f"Using {np.sum(stage1_mask)} modes with 450 ≤ freq ≤ 1800 cm⁻¹")
print("="*70)

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

# These are the κ_k for 3rd order correlations
kappa_k = c_k_physical.copy()

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

n_active_linear = np.sum(c_k_physical > 1e-10)

print(f"\nLinear model results:")
print(f"  R² = {R2_linear:.4f}")
print(f"  ∫J_linear = {J_linear_integral:.6e} cm⁻¹² ({100*J_linear_integral/J_total_integral:.1f}%)")
print(f"  ∫J_residual = {J_residual_1_integral:.6e} cm⁻¹² ({100*J_residual_1_integral/J_total_integral:.1f}%)")
print(f"  Active modes: {n_active_linear}/{N_modes}")

# Top 10 linear modes (these will be used for κ_k)
top_linear_idx = np.argsort(kappa_k)[-15:][::-1]
print(f"\nTop 15 linear modes (κ_k):")
for rank, idx in enumerate(top_linear_idx, 1):
    print(f"  {rank:2d}. Mode {mode_indices[idx]:2d} ({mode_freqs[idx]:6.1f} cm⁻¹): κ_k = {kappa_k[idx]:.6e}")

# ============================================================================
# Stage 2: 3rd Order Correlations
# ============================================================================
print("\n" + "="*70)
print("Stage 2: 3rd Order Correlations on Residual")
print("="*70)

# Select modes for 3rd order correlations
# Extended mode selection: frequency range 0-1600 cm⁻¹
target_freqs = [524, 800, 1092, 1231]  # cm⁻¹
tolerance = 50.0

# Extended mode selection: 0-1600 cm⁻¹
freq_min_cross, freq_max_cross = 0, 1600
mask_cross = (mode_freqs >= freq_min_cross) & (mode_freqs <= freq_max_cross)

# Use frequency-based selection
idx_cross = np.where(mask_cross)[0]
n_cross_modes = len(idx_cross)

print(f"\nModes selected for 3rd order: {n_cross_modes}")
print(f"  Frequency range {freq_min_cross}-{freq_max_cross} cm⁻¹: {np.sum(mask_cross)}")

# Build mode pairs
print(f"\nBuilding mode pairs...")
mode_pairs = []
pair_info = {}

for i_idx, i in enumerate(idx_cross):
    for j_idx, j in enumerate(idx_cross):
        if j <= i:
            continue
        
        freq_i = mode_freqs[i]
        freq_j = mode_freqs[j]
        sum_freq = freq_i + freq_j
        diff_freq = abs(freq_i - freq_j)
        
        # Check target matching
        target_matched = 0
        match_type = 'none'
        for target in target_freqs:
            if abs(sum_freq - target) <= tolerance:
                target_matched = target
                match_type = 'sum'
                break
            elif abs(diff_freq - target) <= tolerance:
                target_matched = target
                match_type = 'diff'
                break
        
        mode_pairs.append((i, j))
        pair_info[(i, j)] = {
            'mode_i': mode_indices[i], 'mode_j': mode_indices[j],
            'freq_i': freq_i, 'freq_j': freq_j,
            'kappa_i': kappa_k[i], 'kappa_j': kappa_k[j],
            'target': target_matched, 'type': match_type
        }

N_pairs = len(mode_pairs)
print(f"Total pairs: {N_pairs}")

# Count pairs with non-zero kappa
pairs_with_kappa = sum(1 for (i,j) in mode_pairs if kappa_k[i] > 0 or kappa_k[j] > 0)
print(f"Pairs with at least one non-zero κ: {pairs_with_kappa}")

# Target frequency coverage
print(f"\nTarget frequency coverage:")
for target in target_freqs:
    count = sum(1 for p, info in pair_info.items() if info['target'] == target)
    print(f"  {target} cm⁻¹: {count} pairs")

# ============================================================================
# Compute 3rd order correlation spectral densities
# ============================================================================
cache_file_3rd = 'spectral_cache_3rd_order_v4.npz'

# First, try to load from cache and filter to v3-style pairs
cache_3rd_loaded = False
if os.path.exists(cache_file_3rd):
    print(f"\nLoading cached 3rd order spectral densities from {cache_file_3rd}...")
    cache_3rd = np.load(cache_file_3rd, allow_pickle=True)
    
    cached_pairs = [tuple(p) for p in cache_3rd['pairs']]
    
    # Filter cached data to only include v3-style pairs
    print(f"  Total cached pairs: {len(cached_pairs)}")
    print(f"  Filtering to v3-style mode selection (0-800 cm⁻¹)...")
    
    # Find which cached pairs match our v3-style selection
    pair_to_cache_idx = {p: idx for idx, p in enumerate(cached_pairs)}
    valid_pairs = []
    valid_cache_idx = []
    
    for pair in mode_pairs:
        if pair in pair_to_cache_idx:
            valid_pairs.append(pair)
            valid_cache_idx.append(pair_to_cache_idx[pair])
    
    if len(valid_pairs) > 0:
        # Extract only the needed data
        J_3rd_type1 = cache_3rd['J_3rd_type1'][valid_cache_idx, :]
        J_3rd_type2 = cache_3rd['J_3rd_type2'][valid_cache_idx, :]
        J_3rd_type3 = cache_3rd['J_3rd_type3'][valid_cache_idx, :]
        J_3rd_type4 = cache_3rd['J_3rd_type4'][valid_cache_idx, :]
        
        # Update mode_pairs and pair_info to match filtered data
        filtered_pair_info = {p: pair_info[p] for p in valid_pairs}
        mode_pairs = valid_pairs
        pair_info = filtered_pair_info
        N_pairs = len(mode_pairs)
        
        print(f"  Filtered to {N_pairs} v3-style pairs")
        cache_3rd_loaded = True
    else:
        print(f"  No matching pairs found in cache, will recompute...")
        cache_3rd_loaded = False
else:
    cache_3rd_loaded = False

if not cache_3rd_loaded:
    print(f"\nComputing 3rd order correlations for {N_pairs} pairs...")
    print(f"  Type 1: <q_k(t)q_l(t)q_k(0)>")
    print(f"  Type 2: <q_k(t)q_l(t)q_l(0)>")
    print(f"  Type 3: <q_k(t)q_k(0)q_l(0)>")
    print(f"  Type 4: <q_l(t)q_l(0)q_k(0)>")
    
    J_3rd_type1 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type2 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type3 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type4 = np.zeros((N_pairs, len(freq_grid)))
    
    t_arr = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag))
    
    for pair_idx, (i, j) in enumerate(mode_pairs):
        q_k = mode_coords[:, i]
        q_l = mode_coords[:, j]
        
        # Center
        q_k_c = q_k - np.mean(q_k)
        q_l_c = q_l - np.mean(q_l)
        
        # Product at same time
        product_kl = q_k_c * q_l_c
        
        # Type 1: <q_k(t)q_l(t)q_k(0)>
        corr_type1 = np.zeros(max_lag)
        for lag in range(max_lag):
            corr_type1[lag] = np.mean(product_kl[lag:] * q_k_c[:N_frames-lag])
        corr_type1 *= window
        J_3rd_type1[pair_idx, :] = compute_spectral_from_corr(corr_type1, freq_grid)
        
        # Type 2: <q_k(t)q_l(t)q_l(0)>
        corr_type2 = np.zeros(max_lag)
        for lag in range(max_lag):
            corr_type2[lag] = np.mean(product_kl[lag:] * q_l_c[:N_frames-lag])
        corr_type2 *= window
        J_3rd_type2[pair_idx, :] = compute_spectral_from_corr(corr_type2, freq_grid)
        
        # Type 3: <q_k(t)q_k(0)q_l(0)>
        corr_type3 = np.zeros(max_lag)
        for lag in range(max_lag):
            corr_type3[lag] = np.mean(q_k_c[lag:] * product_kl[:N_frames-lag])
        corr_type3 *= window
        J_3rd_type3[pair_idx, :] = compute_spectral_from_corr(corr_type3, freq_grid)
        
        # Type 4: <q_l(t)q_l(0)q_k(0)>
        corr_type4 = np.zeros(max_lag)
        for lag in range(max_lag):
            corr_type4[lag] = np.mean(q_l_c[lag:] * product_kl[:N_frames-lag])
        corr_type4 *= window
        J_3rd_type4[pair_idx, :] = compute_spectral_from_corr(corr_type4, freq_grid)
        
        if (pair_idx + 1) % 100 == 0 or pair_idx == 0:
            print(f"    Progress: {pair_idx+1}/{N_pairs} pairs")
    
    # Save cache
    np.savez_compressed(cache_file_3rd,
                        J_3rd_type1=J_3rd_type1,
                        J_3rd_type2=J_3rd_type2,
                        J_3rd_type3=J_3rd_type3,
                        J_3rd_type4=J_3rd_type4,
                        pairs=np.array(mode_pairs))
    print(f"\n  Saved to {cache_file_3rd}")
    print(f"  File size: {os.path.getsize(cache_file_3rd) / 1024**2:.1f} MB")

print(f"\n3rd order correlation statistics:")
print(f"  J_3rd_type1 RMS: {np.sqrt(np.mean(J_3rd_type1**2)):.6e}")
print(f"  J_3rd_type2 RMS: {np.sqrt(np.mean(J_3rd_type2**2)):.6e}")
print(f"  J_3rd_type3 RMS: {np.sqrt(np.mean(J_3rd_type3**2)):.6e}")
print(f"  J_3rd_type4 RMS: {np.sqrt(np.mean(J_3rd_type4**2)):.6e}")

# ============================================================================
# Build design matrix for 3rd order fitting
# ============================================================================
print("\n" + "="*70)
print("Building 3rd order design matrix")
print("="*70)

# The 3rd order contribution for pair (k,l) with coefficient d_kl is:
#   κ_k d_kl * J_type1 + κ_l d_kl * J_type2 + κ_l d_kl * J_type3 + κ_k d_kl * J_type4
# = d_kl * (κ_k * (J_type1 + J_type4) + κ_l * (J_type2 + J_type3))

# Get kappa values for each pair
kappa_i = np.array([kappa_k[i] for i, j in mode_pairs])
kappa_j = np.array([kappa_k[j] for i, j in mode_pairs])

# Combined 3rd order basis: κ_k*(Type1+Type4) + κ_l*(Type2+Type3)
J_3rd_combined = (kappa_i[:, np.newaxis] * (J_3rd_type1 + J_3rd_type4) + 
                  kappa_j[:, np.newaxis] * (J_3rd_type2 + J_3rd_type3))

print(f"Design matrix shape: {J_3rd_combined.T.shape}")
print(f"J_3rd_combined RMS: {np.sqrt(np.mean(J_3rd_combined**2)):.6e}")

# Check how many pairs have non-zero contribution
nonzero_pairs = np.sum(np.any(np.abs(J_3rd_combined) > 1e-15, axis=1))
print(f"Pairs with non-zero 3rd order contribution: {nonzero_pairs}/{N_pairs}")

# ============================================================================
# Fit 3rd order terms to residual
# ============================================================================
print("\n" + "="*70)
print("Fitting 3rd order terms to residual")
print("="*70)

# Normalize for numerical stability
residual_rms = np.sqrt(np.mean(J_residual_1**2))
J_3rd_norm = J_3rd_combined / residual_rms
residual_norm = J_residual_1 / residual_rms

print(f"Normalization factor: {residual_rms:.6e}")

# Method 1: Regularized least squares
print("\nMethod 1: Ridge regression...")
from sklearn.linear_model import Ridge

alphas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1]
best_R2 = -np.inf
best_alpha = None
best_d_kl = None

for alpha in alphas:
    ridge = Ridge(alpha=alpha, fit_intercept=False)
    ridge.fit(J_3rd_norm.T, residual_norm)
    d_kl_ridge = ridge.coef_
    
    J_3rd_model = J_3rd_combined.T @ d_kl_ridge * residual_rms
    residual_after = J_residual_1 - J_3rd_model
    
    ss_res = np.sum(residual_after**2)
    R2 = 1 - ss_res / ss_tot
    
    if R2 > best_R2:
        best_R2 = R2
        best_alpha = alpha
        best_d_kl = d_kl_ridge
    
    print(f"  α={alpha:.0e}: R²={R2:.4f}, ||d_kl||={np.linalg.norm(d_kl_ridge):.3e}")

print(f"\nBest: α={best_alpha:.0e}, R²={best_R2:.4f}")

# Use best result
d_kl = best_d_kl

# Compute final 3rd order model
J_3rd_model = J_3rd_combined.T @ d_kl

# Total model
J_total_model = J_linear + J_3rd_model
J_residual_2 = J_total - J_total_model

# Final statistics
ss_res_2 = np.sum(J_residual_2**2)
R2_total = 1 - ss_res_2 / ss_tot

J_3rd_integral = np.trapezoid(J_3rd_model, freq_grid)
J_total_model_integral = np.trapezoid(J_total_model, freq_grid)
J_residual_2_integral = np.trapezoid(J_residual_2, freq_grid)

print(f"\n" + "="*70)
print("Fit Quality")
print("="*70)
print(f"\n  Stage 1 (Linear): R² = {R2_linear:.4f}")
print(f"  Stage 2 (3rd order): R² = {R2_total:.4f}")
print(f"  Improvement: ΔR² = {R2_total - R2_linear:.4f}")

print(f"\n  Spectral capture:")
print(f"    Linear: {100*J_linear_integral/J_total_integral:.1f}%")
print(f"    3rd order: {100*J_3rd_integral/J_total_integral:+.1f}%")
print(f"    Total: {100*J_total_model_integral/J_total_integral:.1f}%")

print(f"\n  Residual:")
print(f"    After Stage 1: {np.sqrt(np.mean(J_residual_1**2)):.6e} (RMS)")
print(f"    After Stage 2: {np.sqrt(np.mean(J_residual_2**2)):.6e} (RMS)")
print(f"    Reduction: {100*(1 - np.sqrt(np.mean(J_residual_2**2))/np.sqrt(np.mean(J_residual_1**2))):.1f}%")

# ============================================================================
# Fit quality at target frequencies
# ============================================================================
print("\n" + "="*70)
print("Fit Quality at Target Frequencies")
print("="*70)

target_tolerance = 30.0
print(f"\n{'Target':>8} {'J_total':>12} {'J_model':>12} {'J_residual':>12} {'Capture':>10}")
print(f"{'-'*60}")

for target in target_freqs:
    freq_window = (freq_grid >= target - target_tolerance) & (freq_grid <= target + target_tolerance)
    
    J_total_at_target = np.mean(J_total[freq_window])
    J_model_at_target = np.mean(J_total_model[freq_window])
    J_residual_at_target = np.mean(J_residual_2[freq_window])
    capture_pct = 100 * J_model_at_target / J_total_at_target if J_total_at_target > 0 else 0
    
    print(f"{target:8.0f} {J_total_at_target:12.3e} {J_model_at_target:12.3e} "
          f"{J_residual_at_target:12.3e} {capture_pct:9.1f}%")

# ============================================================================
# Top mode pairs
# ============================================================================
print("\n" + "="*70)
print("Top 20 mode pairs by |d_kl|")
print("="*70)

top_idx = np.argsort(np.abs(d_kl))[-20:][::-1]
print(f"\n{'Rank':>4} {'Mode k':>7} {'Mode l':>7} {'ωk':>7} {'ωl':>7} {'κk':>10} {'κl':>10} {'d_kl':>12} {'Target':>8}")
print("-"*95)

for rank, idx in enumerate(top_idx, 1):
    i, j = mode_pairs[idx]
    info = pair_info[(i, j)]
    sign = '+' if d_kl[idx] > 0 else '-'
    print(f"{rank:4d} {info['mode_i']:7d} {info['mode_j']:7d} {info['freq_i']:7.1f} {info['freq_j']:7.1f} "
          f"{info['kappa_i']:10.2e} {info['kappa_j']:10.2e} {sign}{abs(d_kl[idx]):11.3e} {info['target']:8d}")

# ============================================================================
# Analyze 3rd order contributions by target frequency
# ============================================================================
print("\n" + "="*70)
print("3rd Order Contributions by Target Frequency")
print("="*70)

for target in target_freqs:
    freq_window = (freq_grid >= target - target_tolerance) & (freq_grid <= target + target_tolerance)
    
    print(f"\n  Target {target} cm⁻¹:")
    
    contributions = []
    for pair_idx, (i, j) in enumerate(mode_pairs):
        info = pair_info[(i, j)]
        if info['target'] == target and abs(d_kl[pair_idx]) > 1e-15:
            # Contribution at target frequency
            J_contrib = d_kl[pair_idx] * J_3rd_combined[pair_idx, :]
            J_at_target = np.mean(J_contrib[freq_window])
            
            contributions.append({
                'pair_idx': pair_idx,
                'mode_i': info['mode_i'], 'mode_j': info['mode_j'],
                'freq_i': info['freq_i'], 'freq_j': info['freq_j'],
                'kappa_i': info['kappa_i'], 'kappa_j': info['kappa_j'],
                'd_kl': d_kl[pair_idx],
                'J_at_target': J_at_target,
                'type': info['type']
            })
    
    # Sort by absolute contribution
    contributions_sorted = sorted(contributions, key=lambda x: abs(x['J_at_target']), reverse=True)
    
    if contributions_sorted:
        print(f"  {'Rank':>4} {'Modes':>10} {'κk':>10} {'κl':>10} {'d_kl':>12} {'J(target)':>12} {'Type':>6}")
        print(f"  {'-'*75}")
        for rank, c in enumerate(contributions_sorted[:10], 1):
            sign = '+' if c['J_at_target'] > 0 else '-'
            print(f"  {rank:4d} {c['mode_i']:4d}×{c['mode_j']:<4d} {c['kappa_i']:10.2e} {c['kappa_j']:10.2e} "
                  f"{c['d_kl']:12.3e} {sign}{abs(c['J_at_target']):11.3e} {c['type']:>6}")
    else:
        print(f"    No significant contributions")

# ============================================================================
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Generating visualization...")
print("="*70)

fig, axes = plt.subplots(3, 2, figsize=(14, 12))

# Panel 1: Overall fit
ax1 = axes[0, 0]
ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
ax1.plot(freq_grid, J_linear, 'b--', lw=1, label=f'Stage 1: Linear ({100*J_linear_integral/J_total_integral:.1f}%)', alpha=0.6)
ax1.plot(freq_grid, J_total_model, 'r-', lw=1.5, label=f'Stage 2: +3rd order (R²={R2_total:.4f})', alpha=0.8)
ax1.set_ylabel('J(ν) (cm⁻¹²)')
ax1.set_title('Two-Stage Spectral Fitting: Linear + 3rd Order', fontweight='bold')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Panel 2: Residual comparison
ax2 = axes[0, 1]
ax2.plot(freq_grid, J_residual_1, 'b-', lw=1, label='After Stage 1', alpha=0.6)
ax2.plot(freq_grid, J_residual_2, 'r-', lw=1.5, label='After Stage 2 (3rd order)', alpha=0.8)
ax2.axhline(0, color='k', ls='--', lw=0.5)
for target in target_freqs:
    ax2.axvline(target, color='gray', ls=':', alpha=0.5)
ax2.set_ylabel('Residual (cm⁻¹²)')
ax2.set_title('Residual Comparison', fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# Panel 3: 3rd order contribution
ax3 = axes[1, 0]
ax3.plot(freq_grid, J_3rd_model, 'purple', lw=1.5, alpha=0.8)
ax3.axhline(0, color='k', ls='--', lw=0.5)
for target in target_freqs:
    ax3.axvline(target, color='red', ls=':', alpha=0.5)
ax3.set_ylabel('J(ν) (cm⁻¹²)')
ax3.set_title(f'3rd Order Contribution ({100*J_3rd_integral/J_total_integral:+.1f}%)', fontweight='bold')
ax3.grid(True, alpha=0.3)

# Panel 4: d_kl distribution
ax4 = axes[1, 1]
ax4.hist(d_kl, bins=50, color='purple', alpha=0.7, edgecolor='black')
ax4.axvline(0, color='k', ls='--', lw=0.5)
ax4.set_xlabel('d_kl')
ax4.set_ylabel('Count')
ax4.set_title(f'd_kl Distribution ({N_pairs} pairs)', fontweight='bold')
ax4.grid(True, alpha=0.3)

# Panel 5: Top linear modes (κ_k)
ax5 = axes[2, 0]
top_lin = np.argsort(kappa_k)[-10:][::-1]
y_pos = np.arange(10)
ax5.barh(y_pos, kappa_k[top_lin], color='blue', alpha=0.7)
ax5.set_yticks(y_pos)
ax5.set_yticklabels([f"Mode {mode_indices[i]} ({mode_freqs[i]:.0f})" for i in top_lin], fontsize=8)
ax5.set_xlabel('κ_k (linear coefficient)')
ax5.set_title('Top 10 Linear Modes', fontweight='bold')
ax5.grid(True, alpha=0.3, axis='x')

# Panel 6: Top d_kl coefficients
ax6 = axes[2, 1]
colors = ['green' if d > 0 else 'red' for d in d_kl[top_idx[:10]]]
ax6.barh(np.arange(10), d_kl[top_idx[:10]], color=colors, alpha=0.7)
ax6.set_yticks(np.arange(10))
labels = []
for idx in top_idx[:10]:
    info = pair_info[mode_pairs[idx]]
    labels.append(f"{info['mode_i']}×{info['mode_j']}")
ax6.set_yticklabels(labels, fontsize=9)
ax6.axvline(0, color='k', ls='-', lw=0.5)
ax6.set_xlabel('d_kl')
ax6.set_title('Top 10 Mode Pairs (3rd order)', fontweight='bold')
ax6.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('spectral_fit_3rd_order_v4.png', dpi=150, bbox_inches='tight')
print(f"Saved: spectral_fit_3rd_order_v4.png")
plt.close()

# ============================================================================
# Save results
# ============================================================================
print("\nSaving results...")

# Save d_kl coefficients
with open('d_kl_3rd_order_v4.dat', 'w') as f:
    f.write("# 3rd order coupling coefficients d_kl\n")
    f.write(f"# Ridge regularization α = {best_alpha}\n")
    f.write(f"# R² (total) = {R2_total:.6f}\n")
    f.write(f"# 3rd order contribution = {100*J_3rd_integral/J_total_integral:+.2f}%\n")
    f.write("#\n")
    f.write("# mode_k  mode_l  freq_k  freq_l  kappa_k  kappa_l  d_kl  target  type\n")
    
    sorted_idx = np.argsort(np.abs(d_kl))[::-1]
    for idx in sorted_idx:
        i, j = mode_pairs[idx]
        info = pair_info[(i, j)]
        f.write(f"{info['mode_i']:4d}  {info['mode_j']:4d}  {info['freq_i']:8.2f}  {info['freq_j']:8.2f}  "
                f"{info['kappa_i']:12.6e}  {info['kappa_j']:12.6e}  {d_kl[idx]:14.6e}  "
                f"{info['target']:6d}  {info['type']}\n")

print(f"Saved: d_kl_3rd_order_v4.dat")

# Save spectral densities
with open('spectral_density_3rd_order_v4.dat', 'w') as f:
    f.write("# Frequency(cm⁻¹)  J_total  J_linear  J_3rd  J_model_total  Residual_1  Residual_2\n")
    for i in range(len(freq_grid)):
        f.write(f"{freq_grid[i]:10.4f}  {J_total[i]:14.6e}  {J_linear[i]:14.6e}  "
                f"{J_3rd_model[i]:14.6e}  {J_total_model[i]:14.6e}  "
                f"{J_residual_1[i]:14.6e}  {J_residual_2[i]:14.6e}\n")

print(f"Saved: spectral_density_3rd_order_v4.dat")

# ============================================================================
# Peak Decomposition Analysis
# ============================================================================
print("\n" + "="*70)
print("Peak Decomposition Analysis")
print("="*70)

# Define target regions
peak_targets = [
    {'center': 524.0, 'width': 40.0, 'name': '524'},
    {'center': 800.0, 'width': 40.0, 'name': '800'},
    {'center': 1092.0, 'width': 40.0, 'name': '1092'},
    {'center': 1231.0, 'width': 40.0, 'name': '1231'},
]

# Store all peak data
all_peak_data = {}

for peak_info in peak_targets:
    peak_center = peak_info['center']
    peak_width = peak_info['width']
    peak_name = peak_info['name']
    
    # Collect linear contributions (modes directly in peak region)
    peak_linear = []
    for i in range(N_modes):
        if abs(mode_freqs[i] - peak_center) < peak_width and c_k_physical[i] > 1e-10:
            J_mode = c_k_physical[i] * J_k_matrix[:, i]
            peak_linear.append((i, mode_freqs[i], c_k_physical[i], J_mode))
    
    # Collect 3rd order contributions - ALL pairs that contribute at this peak
    peak_3rd = []
    for pair_idx, (i, j) in enumerate(mode_pairs):
        if abs(d_kl[pair_idx]) > 1e-10:
            freq_i = mode_freqs[i]
            freq_j = mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            
            # 3rd order contribution for this pair
            J_3rd_single = d_kl[pair_idx] * J_3rd_combined[pair_idx, :]
            
            # Calculate actual contribution at this peak region
            peak_mask = (freq_grid >= peak_center - peak_width) & (freq_grid <= peak_center + peak_width)
            J_at_peak = np.mean(J_3rd_single[peak_mask])
            
            # Include if it has significant contribution at this peak (not based on sum/diff frequency)
            if abs(J_at_peak) > 1e-12:
                peak_3rd.append({
                    'pair_idx': pair_idx,
                    'i': i, 'j': j,
                    'mode_i': mode_indices[i], 'mode_j': mode_indices[j],
                    'freq_i': freq_i, 'freq_j': freq_j,
                    'sum_freq': sum_freq, 'diff_freq': diff_freq,
                    'kappa_i': kappa_k[i], 'kappa_j': kappa_k[j],
                    'd_kl': d_kl[pair_idx],
                    'J_3rd': J_3rd_single,
                    'J_at_peak': J_at_peak
                })
    
    all_peak_data[peak_name] = {
        'center': peak_center,
        'width': peak_width,
        'linear': peak_linear,
        '3rd': peak_3rd
    }
    
    print(f"\nPeak {peak_name} cm⁻¹ ({peak_center-peak_width:.0f}-{peak_center+peak_width:.0f} cm⁻¹):")
    print(f"  Linear contributions: {len(peak_linear)}")
    print(f"  3rd order contributions: {len(peak_3rd)}")

# Generate decomposition plots for all peaks
for peak_name, peak_data in all_peak_data.items():
    peak_center = peak_data['center']
    peak_width = peak_data['width']
    peak_linear = peak_data['linear']
    peak_3rd = peak_data['3rd']
    
    print(f"\nGenerating peak {peak_name} cm⁻¹ decomposition plot...")
    
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)
    
    # Panel 1: Overall comparison
    ax1 = fig.add_subplot(gs[0, :])
    mask_region = (freq_grid >= peak_center - peak_width) & (freq_grid <= peak_center + peak_width)
    ax1.plot(freq_grid[mask_region], J_total[mask_region], 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid[mask_region], J_total_model[mask_region], 'r--', lw=1.5, label='Full model', alpha=0.8)
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Spectral Density around {peak_center:.0f} cm⁻¹', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Linear contributions (peak-region modes only)
    ax2 = fig.add_subplot(gs[1, 0])
    J_linear_sum = np.zeros_like(freq_grid)
    for i, freq, coeff, J_mode in peak_linear:
        ax2.plot(freq_grid[mask_region], J_mode[mask_region], '-', alpha=0.6, lw=1, 
                label=f'Mode {mode_indices[i]} ({freq:.1f} cm⁻¹)')
        J_linear_sum += J_mode
    ax2.plot(freq_grid[mask_region], J_linear_sum[mask_region], 'b-', lw=2, alpha=0.8, label='Sum')
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax2.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax2.set_title(f'Linear ({len(peak_linear)} modes in region)', fontsize=11, fontweight='bold')
    if len(peak_linear) <= 10:
        ax2.legend(fontsize=8, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: 3rd order contributions (pairs with actual contribution at this peak)
    ax3 = fig.add_subplot(gs[1, 1])
    colors_3rd = plt.cm.tab10(np.linspace(0, 1, min(10, max(1, len(peak_3rd)))))
    
    # Sort by absolute contribution at this peak region (already computed)
    peak_3rd_sorted = sorted(peak_3rd, key=lambda x: abs(x['J_at_peak']), reverse=True)
    
    # Sum of all contributing pairs
    J_3rd_peak_sum = np.zeros_like(freq_grid)
    for item in peak_3rd:
        J_3rd_peak_sum += item['J_3rd']
    
    # Plot top 10 individual contributions
    for idx, item in enumerate(peak_3rd_sorted[:10]):
        fi, fj = item['freq_i'], item['freq_j']
        sign = '+' if item['J_at_peak'] > 0 else '-'
        ax3.plot(freq_grid[mask_region], item['J_3rd'][mask_region], '-', alpha=0.6, lw=1,
                color=colors_3rd[idx % 10],
                label=f'{item["mode_i"]}×{item["mode_j"]} ({fi:.0f},{fj:.0f}) {sign}')
    
    # Show sum of contributing pairs vs full J_3rd_model
    ax3.plot(freq_grid[mask_region], J_3rd_peak_sum[mask_region], 'purple', lw=2, alpha=0.8, 
             label=f'Sum ({len(peak_3rd)} pairs)')
    ax3.plot(freq_grid[mask_region], J_3rd_model[mask_region], 'k--', lw=2, alpha=0.8, 
             label='J_3rd (full model)')
    
    # Calculate coverage statistics at peak
    J_total_at_peak = np.mean(J_total[mask_region])
    J_3rd_peak_at_peak = np.mean(J_3rd_peak_sum[mask_region])
    J_3rd_full_at_peak = np.mean(J_3rd_model[mask_region])
    full_ratio = 100 * J_3rd_full_at_peak / J_total_at_peak if J_total_at_peak > 0 else 0
    
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax3.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax3.set_title(f'3rd Order: J_3rd/J_total = {full_ratio:.1f}%', fontsize=11, fontweight='bold')
    if len(peak_3rd) > 0:
        ax3.legend(fontsize=7, loc='best', ncol=2)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Full model decomposition
    ax4 = fig.add_subplot(gs[2, :])
    ax4.plot(freq_grid[mask_region], J_total[mask_region], 'k-', lw=2, label='J_total (data)', alpha=0.8, zorder=3)
    ax4.fill_between(freq_grid[mask_region], 0, J_linear[mask_region], 
                     alpha=0.5, color='blue', label=f'Linear (full)', zorder=1)
    ax4.fill_between(freq_grid[mask_region], J_linear[mask_region], 
                     J_linear[mask_region] + J_3rd_model[mask_region],
                     alpha=0.5, color='purple', label=f'3rd Order (full)', zorder=2)
    ax4.plot(freq_grid[mask_region], J_total_model[mask_region], 
            'g--', lw=1.5, label='Full model', alpha=0.8, zorder=3)
    
    # Calculate full model coverage at peak
    J_model_at_peak = np.mean(J_total_model[mask_region])
    J_linear_at_peak = np.mean(J_linear[mask_region])
    model_ratio = 100 * J_model_at_peak / J_total_at_peak if J_total_at_peak > 0 else 0
    linear_ratio = 100 * J_linear_at_peak / J_total_at_peak if J_total_at_peak > 0 else 0
    third_ratio = 100 * J_3rd_full_at_peak / J_total_at_peak if J_total_at_peak > 0 else 0
    
    ax4.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax4.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax4.set_title(f'Full Model at {peak_center:.0f} cm⁻¹: Linear {linear_ratio:.1f}% + 3rd {third_ratio:+.1f}% = {model_ratio:.1f}%', 
                  fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.savefig(f'peak_decomposition_3rd_{peak_name}cm.png', dpi=150, bbox_inches='tight')
    print(f"Saved: peak_decomposition_3rd_{peak_name}cm.png")
    plt.close()

# Print detailed contributions for all peaks
for peak_name, peak_data in all_peak_data.items():
    peak_center = peak_data['center']
    peak_linear = peak_data['linear']
    peak_3rd = peak_data['3rd']
    
    print("\n" + "="*70)
    print(f"Detailed contributions to {peak_center:.0f} cm⁻¹ peak:")
    print("="*70)
    
    if peak_linear:
        print(f"\nLinear terms ({len(peak_linear)}):")
        for i, freq, coeff, J_mode in sorted(peak_linear, key=lambda x: x[2], reverse=True)[:15]:
            integral = np.trapezoid(J_mode, freq_grid)
            print(f"  Mode {mode_indices[i]:2d} ({freq:6.1f} cm⁻¹): c_k = {coeff:.6e}, ∫J = {integral:.6e}")
    else:
        print("\n  No direct linear contributions (no modes in frequency range)")
    
    if peak_3rd:
        print(f"\n3rd order terms (top 15 of {len(peak_3rd)}):")
        peak_3rd_sorted = sorted(peak_3rd, key=lambda x: abs(x['d_kl']), reverse=True)
        for item in peak_3rd_sorted[:15]:
            integral = np.trapezoid(item['J_3rd'], freq_grid)
            peak_type = "Sum" if abs(item['sum_freq'] - peak_center) < abs(item['diff_freq'] - peak_center) else "Diff"
            target_f = item['sum_freq'] if peak_type == "Sum" else item['diff_freq']
            print(f"  Mode {item['mode_i']:2d}×{item['mode_j']:2d} ({item['freq_i']:6.1f}×{item['freq_j']:6.1f}): "
                  f"{peak_type}={target_f:6.1f}, κ_i={item['kappa_i']:.2e}, κ_j={item['kappa_j']:.2e}, "
                  f"d_kl={item['d_kl']:.2e}, ∫J={integral:.2e}")
    else:
        print("\n  No active 3rd order contributions to this peak")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("SUMMARY: v4 (Linear + 3rd Order Only)")
print("="*70)
print(f"\nTwo-Stage Fitting:")
print(f"  Stage 1 (Linear NNLS):  R²={R2_linear:.4f}, capture={100*J_linear_integral/J_total_integral:.1f}%")
print(f"  Stage 2 (3rd order):    ΔR²={R2_total-R2_linear:+.4f}, capture={100*J_3rd_integral/J_total_integral:+.1f}%")
print(f"  Final:                  R²={R2_total:.4f}, total={100*J_total_model_integral/J_total_integral:.1f}%")
print(f"\nMode pairs: {N_pairs}")
print(f"  With non-zero κ: {pairs_with_kappa}")
print(f"  Active (|d_kl| > 1e-10): {np.sum(np.abs(d_kl) > 1e-10)}")
print(f"\nNote: 3rd order requires both κ_k (linear coefficient) and d_kl")
print(f"      to be non-zero for significant contribution.")
