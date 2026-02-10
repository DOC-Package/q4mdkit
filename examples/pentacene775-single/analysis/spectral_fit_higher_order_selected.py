#!/usr/bin/env python3
"""
Fit residual from 2nd-order spectral fit using 3rd and 4th order correlations.
Uses user-specified modes from v3 cross term analysis.

Selected modes: 16,17,19,26,32,34,35,36,37,38,39,70,71,79,85
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os.path

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("Higher-Order Spectral Fitting: 3rd + 4th Order Correlations")
print("Using selected modes from v3 cross term analysis")
print("="*70)

# ============================================================================
# Load data
# ============================================================================
print("\nLoading data...")

# Load residual from 2nd-order fit
residual_data = np.loadtxt('spectral_density_quadratic_form.dat')
freq_grid = residual_data[:, 0]
J_total = residual_data[:, 1]
J_model_2nd = residual_data[:, 2]
residual_2nd = residual_data[:, 5]

print(f"  Frequency grid: {len(freq_grid)} points ({freq_grid[0]:.0f}-{freq_grid[-1]:.0f} cm⁻¹)")
print(f"  Residual RMS: {np.sqrt(np.mean(residual_2nd**2)):.6e}")

# Load mode information
with open('g_k_coefficients_below_2000cm.dat') as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith('#') and 'R²' not in l]
    data = [l.split() for l in lines if len(l.split()) >= 3]

all_mode_indices = np.array([int(d[0]) for d in data])
all_mode_freqs = np.array([float(d[1]) for d in data])

# Load κ_k from 2nd-order fit (g_k coefficients)
kappa = {}
with open('g_k_quadratic_form.dat') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            mode_idx = int(parts[0])
            g_k = float(parts[2])
            kappa[mode_idx] = g_k

print(f"  Loaded κ_k for {len(kappa)} modes")

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords.txt')
q_k_all = q_k_data[:, 1:]

N_frames = len(q_k_all)
print(f"  Frames: {N_frames}")

# ============================================================================
# User-specified mode selection
# ============================================================================
print("\n" + "="*70)
print("User-specified mode selection")
print("="*70)

# Selected modes from v3 cross term analysis
selected_modes = [16, 17, 19, 26, 32, 34, 35, 36, 37, 38, 39, 70, 71, 79, 85]
print(f"\nSelected modes: {selected_modes}")

# Find local indices for these modes
selected_local_indices = []
mode_indices = []
mode_freqs = []

for mode_idx in selected_modes:
    # Find in all_mode_indices
    matches = np.where(all_mode_indices == mode_idx)[0]
    if len(matches) > 0:
        local_idx = matches[0]
        selected_local_indices.append(local_idx)
        mode_indices.append(mode_idx)
        mode_freqs.append(all_mode_freqs[local_idx])
        print(f"  Mode {mode_idx}: freq = {all_mode_freqs[local_idx]:.1f} cm⁻¹")
    else:
        print(f"  Mode {mode_idx}: NOT FOUND in data")

mode_indices = np.array(mode_indices)
mode_freqs = np.array(mode_freqs)
N_selected = len(selected_local_indices)

print(f"\nTotal: {N_selected} modes found")
print(f"  Frequency range: {min(mode_freqs):.1f} - {max(mode_freqs):.1f} cm⁻¹")

# Load mode coordinates for selected modes
mode_coords = q_k_all[:, mode_indices - 7]  # Adjust indices (mode 7 is column 0)

# Target frequencies
target_freqs = [524, 800, 1092, 1231]  # cm⁻¹
print(f"\nTarget residual frequencies: {target_freqs} cm⁻¹")

# Build ALL pairs from selected modes
mode_pairs = []
pair_info = {}

for idx_i in range(N_selected):
    for idx_j in range(idx_i + 1, N_selected):
        freq_i = mode_freqs[idx_i]
        freq_j = mode_freqs[idx_j]
        mode_i = mode_indices[idx_i]
        mode_j = mode_indices[idx_j]
        
        # Check which target frequency this pair might contribute to
        freq_sum = freq_i + freq_j
        freq_diff = abs(freq_i - freq_j)
        target_matched = 0
        match_type = 'none'
        
        tolerance = 50
        for target in target_freqs:
            if abs(freq_sum - target) <= tolerance:
                target_matched = target
                match_type = 'sum'
                break
            elif abs(freq_diff - target) <= tolerance:
                target_matched = target
                match_type = 'diff'
                break
        
        mode_pairs.append((idx_i, idx_j))
        pair_info[(idx_i, idx_j)] = {
            'mode_i': mode_i, 'mode_j': mode_j,
            'freq_i': freq_i, 'freq_j': freq_j,
            'target': target_matched, 'type': match_type,
            'result': freq_diff if match_type == 'diff' else freq_sum
        }

N_pairs = len(mode_pairs)
print(f"\nGenerated {N_pairs} mode pairs (all pairs from {N_selected} selected modes)")
print(f"  Expected: {N_selected}×({N_selected}-1)/2 = {N_selected*(N_selected-1)//2}")

# Count pairs matching each target frequency
print(f"\nTarget frequency coverage (pairs within ±50 cm⁻¹):")
for target in target_freqs:
    count = sum(1 for p, info in pair_info.items() if info['target'] == target)
    print(f"  {target} cm⁻¹: {count} pairs")

# Show pairs matching targets
print(f"\nPairs matching target frequencies:")
for target in target_freqs:
    matching = [(p, info) for p, info in pair_info.items() if info['target'] == target]
    if matching:
        print(f"\n  {target} cm⁻¹:")
        for (i, j), info in matching[:5]:
            print(f"    Mode {info['mode_i']:2d} ({info['freq_i']:.1f}) × Mode {info['mode_j']:2d} ({info['freq_j']:.1f}) "
                  f"-> {info['type']} = {info['result']:.1f}")
        if len(matching) > 5:
            print(f"    ... and {len(matching)-5} more")

# ============================================================================
# Compute 3rd and 4th order correlation functions
# ============================================================================
print("\n" + "="*70)
print("Computing higher-order correlations")
print("="*70)

cache_file = 'spectral_cache_higher_order_selected.npz'

def compute_spectral_from_corr(corr, freq):
    """Compute spectral density from correlation via Fourier transform"""
    t = np.arange(len(corr)) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = corr * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    # Multiply by prefactor
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

cache_valid = False
if os.path.exists(cache_file):
    print(f"\nChecking cached correlations from {cache_file}...")
    cache_data = np.load(cache_file, allow_pickle=True)
    
    # Check if cache has all required data
    has_all_types = all(key in cache_data for key in ['J_3rd_type1', 'J_3rd_type2', 'J_3rd_type3', 'J_3rd_type4', 'J_4th'])
    
    if has_all_types:
        cached_pairs = [tuple(p) for p in cache_data['pairs']]
        cached_modes = list(cache_data['selected_modes']) if 'selected_modes' in cache_data else []
        
        # Check if modes and pairs match
        if cached_modes == selected_modes and cached_pairs == mode_pairs:
            J_3rd_type1 = cache_data['J_3rd_type1']
            J_3rd_type2 = cache_data['J_3rd_type2']
            J_3rd_type3 = cache_data['J_3rd_type3']
            J_3rd_type4 = cache_data['J_3rd_type4']
            J_4th = cache_data['J_4th']
            print(f"  Loaded {len(cached_pairs)} pairs with all 4 types of 3rd order")
            cache_valid = True
        else:
            print(f"  Cache mismatch (modes or pairs differ), recomputing...")
    else:
        print(f"  Cache missing required data, recomputing...")

if not cache_valid:
    print(f"\nComputing 3rd and 4th order correlations for {N_pairs} pairs...")
    print(f"  3rd order: 4 types per pair")
    print(f"    Type 1: <q_k(t)q_l(t)q_k(0)> with coeff κ_k d_kl")
    print(f"    Type 2: <q_k(t)q_l(t)q_l(0)> with coeff κ_l d_kl")
    print(f"    Type 3: <q_k(t)q_k(0)q_l(0)> with coeff κ_l d_kl")
    print(f"    Type 4: <q_l(t)q_l(0)q_k(0)> with coeff κ_k d_kl")
    print(f"  4th order: <q_k(t)q_l(t)q_k(0)q_l(0)> with coeff d_kl²")
    
    J_3rd_type1 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type2 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type3 = np.zeros((N_pairs, len(freq_grid)))
    J_3rd_type4 = np.zeros((N_pairs, len(freq_grid)))
    J_4th = np.zeros((N_pairs, len(freq_grid)))
    
    for pair_idx, (i, j) in enumerate(mode_pairs):
        q_k = mode_coords[:, i]
        q_l = mode_coords[:, j]
        
        # Center the coordinates
        q_k_c = q_k - np.mean(q_k)
        q_l_c = q_l - np.mean(q_l)
        
        # Products
        product_kl = q_k_c * q_l_c
        
        # Window
        t_arr = np.arange(max_lag)
        window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag))
        
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
        
        # 4th order: <q_k(t)q_l(t)q_k(0)q_l(0)>
        corr_4th = np.zeros(max_lag)
        for lag in range(max_lag):
            corr_4th[lag] = np.mean(product_kl[lag:] * product_kl[:N_frames-lag])
        corr_4th *= window
        J_4th[pair_idx, :] = compute_spectral_from_corr(corr_4th, freq_grid)
        
        if (pair_idx + 1) % 20 == 0 or pair_idx == 0:
            print(f"    Progress: {pair_idx+1}/{N_pairs} pairs")
    
    # Save to cache
    np.savez_compressed(cache_file, 
                        J_3rd_type1=J_3rd_type1, 
                        J_3rd_type2=J_3rd_type2,
                        J_3rd_type3=J_3rd_type3,
                        J_3rd_type4=J_3rd_type4,
                        J_4th=J_4th,
                        pairs=np.array(mode_pairs),
                        selected_modes=np.array(selected_modes))
    print(f"\n  Saved to {cache_file}")
    print(f"  File size: {os.path.getsize(cache_file) / 1024**2:.1f} MB")

print(f"\nCorrelation statistics:")
print(f"  J_3rd_type1 RMS: {np.sqrt(np.mean(J_3rd_type1**2)):.6e}")
print(f"  J_3rd_type2 RMS: {np.sqrt(np.mean(J_3rd_type2**2)):.6e}")
print(f"  J_3rd_type3 RMS: {np.sqrt(np.mean(J_3rd_type3**2)):.6e}")
print(f"  J_3rd_type4 RMS: {np.sqrt(np.mean(J_3rd_type4**2)):.6e}")
print(f"  J_4th RMS: {np.sqrt(np.mean(J_4th**2)):.6e}")

# ============================================================================
# Setup fitting
# ============================================================================
print("\n" + "="*70)
print("Setting up optimization")
print("="*70)

# Get κ values for each pair
kappa_k = np.array([kappa.get(mode_indices[i], 0.0) for i, j in mode_pairs])
kappa_l = np.array([kappa.get(mode_indices[j], 0.0) for i, j in mode_pairs])

print(f"\nκ statistics:")
print(f"  Mean |κ_k|: {np.mean(np.abs(kappa_k)):.6e}")
print(f"  Mean |κ_l|: {np.mean(np.abs(kappa_l)):.6e}")

# Normalize for numerical stability
residual_rms = np.sqrt(np.mean(residual_2nd**2))
J_3rd_type1_norm = J_3rd_type1 / residual_rms
J_3rd_type2_norm = J_3rd_type2 / residual_rms
J_3rd_type3_norm = J_3rd_type3 / residual_rms
J_3rd_type4_norm = J_3rd_type4 / residual_rms
J_4th_norm = J_4th / residual_rms
residual_norm = residual_2nd / residual_rms

print(f"\nNormalization:")
print(f"  Scale factor: {residual_rms:.6e}")
print(f"  Normalized residual RMS: {np.sqrt(np.mean(residual_norm**2)):.4f}")

# ============================================================================
# Optimization
# ============================================================================
print("\n" + "="*70)
print("Optimizing d_kl to fit residual with 3rd + 4th order terms")
print("="*70)

# Precompute combined 3rd order basis (κ_k * type1 + κ_l * type2 + κ_l * type3 + κ_k * type4)
J_3rd_combined_norm = (kappa_k[:, None] * J_3rd_type1_norm + 
                       kappa_l[:, None] * J_3rd_type2_norm +
                       kappa_l[:, None] * J_3rd_type3_norm + 
                       kappa_k[:, None] * J_3rd_type4_norm)

def compute_model(d_kl, include_4th=True):
    """Compute model spectral density from d_kl coefficients."""
    # 3rd order contribution (all 4 types combined)
    J_3rd = np.einsum('p,pw->w', d_kl, J_3rd_combined_norm)
    
    if include_4th:
        J_4th_contrib = np.einsum('p,pw->w', d_kl**2, J_4th_norm)
    else:
        J_4th_contrib = np.zeros_like(J_3rd)
    
    return J_3rd + J_4th_contrib

def objective(d_kl, include_4th=True):
    """Objective function: 0.5 * ||residual - model||²"""
    J_model = compute_model(d_kl, include_4th)
    r = residual_norm - J_model
    return 0.5 * np.sum(r**2)

def gradient(d_kl, include_4th=True):
    """Gradient of objective function."""
    J_model = compute_model(d_kl, include_4th)
    r = residual_norm - J_model  # shape: (N_freq,)
    
    # d(objective)/d(d_kl[p]) = -sum_w r[w] * d(J_model)/d(d_kl[p])
    # d(J_model)/d(d_kl[p]) = J_3rd_combined[p,w] + 2*d_kl[p]*J_4th[p,w]
    grad = -np.einsum('w,pw->p', r, J_3rd_combined_norm)
    if include_4th:
        grad -= 2 * np.einsum('w,pw,p->p', r, J_4th_norm, d_kl)
    return grad

# Initial guess
np.random.seed(42)
d_kl_init = np.random.randn(N_pairs) * 0.01
bound_limit = 1000.0
bounds = [(-bound_limit, bound_limit) for _ in range(N_pairs)]

print(f"\nOptimization settings:")
print(f"  Parameters: {N_pairs} (d_kl)")
print(f"  Data points: {len(freq_grid)}")
print(f"  Method: L-BFGS-B")
print(f"  Initial ||d_kl||: {np.linalg.norm(d_kl_init):.6e}")

# Step 1: Get initial guess from 3rd order only (linear problem)
print("\n  Step 1: Linear fit using 3rd order only...")
result_3rd = minimize(
    lambda x: objective(x, include_4th=False),
    d_kl_init,
    method='L-BFGS-B',
    jac=lambda x: gradient(x, include_4th=False),
    bounds=bounds,
    options={'maxiter': 1000, 'ftol': 1e-10, 'gtol': 1e-8}
)
d_kl_init_refined = result_3rd.x
R2_3rd = 1 - 2 * result_3rd.fun / np.sum(residual_norm**2)
print(f"    3rd-order only R²: {R2_3rd:.4f}")
print(f"    Max |d_kl| from 3rd-only: {np.max(np.abs(d_kl_init_refined)):.6e}")
print(f"    nfev: {result_3rd.nfev}")

# Step 2: Full nonlinear optimization with 4th order
print("\n  Step 2: Nonlinear fit with 3rd + 4th order...")
result = minimize(
    lambda x: objective(x, include_4th=True),
    d_kl_init_refined,
    method='L-BFGS-B',
    jac=lambda x: gradient(x, include_4th=True),
    bounds=bounds,
    options={'maxiter': 5000, 'ftol': 1e-10, 'gtol': 1e-8, 'disp': True}
)

d_kl_opt = result.x

print(f"\nOptimization result:")
print(f"  Success: {result.success}")
print(f"  Message: {result.message}")
print(f"  ||d_kl||: {np.linalg.norm(d_kl_opt):.6e}")
print(f"  Mean |d_kl|: {np.mean(np.abs(d_kl_opt)):.6e}")
print(f"  Max |d_kl|: {np.max(np.abs(d_kl_opt)):.6e}")

# ============================================================================
# Compute final model and fit quality
# ============================================================================
print("\n" + "="*70)
print("Fit Quality")
print("="*70)

# Compute contributions (in original units)
J_3rd_contrib = np.einsum('p,pw->w', kappa_k * d_kl_opt, J_3rd_type1)
J_3rd_contrib += np.einsum('p,pw->w', kappa_l * d_kl_opt, J_3rd_type2)
J_3rd_contrib += np.einsum('p,pw->w', kappa_l * d_kl_opt, J_3rd_type3)
J_3rd_contrib += np.einsum('p,pw->w', kappa_k * d_kl_opt, J_3rd_type4)
J_4th_contrib = np.einsum('p,pw->w', d_kl_opt**2, J_4th)
J_higher_model = J_3rd_contrib + J_4th_contrib

# New residual
new_residual = residual_2nd - J_higher_model

# Total model including 2nd order
J_model_total = J_model_2nd + J_higher_model

# Fit quality for residual
ss_res = np.sum(new_residual**2)
ss_tot = np.sum(residual_2nd**2)
R2_residual = 1 - ss_res / ss_tot

# Fit quality for total
ss_res_total = np.sum((J_total - J_model_total)**2)
ss_tot_total = np.sum((J_total - np.mean(J_total))**2)
R2_total = 1 - ss_res_total / ss_tot_total

# Integrals
J_total_integral = np.trapezoid(J_total, freq_grid)
J_model_total_integral = np.trapezoid(J_model_total, freq_grid)
J_3rd_integral = np.trapezoid(J_3rd_contrib, freq_grid)
J_4th_integral = np.trapezoid(J_4th_contrib, freq_grid)

print(f"\n  Residual fitting:")
print(f"    R² (residual fit): {R2_residual:.4f}")
print(f"    Original residual RMS: {np.sqrt(np.mean(residual_2nd**2)):.6e}")
print(f"    New residual RMS: {np.sqrt(np.mean(new_residual**2)):.6e}")
print(f"    Reduction: {100*(1 - np.sqrt(np.mean(new_residual**2))/np.sqrt(np.mean(residual_2nd**2))):.1f}%")

print(f"\n  Total fitting (2nd + 3rd + 4th order):")
print(f"    R² (total): {R2_total:.4f}")
print(f"    ∫J_total = {J_total_integral:.6e} cm⁻¹²")
print(f"    ∫J_model_total = {J_model_total_integral:.6e} cm⁻¹²")
print(f"    Capture: {100*J_model_total_integral/J_total_integral:.1f}%")

print(f"\n  Contribution breakdown:")
print(f"    2nd order: {100*np.trapezoid(J_model_2nd, freq_grid)/J_total_integral:.1f}%")
print(f"    3rd order: {100*J_3rd_integral/J_total_integral:.1f}%")
print(f"    4th order: {100*J_4th_integral/J_total_integral:.1f}%")

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
    J_model_at_target = np.mean(J_model_total[freq_window])
    J_residual_at_target = np.mean(new_residual[freq_window])
    capture_pct = 100 * J_model_at_target / J_total_at_target if J_total_at_target > 0 else 0
    
    print(f"{target:8.0f} {J_total_at_target:12.3e} {J_model_at_target:12.3e} "
          f"{J_residual_at_target:12.3e} {capture_pct:9.1f}%")

# ============================================================================
# Top mode pairs
# ============================================================================
print("\n" + "="*70)
print(f"All {N_pairs} mode pairs by |d_kl|")
print("="*70)

top_idx = np.argsort(np.abs(d_kl_opt))[::-1]
print(f"\n{'Rank':>4} {'Mode k':>7} {'Mode l':>7} {'ω_k':>8} {'ω_l':>8} {'d_kl':>12} {'Target':>8} {'Type':>6}")
print("-"*75)
for rank, idx in enumerate(top_idx, 1):
    i, j = mode_pairs[idx]
    info = pair_info[(i, j)]
    sign = '+' if d_kl_opt[idx] > 0 else '-'
    print(f"{rank:4d} {info['mode_i']:7d} {info['mode_j']:7d} {info['freq_i']:8.1f} {info['freq_j']:8.1f} "
          f"{sign}{abs(d_kl_opt[idx]):11.6e} {info['target']:8d} {info['type']:>6}")

# ============================================================================
# Visualization
# ============================================================================
print("\nGenerating visualization...")

fig, axes = plt.subplots(3, 2, figsize=(14, 12))

# Panel 1: Original vs Total model
ax1 = axes[0, 0]
ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
ax1.plot(freq_grid, J_model_total, 'r-', lw=1.5, label=f'J_model (R²={R2_total:.4f})', alpha=0.8)
ax1.set_ylabel('J(ν) (cm⁻¹²)')
ax1.set_title(f'Total Fit: 2nd + 3rd + 4th Order\nModes: {selected_modes}')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Residual comparison
ax2 = axes[0, 1]
ax2.plot(freq_grid, residual_2nd, 'b-', lw=1.5, label='2nd order residual', alpha=0.7)
ax2.plot(freq_grid, new_residual, 'g-', lw=1.5, label='After 3rd+4th fit', alpha=0.7)
ax2.axhline(0, color='k', ls='--', lw=0.5)
for target in target_freqs:
    ax2.axvline(target, color='r', ls=':', alpha=0.5)
ax2.set_ylabel('Residual (cm⁻¹²)')
ax2.set_title(f'Residual Comparison (R²={R2_residual:.4f})')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Panel 3: 3rd order contribution
ax3 = axes[1, 0]
ax3.plot(freq_grid, J_3rd_contrib, 'm-', lw=1.5, alpha=0.8)
ax3.axhline(0, color='k', ls='--', lw=0.5)
for target in target_freqs:
    ax3.axvline(target, color='r', ls=':', alpha=0.5)
ax3.set_ylabel('J(ν) (cm⁻¹²)')
ax3.set_title(f'3rd Order Contribution ({100*J_3rd_integral/J_total_integral:.1f}%)')
ax3.grid(True, alpha=0.3)

# Panel 4: 4th order contribution
ax4 = axes[1, 1]
ax4.plot(freq_grid, J_4th_contrib, 'c-', lw=1.5, alpha=0.8)
ax4.axhline(0, color='k', ls='--', lw=0.5)
for target in target_freqs:
    ax4.axvline(target, color='r', ls=':', alpha=0.5)
ax4.set_ylabel('J(ν) (cm⁻¹²)')
ax4.set_title(f'4th Order Contribution ({100*J_4th_integral/J_total_integral:.1f}%)')
ax4.grid(True, alpha=0.3)

# Panel 5: Stacked contributions
ax5 = axes[2, 0]
ax5.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
ax5.fill_between(freq_grid, 0, J_model_2nd, alpha=0.4, color='blue', label='2nd order')
ax5.fill_between(freq_grid, J_model_2nd, J_model_2nd + J_3rd_contrib, 
                 alpha=0.4, color='magenta', label='3rd order')
ax5.fill_between(freq_grid, J_model_2nd + J_3rd_contrib, J_model_total,
                 alpha=0.4, color='cyan', label='4th order')
ax5.set_xlabel('Frequency (cm⁻¹)')
ax5.set_ylabel('J(ν) (cm⁻¹²)')
ax5.set_title('Stacked Contributions')
ax5.legend()
ax5.grid(True, alpha=0.3)

# Panel 6: d_kl bar chart
ax6 = axes[2, 1]
x_labels = [f"{pair_info[p]['mode_i']}×{pair_info[p]['mode_j']}" for p in mode_pairs]
colors = ['green' if d > 0 else 'red' for d in d_kl_opt]
bars = ax6.bar(range(N_pairs), d_kl_opt, color=colors, alpha=0.7)
ax6.axhline(0, color='k', ls='-', lw=0.5)
ax6.set_xticks(range(N_pairs))
ax6.set_xticklabels(x_labels, rotation=90, fontsize=7)
ax6.set_xlabel('Mode pair')
ax6.set_ylabel('d_kl')
ax6.set_title(f'd_kl Coefficients ({N_pairs} pairs)')
ax6.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('spectral_fit_higher_order_selected.png', dpi=150, bbox_inches='tight')
print(f"Saved: spectral_fit_higher_order_selected.png")
plt.close()

# ============================================================================
# Save results
# ============================================================================
print("\nSaving results...")

# Save d_kl coefficients
with open('d_kl_higher_order_selected.dat', 'w') as f:
    f.write("# d_kl coefficients from 3rd+4th order fit (selected modes)\n")
    f.write(f"# Selected modes: {selected_modes}\n")
    f.write(f"# Total pairs: {N_pairs}\n")
    f.write(f"# R² (residual): {R2_residual:.6f}\n")
    f.write(f"# R² (total): {R2_total:.6f}\n")
    f.write("#\n")
    f.write("# mode_k  mode_l  freq_k  freq_l  d_kl  target_freq  type\n")
    for idx, (i, j) in enumerate(mode_pairs):
        info = pair_info[(i, j)]
        f.write(f"{info['mode_i']:4d}  {info['mode_j']:4d}  {info['freq_i']:8.2f}  {info['freq_j']:8.2f}  "
                f"{d_kl_opt[idx]:14.6e}  {info['target']:6d}  {info['type']}\n")
print(f"Saved: d_kl_higher_order_selected.dat")

# Save spectral densities
with open('spectral_density_higher_order_selected.dat', 'w') as f:
    f.write("# Frequency(cm⁻¹)  J_total  J_model_2nd  J_3rd  J_4th  J_model_total  Residual_2nd  Residual_final\n")
    for i in range(len(freq_grid)):
        f.write(f"{freq_grid[i]:10.4f}  {J_total[i]:14.6e}  {J_model_2nd[i]:14.6e}  "
                f"{J_3rd_contrib[i]:14.6e}  {J_4th_contrib[i]:14.6e}  {J_model_total[i]:14.6e}  "
                f"{residual_2nd[i]:14.6e}  {new_residual[i]:14.6e}\n")
print(f"Saved: spectral_density_higher_order_selected.dat")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nHigher-Order Spectral Fitting (Selected Modes):")
print(f"  Selected modes: {selected_modes}")
print(f"  Mode pairs used: {N_pairs}")
print(f"  Target frequencies: {target_freqs} cm⁻¹")
print(f"\nFit Quality:")
print(f"  2nd order R²: {1 - np.sum(residual_2nd**2)/ss_tot_total:.4f}")
print(f"  Total R² (2nd+3rd+4th): {R2_total:.4f}")
print(f"  Improvement: {R2_total - (1 - np.sum(residual_2nd**2)/ss_tot_total):.4f}")
print(f"\nContributions to ∫J:")
print(f"  2nd order: {100*np.trapezoid(J_model_2nd, freq_grid)/J_total_integral:.1f}%")
print(f"  3rd order: {100*J_3rd_integral/J_total_integral:.1f}%")
print(f"  4th order: {100*J_4th_integral/J_total_integral:.1f}%")
print(f"  Total: {100*J_model_total_integral/J_total_integral:.1f}%")
print(f"\nResidual reduction:")
print(f"  Before (2nd only): RMS = {np.sqrt(np.mean(residual_2nd**2)):.6e}")
print(f"  After (2nd+3rd+4th): RMS = {np.sqrt(np.mean(new_residual**2)):.6e}")
print(f"  Reduction: {100*(1 - np.sqrt(np.mean(new_residual**2))/np.sqrt(np.mean(residual_2nd**2))):.1f}%")
