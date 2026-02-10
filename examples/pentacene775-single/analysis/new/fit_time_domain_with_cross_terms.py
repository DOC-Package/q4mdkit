#!/usr/bin/env python3
"""
Time-domain fitting with cross terms:
  ΔE(t) = Σ_k κ_k Δq_k(t) + Σ_{k≤l} d_{kl} Δq_k(t)Δq_l(t)

This is analogous to spectral_fit_with_cross_terms_v3.py but in time domain.

Stage 1: Fit linear terms using NNLS (κ_k ≥ 0)
Stage 2: Fit quadratic terms to residual using NNLS (d_{kl} ≥ 0)
"""
import numpy as np
from scipy.optimize import nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs

print("="*70)
print("Time-Domain Fitting with Cross Terms")
print("="*70)

# Parameters
n_frames_fit = None  # Use all frames (or set to e.g. 5000 for subset)

# ============================================================================
# Load data
# ============================================================================
print("\nLoading data...")

# Load energy difference
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff_full = energy_diff_data[:, 4]  # a.u.

# Load mode information from stats file
stats_data = []
with open('mode_coords_new_stats.txt') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            stats_data.append([int(parts[0]), float(parts[1])])  # Mode, Freq

stats_data = np.array(stats_data)
all_mode_indices = stats_data[:, 0].astype(int)
all_mode_freqs = stats_data[:, 1]

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords_new.txt')
q_k_all = q_k_data[:, 1:]

N_frames_full = min(len(energy_diff_full), len(q_k_all))
N_modes_total = len(all_mode_freqs)

# Map mode indices to column indices
col_indices_all = all_mode_indices - all_mode_indices[0]
mode_coords_all = q_k_all[:N_frames_full, :]

print(f"  Total frames: {N_frames_full}")
print(f"  Total modes: {N_modes_total}")
print(f"  Frequency range: {all_mode_freqs.min():.1f} - {all_mode_freqs.max():.1f} cm⁻¹")

# ============================================================================
# Filter modes for Stage 1: 450 ≤ freq ≤ 1800 cm⁻¹
# ============================================================================
stage1_mask = (all_mode_freqs >= 450) & (all_mode_freqs <= 1800)
print(f"\nStage 1 modes (450-1800 cm⁻¹): {np.sum(stage1_mask)}")

# ============================================================================
# Select subset for fitting
# ============================================================================
if n_frames_fit is None:
    n_frames_fit = N_frames_full
else:
    n_frames_fit = min(n_frames_fit, N_frames_full)
print(f"\nUsing {n_frames_fit} frames for fitting ({n_frames_fit * dt / 1000:.1f} ps)")

energy_diff = energy_diff_full[:n_frames_fit]
mode_coords = mode_coords_all[:n_frames_fit, :]

# ============================================================================
# Subtract means
# ============================================================================
print("\nSubtracting means...")

energy_mean = np.mean(energy_diff)
energy_diff_centered = energy_diff - energy_mean

mode_means = np.mean(mode_coords, axis=0)
mode_coords_centered = mode_coords - mode_means

print(f"  Energy mean: {energy_mean:.6e} a.u.")
print(f"  Energy std: {np.std(energy_diff_centered):.6e} a.u.")

# ============================================================================
# Stage 1: Linear terms with NNLS (450-1800 cm⁻¹ modes)
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, κ_k ≥ 0)")
print(f"Using {np.sum(stage1_mask)} modes with 450 ≤ freq ≤ 1800 cm⁻¹")
print("="*70)

# Design matrix for linear terms (filtered modes)
X_linear_stage1 = mode_coords_centered[:, stage1_mask]
y = energy_diff_centered

# NNLS fit for non-negative constraint
print("\nFitting linear terms with NNLS...")
kappa_stage1, residual_norm = nnls(X_linear_stage1, y)

# Expand to full mode array
kappa_linear = np.zeros(N_modes_total)
kappa_linear[stage1_mask] = kappa_stage1

# Predictions
y_pred_linear = mode_coords_centered @ kappa_linear

# Metrics
residuals_stage1 = y - y_pred_linear
ss_tot = np.sum(y**2)  # Already mean-subtracted
ss_res_stage1 = np.sum(residuals_stage1**2)
R2_linear = 1 - ss_res_stage1 / ss_tot

n_active_linear = np.sum(kappa_linear > 1e-15)

print(f"\nStage 1 Results:")
print(f"  R² = {R2_linear:.4f}")
print(f"  RMSE: {np.sqrt(np.mean(residuals_stage1**2)):.6e} a.u.")
print(f"  Active modes: {n_active_linear}/{np.sum(stage1_mask)}")

# ============================================================================
# Stage 2: Cross terms on residual
# ============================================================================
print("\n" + "="*70)
print("Stage 2: Cross Terms on Residual (NNLS, d_{kl} ≥ 0)")
print("="*70)

# Select modes for cross terms (0-800 cm⁻¹)
freq_min, freq_max = 0, 800
mask_cross = (all_mode_freqs >= freq_min) & (all_mode_freqs <= freq_max)
idx_cross = np.where(mask_cross)[0]
n_cross_modes = len(idx_cross)

print(f"\nModes in {freq_min}-{freq_max} cm⁻¹: {n_cross_modes}")
n_cross_terms = n_cross_modes * (n_cross_modes + 1) // 2
print(f"Number of quadratic terms: {n_cross_terms} (diagonal: {n_cross_modes}, off-diagonal: {n_cross_modes * (n_cross_modes - 1) // 2})")

# Target frequencies for additional pairs
target_freqs = [524, 800, 1092, 1231]  # cm⁻¹
tolerance = 50.0  # cm⁻¹

# ============================================================================
# Cache for cross term design matrix
# ============================================================================
cache_file_cross = 'time_domain_cross_cache.npz'

if os.path.exists(cache_file_cross):
    print(f"\nLoading cached cross term design matrix from {cache_file_cross}...")
    cache_data = np.load(cache_file_cross, allow_pickle=True)
    
    # Verify compatibility
    cached_n_frames = cache_data['n_frames']
    cached_idx_cross = cache_data['idx_cross']
    cached_target_freqs = cache_data['target_freqs']
    
    if (cached_n_frames == n_frames_fit and 
        np.array_equal(cached_idx_cross, idx_cross) and 
        np.array_equal(cached_target_freqs, target_freqs)):
        X_cross = cache_data['X_cross']
        cross_pairs = [tuple(p) for p in cache_data['cross_pairs']]
        n_cross = len(cross_pairs)
        print(f"  Loaded design matrix with {n_cross} cross terms")
        cache_loaded = True
    else:
        print(f"  Cache parameters mismatch, recomputing...")
        os.remove(cache_file_cross)
        cache_loaded = False
else:
    cache_loaded = False

if not cache_loaded and n_cross_modes >= 2:
    print("\nBuilding cross term design matrix...")
    cross_pairs = []
    X_cross_list = []
    
    for i_idx, i in enumerate(idx_cross):
        q_i = mode_coords_centered[:, i]
        
        # Diagonal term: q_i^2
        q_ii = q_i ** 2
        q_ii_centered = q_ii - np.mean(q_ii)  # Center the squared term
        X_cross_list.append(q_ii_centered)
        cross_pairs.append((i, i))
        
        # Off-diagonal terms: q_i * q_j (i < j)
        for j in idx_cross[i_idx+1:]:
            q_j = mode_coords_centered[:, j]
            q_ij = q_i * q_j
            q_ij_centered = q_ij - np.mean(q_ij)  # Center the cross term
            X_cross_list.append(q_ij_centered)
            cross_pairs.append((i, j))
        
        if (i_idx + 1) % 5 == 0:
            print(f"  Progress: {i_idx+1}/{n_cross_modes} modes")
    
    n_cross = len(cross_pairs)
    print(f"\nBase cross terms: {n_cross}")
    
    # ========================================================================
    # Add high-frequency mode combinations for target peaks
    # ========================================================================
    print("\n" + "-"*50)
    print("Adding mode pairs for target residual peaks")
    print("-"*50)
    
    print(f"\nTarget frequencies: {target_freqs} cm⁻¹")
    print(f"Tolerance: ±{tolerance:.0f} cm⁻¹")
    
    # Find mode pairs that generate target frequencies via sum or difference
    additional_pairs_by_target = {t: [] for t in target_freqs}
    
    for i in range(N_modes_total):
        for j in range(i+1, N_modes_total):
            freq_i = all_mode_freqs[i]
            freq_j = all_mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            
            for target in target_freqs:
                if abs(sum_freq - target) < tolerance:
                    if (i, j) not in cross_pairs:
                        additional_pairs_by_target[target].append((i, j, sum_freq, 'sum'))
                
                if abs(diff_freq - target) < tolerance:
                    if (i, j) not in cross_pairs:
                        additional_pairs_by_target[target].append((i, j, diff_freq, 'diff'))
    
    # Report and collect unique pairs
    additional_pairs_set = set()
    print("\nMode pairs matching target frequencies:")
    for target in target_freqs:
        pairs = additional_pairs_by_target[target]
        print(f"\n  Target {target} cm⁻¹: {len(pairs)} pairs found")
        
        pairs_sorted = sorted(pairs, key=lambda x: abs(x[2] - target))
        
        if len(pairs_sorted) > 0:
            print(f"    Top 5 pairs:")
            for rank, (i, j, freq, ptype) in enumerate(pairs_sorted[:5], 1):
                freq_i = all_mode_freqs[i]
                freq_j = all_mode_freqs[j]
                if ptype == 'sum':
                    print(f"      {rank}. Mode {all_mode_indices[i]:2d} ({freq_i:6.1f}) + Mode {all_mode_indices[j]:2d} ({freq_j:6.1f}) = {freq:.1f} cm⁻¹")
                else:
                    print(f"      {rank}. Mode {all_mode_indices[i]:2d} ({freq_i:6.1f}) - Mode {all_mode_indices[j]:2d} ({freq_j:6.1f}) = {freq:.1f} cm⁻¹")
        
        for i, j, freq, ptype in pairs:
            additional_pairs_set.add((i, j))
    
    additional_pairs = [(i, j) for (i, j) in additional_pairs_set if (i, j) not in cross_pairs]
    print(f"\nTotal unique additional pairs: {len(additional_pairs)}")
    
    if len(additional_pairs) > 0:
        print(f"\nComputing cross terms for additional pairs...")
        for idx, (i, j) in enumerate(additional_pairs):
            q_i = mode_coords_centered[:, i]
            q_j = mode_coords_centered[:, j]
            q_ij = q_i * q_j
            q_ij_centered = q_ij - np.mean(q_ij)
            X_cross_list.append(q_ij_centered)
            cross_pairs.append((i, j))
            
            if (idx + 1) % 50 == 0:
                print(f"  Progress: {idx+1}/{len(additional_pairs)} pairs")
        
        n_cross = len(cross_pairs)
        print(f"\nTotal cross terms (including target frequency pairs): {n_cross}")
    
    X_cross = np.column_stack(X_cross_list)
    
    # Save cache
    print(f"\nSaving cross term cache to {cache_file_cross}...")
    np.savez_compressed(cache_file_cross,
                        X_cross=X_cross,
                        cross_pairs=np.array(cross_pairs),
                        idx_cross=idx_cross,
                        target_freqs=np.array(target_freqs),
                        n_frames=n_frames_fit)

# ============================================================================
# Fit cross terms with NNLS
# ============================================================================
if n_cross_modes >= 2:
    print("\n" + "-"*50)
    print("Fitting Cross Terms to Residual (NNLS)")
    print("-"*50)
    
    # Target: residual from Stage 1
    y_residual = residuals_stage1
    
    # NNLS fit
    d_ij, residual_norm = nnls(X_cross, y_residual)
    
    # Predictions
    y_pred_cross = X_cross @ d_ij
    y_pred_total = y_pred_linear + y_pred_cross
    
    # Final residual
    residuals_stage2 = y - y_pred_total
    ss_res_stage2 = np.sum(residuals_stage2**2)
    R2_total = 1 - ss_res_stage2 / ss_tot
    
    n_active_cross = np.sum(d_ij > 1e-15)
    
    print(f"\nStage 2 Results:")
    print(f"  R² (total) = {R2_total:.4f} (linear: {R2_linear:.4f}, improvement: +{R2_total - R2_linear:.4f})")
    print(f"  RMSE: {np.sqrt(np.mean(residuals_stage2**2)):.6e} a.u.")
    print(f"  Active cross terms: {n_active_cross}/{n_cross}")
    
    # Contribution analysis
    var_linear = np.var(y_pred_linear)
    var_cross = np.var(y_pred_cross)
    var_total = np.var(y)
    
    print(f"\n  Variance contributions:")
    print(f"    Linear: {100 * var_linear / var_total:.1f}%")
    print(f"    Cross:  {100 * var_cross / var_total:.1f}%")
    
    # ========================================================================
    # Analysis: Top contributions
    # ========================================================================
    print("\n" + "="*70)
    print("Top Contributing Terms")
    print("="*70)
    
    # Top linear terms
    print("\nTop 15 linear terms (κ_k):")
    top_linear_idx = np.argsort(kappa_linear)[-15:][::-1]
    for rank, idx in enumerate(top_linear_idx, 1):
        if kappa_linear[idx] > 1e-15:
            print(f"  {rank:2d}. Mode {all_mode_indices[idx]:2d} ({all_mode_freqs[idx]:6.1f} cm⁻¹): κ = {kappa_linear[idx]:.6e}")
    
    # Top cross terms
    print("\nTop 15 cross terms (d_{kl}):")
    top_cross_idx = np.argsort(d_ij)[-15:][::-1]
    for rank, idx in enumerate(top_cross_idx, 1):
        if d_ij[idx] > 1e-15:
            i, j = cross_pairs[idx]
            freq_i = all_mode_freqs[i]
            freq_j = all_mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            if i == j:
                print(f"  {rank:2d}. Mode {all_mode_indices[i]:2d}² ({freq_i:6.1f} cm⁻¹): d = {d_ij[idx]:.6e}")
            else:
                print(f"  {rank:2d}. Mode {all_mode_indices[i]:2d}×{all_mode_indices[j]:2d} ({freq_i:5.1f}×{freq_j:5.1f}, Σ={sum_freq:.0f}, Δ={diff_freq:.0f}): d = {d_ij[idx]:.6e}")
    
    # ========================================================================
    # Fit quality at target frequencies (via correlation analysis)
    # ========================================================================
    print("\n" + "="*70)
    print("Cross Term Contributions at Target Frequencies")
    print(f"Target peaks: {target_freqs} cm⁻¹")
    print("="*70)
    
    target_tolerance = 30.0
    
    for target in target_freqs:
        print(f"\n{'='*50}")
        print(f"Target: {target} cm⁻¹")
        print(f"{'='*50}")
        
        # Find cross terms that contribute at this frequency
        contributions = []
        for pair_idx, (i, j) in enumerate(cross_pairs):
            if d_ij[pair_idx] < 1e-15:
                continue
                
            freq_i = all_mode_freqs[i]
            freq_j = all_mode_freqs[j]
            sum_freq = freq_i + freq_j
            diff_freq = abs(freq_i - freq_j)
            
            match_type = None
            match_freq = 0
            if abs(sum_freq - target) < target_tolerance:
                match_type = 'sum'
                match_freq = sum_freq
            elif abs(diff_freq - target) < target_tolerance:
                match_type = 'diff'
                match_freq = diff_freq
            
            if match_type:
                contributions.append({
                    'pair_idx': pair_idx,
                    'i': i, 'j': j,
                    'freq_i': freq_i, 'freq_j': freq_j,
                    'match_type': match_type,
                    'match_freq': match_freq,
                    'd_ij': d_ij[pair_idx]
                })
        
        contributions_sorted = sorted(contributions, key=lambda x: x['d_ij'], reverse=True)
        
        print(f"\n  Cross terms matching {target} cm⁻¹ ({len(contributions)} total):")
        print(f"  {'Rank':>4} {'Mode k':>7} {'Mode l':>7} {'ωk':>7} {'ωl':>7} {'Type':>5} {'Result':>8} {'d_kl':>12}")
        print(f"  {'-'*70}")
        
        for rank, c in enumerate(contributions_sorted[:15], 1):
            print(f"  {rank:4d} {all_mode_indices[c['i']]:7d} {all_mode_indices[c['j']]:7d} "
                  f"{c['freq_i']:7.1f} {c['freq_j']:7.1f} {c['match_type']:>5} "
                  f"{c['match_freq']:8.1f} {c['d_ij']:12.3e}")
    
    # ========================================================================
    # Visualization
    # ========================================================================
    print("\n" + "="*70)
    print("Creating Plots")
    print("="*70)
    
    # Figure 1: Time series comparison
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    
    t_plot = np.arange(min(2000, n_frames_fit)) * dt / 1000  # ps, first 2000 frames
    n_plot = len(t_plot)
    
    # Plot 1: Stage 1 (linear only)
    ax = axes[0]
    ax.plot(t_plot, y[:n_plot] * 27211.4, 'k-', alpha=0.7, linewidth=0.8, label='Actual ΔE(t)')
    ax.plot(t_plot, y_pred_linear[:n_plot] * 27211.4, 'b-', alpha=0.8, linewidth=1.0, label=f'Linear (R²={R2_linear:.4f})')
    ax.set_ylabel('ΔE (meV)', fontsize=11)
    ax.set_title('Stage 1: Linear Terms Only', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Stage 2 (linear + cross)
    ax = axes[1]
    ax.plot(t_plot, y[:n_plot] * 27211.4, 'k-', alpha=0.7, linewidth=0.8, label='Actual ΔE(t)')
    ax.plot(t_plot, y_pred_total[:n_plot] * 27211.4, 'r-', alpha=0.8, linewidth=1.0, label=f'Linear+Cross (R²={R2_total:.4f})')
    ax.set_ylabel('ΔE (meV)', fontsize=11)
    ax.set_title('Stage 2: Linear + Cross Terms', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Residuals
    ax = axes[2]
    ax.plot(t_plot, residuals_stage1[:n_plot] * 27211.4, 'b-', alpha=0.5, linewidth=0.8, label=f'After Stage 1')
    ax.plot(t_plot, residuals_stage2[:n_plot] * 27211.4, 'r-', alpha=0.7, linewidth=0.8, label=f'After Stage 2')
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Time (ps)', fontsize=11)
    ax.set_ylabel('Residual (meV)', fontsize=11)
    ax.set_title('Residuals', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('fit_time_domain_cross_timeseries.png', dpi=150, bbox_inches='tight')
    print("  Saved: fit_time_domain_cross_timeseries.png")
    plt.close()
    
    # Figure 2: Coefficient analysis
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Panel 1: Linear coefficients vs frequency
    ax1 = fig.add_subplot(gs[0, 0])
    active_linear_mask = kappa_linear > 1e-15
    ax1.scatter(all_mode_freqs[active_linear_mask], kappa_linear[active_linear_mask] * 27211.4, 
                s=30, alpha=0.6, color='blue')
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax1.set_ylabel('κ (meV)', fontsize=11)
    ax1.set_title(f'Linear Coefficients ({n_active_linear} active)', fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Top linear modes
    ax2 = fig.add_subplot(gs[0, 1])
    top_n = min(15, n_active_linear)
    top_idx = np.argsort(kappa_linear)[-top_n:][::-1]
    ax2.barh(range(top_n), kappa_linear[top_idx] * 27211.4, color='blue', alpha=0.7)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels([f"Mode {all_mode_indices[i]} ({all_mode_freqs[i]:.0f})" for i in top_idx], fontsize=8)
    ax2.set_xlabel('κ (meV)', fontsize=11)
    ax2.set_title('Top Linear Modes', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='x')
    
    # Panel 3: Cross coefficient distribution
    ax3 = fig.add_subplot(gs[0, 2])
    d_sorted = np.sort(d_ij[d_ij > 1e-15])[::-1]
    if len(d_sorted) > 0:
        ax3.semilogy(d_sorted * 27211.4, 'o-', markersize=3, lw=1, color='red')
    ax3.set_xlabel('Cross term index (sorted)', fontsize=11)
    ax3.set_ylabel('d_{kl} (meV)', fontsize=11)
    ax3.set_title(f'Cross Coefficients ({n_active_cross} active)', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3, which='both')
    
    # Panel 4: Top cross terms
    ax4 = fig.add_subplot(gs[1, 0])
    top_cross_n = min(15, n_active_cross)
    top_cross_idx_plot = np.argsort(d_ij)[-top_cross_n:][::-1]
    ax4.barh(range(top_cross_n), d_ij[top_cross_idx_plot] * 27211.4, color='red', alpha=0.7)
    ax4.set_yticks(range(top_cross_n))
    labels = []
    for idx in top_cross_idx_plot:
        i, j = cross_pairs[idx]
        if i == j:
            labels.append(f"{all_mode_indices[i]}² ({all_mode_freqs[i]:.0f})")
        else:
            labels.append(f"{all_mode_indices[i]}×{all_mode_indices[j]} ({all_mode_freqs[i]:.0f},{all_mode_freqs[j]:.0f})")
    ax4.set_yticklabels(labels, fontsize=7)
    ax4.set_xlabel('d_{kl} (meV)', fontsize=11)
    ax4.set_title('Top Cross Terms', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='x')
    
    # Panel 5: Scatter plot y vs y_pred
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.scatter(y * 27211.4, y_pred_linear * 27211.4, s=5, alpha=0.3, label='Linear', color='blue')
    ax5.scatter(y * 27211.4, y_pred_total * 27211.4, s=5, alpha=0.3, label='Linear+Cross', color='red')
    lim = max(abs(y.min()), abs(y.max())) * 27211.4 * 1.1
    ax5.plot([-lim, lim], [-lim, lim], 'k--', lw=0.5)
    ax5.set_xlabel('Actual ΔE (meV)', fontsize=11)
    ax5.set_ylabel('Predicted ΔE (meV)', fontsize=11)
    ax5.set_title('Predicted vs Actual', fontsize=11, fontweight='bold')
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3)
    ax5.set_aspect('equal')
    ax5.set_xlim(-lim, lim)
    ax5.set_ylim(-lim, lim)
    
    # Panel 6: R² breakdown
    ax6 = fig.add_subplot(gs[1, 2])
    stages = ['Linear\nOnly', 'Linear+\nCross']
    R2_values = [R2_linear, R2_total]
    colors = ['blue', 'red']
    bars = ax6.bar(stages, R2_values, color=colors, alpha=0.7, edgecolor='black')
    ax6.set_ylabel('R²', fontsize=11)
    ax6.set_title('Model Performance', fontsize=11, fontweight='bold')
    ax6.set_ylim(0, 1)
    for bar, r2 in zip(bars, R2_values):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{r2:.4f}', 
                ha='center', va='bottom', fontsize=10)
    ax6.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('fit_time_domain_cross_coefficients.png', dpi=150, bbox_inches='tight')
    print("  Saved: fit_time_domain_cross_coefficients.png")
    plt.close()
    
    # ========================================================================
    # Save results
    # ========================================================================
    print("\n" + "="*70)
    print("Saving Results")
    print("="*70)
    
    # Linear coefficients
    with open('kappa_time_domain_cross_linear.dat', 'w') as f:
        f.write(f"# Time-domain fitting: Linear coefficients (NNLS, κ_k ≥ 0)\n")
        f.write(f"# Fitted on {n_frames_fit} frames ({n_frames_fit * dt / 1000:.1f} ps)\n")
        f.write(f"# R² (linear only) = {R2_linear:.6f}\n")
        f.write(f"# Active modes: {n_active_linear}\n")
        f.write(f"# Mode_index  Frequency(cm⁻¹)  κ(a.u.)  κ(meV)  Active\n")
        for i in range(N_modes_total):
            active = 'Yes' if kappa_linear[i] > 1e-15 else 'No'
            f.write(f"{all_mode_indices[i]:4d}  {all_mode_freqs[i]:10.2f}  {kappa_linear[i]:16.8e}  {kappa_linear[i] * 27211.4:12.6f}  {active}\n")
    print("  Saved: kappa_time_domain_cross_linear.dat")
    
    # Cross coefficients
    with open('kappa_time_domain_cross_quadratic.dat', 'w') as f:
        f.write(f"# Time-domain fitting: Cross coefficients (NNLS, d_kl ≥ 0)\n")
        f.write(f"# Fitted on {n_frames_fit} frames ({n_frames_fit * dt / 1000:.1f} ps)\n")
        f.write(f"# R² (total) = {R2_total:.6f}\n")
        f.write(f"# R² improvement from cross terms: +{R2_total - R2_linear:.6f}\n")
        f.write(f"# Active cross terms: {n_active_cross}/{n_cross}\n")
        f.write(f"# Mode_i  Freq_i  Mode_j  Freq_j  d_ij(a.u.)  d_ij(meV)\n")
        
        # Sort by magnitude
        sorted_idx = np.argsort(d_ij)[::-1]
        for idx in sorted_idx:
            if d_ij[idx] > 1e-15:
                i, j = cross_pairs[idx]
                f.write(f"{all_mode_indices[i]:4d}  {all_mode_freqs[i]:7.2f}  "
                       f"{all_mode_indices[j]:4d}  {all_mode_freqs[j]:7.2f}  "
                       f"{d_ij[idx]:16.8e}  {d_ij[idx] * 27211.4:12.6f}\n")
    print("  Saved: kappa_time_domain_cross_quadratic.dat")
    
    # Combined results (npz)
    np.savez('fit_time_domain_cross_results.npz',
             mode_indices=all_mode_indices,
             mode_freqs=all_mode_freqs,
             kappa_linear=kappa_linear,
             d_ij=d_ij,
             cross_pairs=np.array(cross_pairs),
             R2_linear=R2_linear,
             R2_total=R2_total,
             n_frames=n_frames_fit,
             dt=dt)
    print("  Saved: fit_time_domain_cross_results.npz")
    
    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTwo-Stage Time-Domain Fitting:")
    print(f"  Stage 1 (Linear NNLS):    R²={R2_linear:.4f}, {n_active_linear} active modes")
    print(f"  Stage 2 (Quadratic NNLS): ΔR²=+{R2_total - R2_linear:.4f}, {n_active_cross} active cross terms")
    print(f"  Final: R²={R2_total:.4f}")
    print(f"\n  Variance captured:")
    print(f"    Linear: {100 * var_linear / var_total:.1f}%")
    print(f"    Cross:  {100 * var_cross / var_total:.1f}%")
    print(f"    Total:  {100 * (var_linear + var_cross) / var_total:.1f}%")

else:
    print("\nNot enough modes in frequency range for cross terms!")
    print("Only linear fitting was performed.")

print("\n" + "="*70)
print("Done!")
print("="*70)
