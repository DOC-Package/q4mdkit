#!/usr/bin/env python3
"""
Build spectral density from fitted g_k (κ_k) coefficients and J_qq matrix.

Spectral density formula:
    J(ω) = Σ_k Σ_l g_k g_l J_qq,kl(ω)
         = g^T J_qq(ω) g

Where:
- g_k: coupling coefficients from time-domain fitting (kappa_time_domain_ols.dat)
- J_qq,kl(ω): mode coordinate spectral density matrix (spectral_cache_J_qq.npz)

This combines time-domain fitted couplings with frequency-domain spectral information.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*70)
print("Build Spectral Density from g_k and J_qq")
print("J(ω) = g^T J_qq(ω) g")
print("="*70)

# ============================================================================
# Load g_k (kappa) coefficients from time-domain fitting
# ============================================================================
print("\nLoading g_k coefficients...")

gk_data = []
with open('kappa_time_domain_ols.dat') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 3:
            mode_idx = int(parts[0])
            freq = float(parts[1])
            kappa = float(parts[2])
            gk_data.append([mode_idx, freq, kappa])

gk_data = np.array(gk_data)
gk_mode_indices = gk_data[:, 0].astype(int)
gk_freqs = gk_data[:, 1]
gk_values = gk_data[:, 2]

print(f"  Loaded {len(gk_values)} coupling coefficients")
print(f"  Mode range: {gk_mode_indices.min()} - {gk_mode_indices.max()}")
print(f"  Frequency range: {gk_freqs.min():.1f} - {gk_freqs.max():.1f} cm⁻¹")

# ============================================================================
# Load J_qq spectral density matrix
# ============================================================================
print("\nLoading J_qq spectral density matrix...")

J_qq_data = np.load('spectral_cache_J_qq.npz')
J_qq_full = J_qq_data['J_qq']  # shape: (N_modes_all, N_modes_all, N_freq)

print(f"  J_qq shape: {J_qq_full.shape}")

# Load mode information to get mapping
stats_data = []
with open('mode_coords_new_stats.txt') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            stats_data.append([int(parts[0]), float(parts[1])])

stats_data = np.array(stats_data)
all_mode_indices = stats_data[:, 0].astype(int)
all_mode_freqs = stats_data[:, 1]

# Frequency grid (should match the one used to compute J_qq)
# From spectral_fit_cross.py: freq_grid = np.linspace(420, 1800, 1381)
freq_grid = np.linspace(420, 1800, J_qq_full.shape[2])
print(f"  Frequency grid: {len(freq_grid)} points ({freq_grid[0]:.1f} - {freq_grid[-1]:.1f} cm⁻¹)")

# ============================================================================
# Map g_k modes to J_qq indices
# ============================================================================
print("\nMapping g_k modes to J_qq indices...")

# J_qq is computed for modes with 400 ≤ freq ≤ 1800 cm⁻¹
jqq_mask = (all_mode_freqs >= 400) & (all_mode_freqs <= 1800)
jqq_mode_indices = all_mode_indices[jqq_mask]
jqq_mode_freqs = all_mode_freqs[jqq_mask]

print(f"  J_qq modes: {len(jqq_mode_indices)} (indices {jqq_mode_indices.min()}-{jqq_mode_indices.max()})")

# Create mapping from mode index to J_qq index
mode_to_jqq_idx = {mode: i for i, mode in enumerate(jqq_mode_indices)}

# Check that all g_k modes are in J_qq
n_gk = len(gk_mode_indices)
gk_jqq_indices = []
valid_gk_mask = []

for i, mode in enumerate(gk_mode_indices):
    if mode in mode_to_jqq_idx:
        gk_jqq_indices.append(mode_to_jqq_idx[mode])
        valid_gk_mask.append(True)
    else:
        print(f"  Warning: Mode {mode} not found in J_qq, skipping")
        valid_gk_mask.append(False)

valid_gk_mask = np.array(valid_gk_mask)
gk_jqq_indices = np.array(gk_jqq_indices)
gk_values_valid = gk_values[valid_gk_mask]
gk_freqs_valid = gk_freqs[valid_gk_mask]
gk_modes_valid = gk_mode_indices[valid_gk_mask]

print(f"  Valid g_k modes: {len(gk_values_valid)}/{n_gk}")

# ============================================================================
# Build spectral density: J(ω) = g^T J_qq(ω) g
# ============================================================================
print("\nBuilding spectral density...")

N_freq = len(freq_grid)
n_modes_gk = len(gk_values_valid)

# Extract relevant submatrix of J_qq
J_qq_sub = np.zeros((n_modes_gk, n_modes_gk, N_freq))
for i, idx_i in enumerate(gk_jqq_indices):
    for j, idx_j in enumerate(gk_jqq_indices):
        J_qq_sub[i, j, :] = J_qq_full[idx_i, idx_j, :]

# Compute J(ω) = g^T J_qq(ω) g
J_reconstructed = np.zeros(N_freq)
for w in range(N_freq):
    J_reconstructed[w] = gk_values_valid @ J_qq_sub[:, :, w] @ gk_values_valid

print(f"  Spectral density computed at {N_freq} frequency points")

# ============================================================================
# Also compute diagonal-only approximation
# ============================================================================
J_diagonal_only = np.zeros(N_freq)
for w in range(N_freq):
    for i, idx_i in enumerate(gk_jqq_indices):
        J_diagonal_only[w] += gk_values_valid[i]**2 * J_qq_full[idx_i, idx_i, w]

# ============================================================================
# Load reference J_total for comparison
# ============================================================================
print("\nLoading reference J_total...")
try:
    J_total = np.load('spectral_cache_J_total.npy')
    has_reference = True
    print(f"  J_total shape: {J_total.shape}")
except FileNotFoundError:
    print("  J_total not found, skipping comparison")
    has_reference = False
    J_total = None

# ============================================================================
# Statistics
# ============================================================================
print("\n" + "="*70)
print("Results")
print("="*70)

print(f"\nReconstructed spectral density:")
print(f"  Max: {np.max(J_reconstructed):.6e}")
print(f"  Integral: {np.trapz(J_reconstructed, freq_grid):.6e}")

print(f"\nDiagonal-only approximation:")
print(f"  Max: {np.max(J_diagonal_only):.6e}")
print(f"  Integral: {np.trapz(J_diagonal_only, freq_grid):.6e}")

if has_reference:
    # Compute R² and RMSE
    mask = (freq_grid >= 420) & (freq_grid <= 1800)
    J_ref = J_total[mask] if len(J_total) > len(freq_grid) else J_total
    
    if len(J_ref) == len(J_reconstructed):
        SS_res = np.sum((J_ref - J_reconstructed)**2)
        SS_tot = np.sum((J_ref - np.mean(J_ref))**2)
        R2_full = 1 - SS_res / SS_tot
        RMSE_full = np.sqrt(np.mean((J_ref - J_reconstructed)**2))
        
        SS_res_diag = np.sum((J_ref - J_diagonal_only)**2)
        R2_diag = 1 - SS_res_diag / SS_tot
        RMSE_diag = np.sqrt(np.mean((J_ref - J_diagonal_only)**2))
        
        # Compute optimal scaling factor
        scale_factor = np.sum(J_ref * J_reconstructed) / np.sum(J_reconstructed**2)
        J_reconstructed_scaled = J_reconstructed * scale_factor
        J_diagonal_scaled = J_diagonal_only * scale_factor
        
        SS_res_scaled = np.sum((J_ref - J_reconstructed_scaled)**2)
        R2_scaled = 1 - SS_res_scaled / SS_tot
        
        print(f"\nComparison with reference J_total:")
        print(f"  Reference max: {np.max(J_ref):.6e}")
        print(f"  Reconstructed max: {np.max(J_reconstructed):.6e}")
        print(f"  Scale factor (J_total / J_reconstructed): {scale_factor:.6e}")
        print(f"\n  Without scaling:")
        print(f"    Full:     R² = {R2_full:.6f}")
        print(f"    Diagonal: R² = {R2_diag:.6f}")
        print(f"\n  With optimal scaling:")
        print(f"    Full:     R² = {R2_scaled:.6f}")
        print(f"  Reference integral: {np.trapz(J_ref, freq_grid):.6e}")
        print(f"  Reconstructed (scaled) integral: {np.trapz(J_reconstructed_scaled, freq_grid):.6e}")

# ============================================================================
# Save results
# ============================================================================
print("\nSaving results...")

# Save spectral density data
output_file = 'spectral_density_from_gk.dat'
with open(output_file, 'w') as f:
    f.write("# Spectral density reconstructed from g_k and J_qq\n")
    f.write("# J(ω) = g^T J_qq(ω) g\n")
    f.write(f"# g_k from: kappa_time_domain_ols.dat\n")
    f.write(f"# J_qq from: spectral_cache_J_qq.npz\n")
    f.write(f"# Number of modes: {n_modes_gk}\n")
    f.write("# Freq(cm⁻¹)  J_full(ω)  J_diagonal(ω)\n")
    for i, freq in enumerate(freq_grid):
        f.write(f"{freq:12.4f} {J_reconstructed[i]:15.8e} {J_diagonal_only[i]:15.8e}\n")
print(f"  -> {output_file}")

# Save as npz for further analysis
np.savez('spectral_density_from_gk.npz',
         freq_grid=freq_grid,
         J_reconstructed=J_reconstructed,
         J_diagonal_only=J_diagonal_only,
         gk_modes=gk_modes_valid,
         gk_freqs=gk_freqs_valid,
         gk_values=gk_values_valid)
print("  -> spectral_density_from_gk.npz")

# ============================================================================
# Plot
# ============================================================================
print("\nGenerating plot...")

fig, axes = plt.subplots(3, 1, figsize=(12, 12))

# Top panel: Scaled comparison
ax1 = axes[0]
if has_reference and len(J_ref) == len(J_reconstructed):
    ax1.plot(freq_grid, J_ref, 'k-', lw=2, label='Reference J_total', alpha=0.8)
    ax1.plot(freq_grid, J_reconstructed_scaled, 'b-', lw=1.5, label='Reconstructed (scaled)', alpha=0.8)
    ax1.plot(freq_grid, J_diagonal_scaled, 'g--', lw=1.5, label='Diagonal only (scaled)', alpha=0.7)
    ax1.set_title(f'Spectral Density Comparison (scaled, R²={R2_scaled:.4f})', fontsize=14)
else:
    ax1.plot(freq_grid, J_reconstructed, 'b-', lw=1.5, label='Reconstructed (full)', alpha=0.8)
    ax1.plot(freq_grid, J_diagonal_only, 'g--', lw=1.5, label='Diagonal only', alpha=0.8)
    ax1.set_title('Reconstructed Spectral Density', fontsize=14)
ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
ax1.set_ylabel('J(ω)', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(420, 1800)

# Middle panel: Unscaled reconstruction (to see shape)
ax2 = axes[1]
ax2.plot(freq_grid, J_reconstructed, 'b-', lw=1.5, label='Reconstructed', alpha=0.8)
ax2.plot(freq_grid, J_diagonal_only, 'g--', lw=1.5, label='Diagonal only', alpha=0.7)
ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
ax2.set_ylabel('g^T J_qq(ω) g', fontsize=12)
ax2.set_title('Reconstructed Spectral Density (unscaled)', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(420, 1800)

# Bottom panel: g_k values as bar plot
ax3 = axes[2]
colors = ['blue' if g > 0 else 'red' for g in gk_values_valid]
ax3.bar(gk_freqs_valid, np.abs(gk_values_valid), width=15, color=colors, alpha=0.7,
        edgecolor='black', linewidth=0.5)
ax3.set_xlabel('Mode Frequency (cm⁻¹)', fontsize=12)
ax3.set_ylabel('|g_k| (a.u.)', fontsize=12)
ax3.set_title('Coupling Coefficients g_k (blue=positive, red=negative)', fontsize=14)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(400, 1850)

plt.tight_layout()
plt.savefig('spectral_density_from_gk.png', dpi=150, bbox_inches='tight')
print("  -> spectral_density_from_gk.png")

# ============================================================================
# Contribution analysis
# ============================================================================
print("\n" + "="*70)
print("Mode Contribution Analysis")
print("="*70)

# Compute contribution of each mode to total spectral density
contributions = np.zeros(n_modes_gk)
for i in range(n_modes_gk):
    # Diagonal contribution
    contributions[i] = np.trapz(gk_values_valid[i]**2 * J_qq_full[gk_jqq_indices[i], gk_jqq_indices[i], :], freq_grid)

# Sort by contribution
sorted_idx = np.argsort(np.abs(contributions))[::-1]

print(f"\nTop 15 contributing modes (diagonal only):")
print(f"{'Rank':>5} {'Mode':>6} {'Freq':>10} {'g_k':>15} {'Contribution':>15} {'%':>8}")
print("-" * 65)
total_diag_int = np.trapz(J_diagonal_only, freq_grid)
for rank, idx in enumerate(sorted_idx[:15]):
    mode = gk_modes_valid[idx]
    freq = gk_freqs_valid[idx]
    gk = gk_values_valid[idx]
    contrib = contributions[idx]
    pct = 100 * contrib / total_diag_int if total_diag_int > 0 else 0
    print(f"{rank+1:5d} {mode:6d} {freq:10.2f} {gk:15.6e} {contrib:15.6e} {pct:8.2f}%")

print("\n" + "="*70)
print("Done!")
print("="*70)
