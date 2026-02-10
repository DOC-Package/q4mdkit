#!/usr/bin/env python3
"""
Direct time-domain fitting: ΔE(t) = Σ_k κ_k Δq_k(t)

Fit coupling coefficients κ_k using linear regression on time-domain data.
Both ΔE(t) and Δq_k(t) are mean-subtracted before fitting.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*70)
print("Direct Time-Domain Fitting")
print("="*70)

# Parameters
n_frames_fit = 250  # Use subset of time series for fitting (125 frames = 500 fs = 0.5 ps)
dt = 4.0  # fs

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
        if len(parts) >= 5:
            stats_data.append([int(parts[0]), float(parts[1])])  # Mode, Freq

stats_data = np.array(stats_data)
all_mode_indices = stats_data[:, 0].astype(int)
all_mode_freqs = stats_data[:, 1]

# Filter modes: 400 ≤ freq ≤ 1800 cm⁻¹
freq_mask = (all_mode_freqs >= 400) & (all_mode_freqs <= 1800)
mode_indices = all_mode_indices[freq_mask]
mode_freqs = all_mode_freqs[freq_mask]

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords_new.txt')
q_k_all = q_k_data[:, 1:]
# Map mode indices to column indices (mode 7 -> column 0, etc.)
col_indices = mode_indices - all_mode_indices[0]  # First mode index offset
mode_coords_full = q_k_all[:, col_indices]

N_frames_full = min(len(energy_diff_full), len(mode_coords_full))
N_modes = len(mode_indices)

print(f"  Total frames: {N_frames_full}")
print(f"  Number of modes: {N_modes}")
print(f"  Frequency range: {mode_freqs.min():.1f} - {mode_freqs.max():.1f} cm⁻¹")

# ============================================================================
# Select subset for fitting
# ============================================================================
n_frames_fit = min(n_frames_fit, N_frames_full)
print(f"\nUsing {n_frames_fit} frames for fitting ({n_frames_fit * dt / 1000:.1f} ps)")

energy_diff = energy_diff_full[:n_frames_fit]
mode_coords = mode_coords_full[:n_frames_fit, :]

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
print(f"  Mode coord mean range: [{mode_means.min():.6e}, {mode_means.max():.6e}]")
print(f"  Mode coord std range: [{np.std(mode_coords_centered, axis=0).min():.6e}, {np.std(mode_coords_centered, axis=0).max():.6e}]")

# ============================================================================
# Linear regression: ΔE(t) = Σ_k κ_k Δq_k(t)
# ============================================================================
print("\n" + "="*70)
print("Linear Regression (OLS)")
print("="*70)

# X: [n_frames, n_modes], y: [n_frames]
X = mode_coords_centered
y = energy_diff_centered

# Ordinary Least Squares
print("\nSolving normal equations...")
XTX = X.T @ X
XTy = X.T @ y
kappa_ols = np.linalg.solve(XTX, XTy)

# Predictions
y_pred_ols = X @ kappa_ols

# Metrics
residuals_ols = y - y_pred_ols
mse_ols = np.mean(residuals_ols**2)
rmse_ols = np.sqrt(mse_ols)
ss_tot = np.sum((y - np.mean(y))**2)
ss_res = np.sum(residuals_ols**2)
r2_ols = 1 - ss_res / ss_tot

print(f"\nOLS Results:")
print(f"  RMSE: {rmse_ols:.6e} a.u.")
print(f"  R²: {r2_ols:.6f}")
print(f"  κ range: [{kappa_ols.min():.6e}, {kappa_ols.max():.6e}]")
print(f"  |κ| mean: {np.mean(np.abs(kappa_ols)):.6e}")

# ============================================================================
# Save results
# ============================================================================
print("\n" + "="*70)
print("Saving Results")
print("="*70)

# Save OLS coefficients
with open('kappa_time_domain_ols.dat', 'w') as f:
    f.write("# Direct time-domain fitting (OLS): ΔE(t) = Σ_k κ_k Δq_k(t)\n")
    f.write(f"# Fitted on {n_frames_fit} frames ({n_frames_fit * dt / 1000:.1f} ps)\n")
    f.write(f"# RMSE: {rmse_ols:.6e} a.u., R²: {r2_ols:.6f}\n")
    f.write("# Mode_index  Frequency(cm⁻¹)  κ(a.u.)\n")
    for idx, freq, k in zip(mode_indices, mode_freqs, kappa_ols):
        f.write(f"{idx:4d}  {freq:10.2f}  {k:16.8e}\n")
print("  Saved: kappa_time_domain_ols.dat")

# ============================================================================
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Creating Plots")
print("="*70)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Time series comparison
ax = axes[0]
t_plot = np.arange(n_frames_fit) * dt / 1000  # ps
ax.plot(t_plot, y * 27211.4, 'k-', alpha=0.7, linewidth=1.0, label='Actual ΔE(t)')
ax.plot(t_plot, y_pred_ols * 27211.4, 'r-', alpha=0.8, linewidth=1.2, label='OLS Fit')
ax.set_xlabel('Time (ps)', fontsize=11)
ax.set_ylabel('ΔE (meV)', fontsize=11)
ax.set_title(f'Time-Domain Fit (R²={r2_ols:.4f})', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 2: Coupling coefficients
ax = axes[1]
ax.scatter(mode_freqs, kappa_ols * 27211.4, s=30, alpha=0.6, color='red')
ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('κ (meV)', fontsize=11)
ax.set_title('Coupling Coefficients vs Frequency', fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fit_time_domain_direct.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_time_domain_direct.png")

# ============================================================================
# Compute spectral density from time-domain coefficients
# ============================================================================
print("\n" + "="*70)
print("Computing Spectral Density from Time-Domain Coefficients")
print("="*70)

# Convert frequency from cm^-1 to angular frequency (rad/fs)
# ω = 2πc * ν, where c = 2.998e10 cm/s = 2.998e-5 cm/fs
c_cm_per_fs = 2.99792458e-5  # cm/fs
omega_k = 2 * np.pi * c_cm_per_fs * mode_freqs  # rad/fs

# ℏ = 2.4188843e-5 a.u.·fs
hbar_au_fs = 2.4188843e-5  # a.u.·fs

# Spectral density contribution for each mode
# In the linear coupling model: ΔE = Σ_k κ_k Δq_k
# The reorganization energy is: λ_k = κ_k² * <Δq_k²> / (2 * ℏω_k)
# where <Δq_k²> is the mean square displacement

# Compute mean square displacement for each mode
mean_sq_disp = np.var(mode_coords_centered, axis=0)  # Already centered

# Reorganization energy
lambda_k = (kappa_ols**2 * mean_sq_disp) / (2 * hbar_au_fs * omega_k)  # a.u.
lambda_k_meV = lambda_k * 27211.4  # meV

# Huang-Rhys factor: S_k = λ_k / (ℏω_k)
S_k = lambda_k / (hbar_au_fs * omega_k)

# Total reorganization energy
lambda_total = np.sum(lambda_k)
lambda_total_meV = lambda_total * 27211.4

print(f"\nSpectral density statistics:")
print(f"  Huang-Rhys factors:")
print(f"    Range: [{S_k.min():.6e}, {S_k.max():.6e}]")
print(f"    Sum: {np.sum(S_k):.6f}")
print(f"  Reorganization energies:")
print(f"    Range: [{lambda_k_meV.min():.3f}, {lambda_k_meV.max():.3f}] meV")
print(f"    Total: {lambda_total_meV:.3f} meV ({lambda_total*1000:.3f} cm⁻¹)")

# Save spectral density
with open('spectral_density_time_domain.dat', 'w') as f:
    f.write("# Spectral density from time-domain fitting\n")
    f.write(f"# Total reorganization energy: {lambda_total_meV:.3f} meV\n")
    f.write(f"# Total Huang-Rhys factor: {np.sum(S_k):.6f}\n")
    f.write("# Mode_index  Frequency(cm⁻¹)  κ(meV)  S_k  λ_k(meV)\n")
    for idx, freq, k, s, lam in zip(mode_indices, mode_freqs, kappa_ols * 27211.4, S_k, lambda_k_meV):
        f.write(f"{idx:4d}  {freq:10.2f}  {k:12.6f}  {s:12.6e}  {lam:12.6f}\n")
print("  Saved: spectral_density_time_domain.dat")

# ============================================================================
# Comparison with frequency-domain results
# ============================================================================
print("\n" + "="*70)
print("Comparison with Frequency-Domain Results")
print("="*70)

try:
    # Try to load previous frequency-domain results
    freq_results = np.loadtxt('g_k_quadratic_form_twostage.dat')
    g_k_freq = freq_results[:, 2]  # Assuming column 2 is g_k
    
    # Compute spectral density from frequency-domain results for comparison
    # Assuming g_k is related to reorganization energy
    lambda_freq = g_k_freq  # May need scaling depending on definition
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Coupling coefficients
    ax = axes[0, 0]
    ax.scatter(mode_freqs, kappa_ols * 27211.4, s=30, alpha=0.6, label='Time-domain', color='red')
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('κ (meV)', fontsize=11)
    ax.set_title('Coupling Coefficients', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Huang-Rhys factors
    ax = axes[0, 1]
    ax.scatter(mode_freqs, S_k, s=30, alpha=0.6, color='red')
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('Huang-Rhys factor S_k', fontsize=11)
    ax.set_title('Huang-Rhys Factors', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Reorganization energies comparison
    ax = axes[1, 0]
    ax.scatter(mode_freqs, lambda_k_meV, s=30, alpha=0.6, label='Time-domain', color='red')
    ax.scatter(mode_freqs, g_k_freq, s=20, alpha=0.6, label='Frequency-domain', color='green', marker='x')
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('Reorganization energy (meV)', fontsize=11)
    ax.set_title('Reorganization Energy Comparison', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Spectral density
    ax = axes[1, 1]
    # J(ω) = π/2 * Σ_k λ_k δ(ω - ω_k)
    # For discrete representation, show λ_k at each frequency
    ax.stem(mode_freqs, lambda_k_meV, linefmt='r-', markerfmt='ro', basefmt=' ', label='Time-domain')
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('λ_k (meV)', fontsize=11)
    ax.set_title(f'Spectral Density (Total λ = {lambda_total_meV:.1f} meV)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('fit_time_vs_frequency_domain.png', dpi=150, bbox_inches='tight')
    print("  Saved: fit_time_vs_frequency_domain.png")
    
    # Compute and print comparison statistics
    print(f"\n  Comparison statistics:")
    print(f"    Time-domain total λ: {lambda_total_meV:.3f} meV")
    print(f"    Frequency-domain total: {np.sum(g_k_freq):.3f} meV")
    print(f"    Ratio: {lambda_total_meV / np.sum(g_k_freq):.4f}")
    
except FileNotFoundError:
    print("  Frequency-domain results not found, skipping comparison")
    
    # Plot time-domain results only
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot 1: Coupling coefficients
    ax = axes[0]
    ax.scatter(mode_freqs, kappa_ols * 27211.4, s=30, alpha=0.6, color='red')
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('κ (meV)', fontsize=11)
    ax.set_title('Coupling Coefficients', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Spectral density
    ax = axes[1]
    ax.stem(mode_freqs, lambda_k_meV, linefmt='r-', markerfmt='ro', basefmt=' ')
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('λ_k (meV)', fontsize=11)
    ax.set_title(f'Spectral Density (Total λ = {lambda_total_meV:.1f} meV)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('spectral_density_time_domain.png', dpi=150, bbox_inches='tight')
    print("  Saved: spectral_density_time_domain.png")
    
except FileNotFoundError:
    print("  Frequency-domain results not found, skipping comparison")

print("\n" + "="*70)
print("Done!")
print("="*70)
