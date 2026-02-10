#!/usr/bin/env python3
"""
Time-domain linear fitting with spectral density analysis:
  ΔE(t) = Σ_k κ_k Δq_k(t)

Fit coupling coefficients κ_k using NNLS (non-negative constraint).
Then compute and compare spectral densities in frequency domain.
"""
import numpy as np
from scipy.optimize import nnls
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
n_frames_fit = 1000  # frames for time-domain fitting (1000 * 4 fs = 4 ps)
max_lag_spec = 5000  # frames for spectral density (5000 * 4 fs = 20 ps)

print("="*70)
print("Time-Domain Linear Fitting with Spectral Density Analysis")
print("="*70)

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

N_frames = min(len(energy_diff_full), len(q_k_all))
N_modes = len(all_mode_freqs)

energy_diff = energy_diff_full[:N_frames]
mode_coords = q_k_all[:N_frames, :]

print(f"  Frames: {N_frames}")
print(f"  Modes: {N_modes}")
print(f"  Frequency range: {all_mode_freqs.min():.1f} - {all_mode_freqs.max():.1f} cm⁻¹")
print(f"  Time-domain fitting: {n_frames_fit} frames ({n_frames_fit * dt / 1000:.1f} ps)")
print(f"  Spectral density from: {max_lag_spec} frames ({max_lag_spec * dt / 1000:.1f} ps)")

# ============================================================================
# Filter modes: 450 ≤ freq ≤ 1800 cm⁻¹
# ============================================================================
freq_mask = (all_mode_freqs >= 450) & (all_mode_freqs <= 1800)
print(f"\nFiltered modes (450-1800 cm⁻¹): {np.sum(freq_mask)}")

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
# Linear fitting with OLS (no constraint)
# ============================================================================
print("\n" + "="*70)
print("Linear Fitting (OLS, no constraint on κ)")
print("="*70)

# Use subset of frames for fitting
X_fit = mode_coords_centered[:n_frames_fit, freq_mask]
y_fit = energy_diff_centered[:n_frames_fit]

# Design matrix (filtered modes)
X = X_fit
y = y_fit

# OLS fit (allows both positive and negative κ)
print("\nFitting with OLS (least squares)...")
XTX = X.T @ X
XTy = X.T @ y
kappa_filtered = np.linalg.solve(XTX, XTy)

# Expand to full mode array
kappa = np.zeros(N_modes)
kappa[freq_mask] = kappa_filtered

# Count positive/negative
n_positive = np.sum(kappa > 1e-15)
n_negative = np.sum(kappa < -1e-15)
print(f"  Positive κ: {n_positive}, Negative κ: {n_negative}")

# Predictions on fitting subset
y_pred = X_fit @ kappa_filtered

# Metrics (on fitting subset)
residuals = y_fit - y_pred
ss_tot = np.sum(y_fit**2)
ss_res = np.sum(residuals**2)
R2 = 1 - ss_res / ss_tot
rmse = np.sqrt(np.mean(residuals**2))

n_active = np.sum(np.abs(kappa) > 1e-15)

print(f"\nResults:")
print(f"  R² = {R2:.4f}")
print(f"  RMSE: {rmse:.6e} a.u. ({rmse * 27211.4:.3f} meV)")
print(f"  Active modes: {n_active}/{np.sum(freq_mask)}")

# ============================================================================
# Correlation functions
# ============================================================================
print("\n" + "="*70)
print("Computing Correlation Functions")
print("="*70)

def compute_autocorr(x, max_lag):
    """Compute autocorrelation (no window)"""
    x_mean = np.mean(x)
    x_centered = x - x_mean
    N = len(x)
    
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(x_centered[:N-lag] * x_centered[lag:])
    
    return acf

# Compute predictions on full data
y_pred_full = mode_coords_centered @ kappa
residuals_full = energy_diff_centered - y_pred_full

print(f"\nComputing C_total(t) = <ΔE(0)ΔE(t)> ({max_lag_spec * dt / 1000:.0f} ps)...")
C_total = compute_autocorr(energy_diff_centered, max_lag_spec)

print(f"Computing C_model(t) = <ΔE_model(0)ΔE_model(t)> ({max_lag_spec * dt / 1000:.0f} ps)...")
C_model = compute_autocorr(y_pred_full, max_lag_spec)

print(f"Computing C_residual(t) = <residual(0)residual(t)> ({max_lag_spec * dt / 1000:.0f} ps)...")
C_residual = compute_autocorr(residuals_full, max_lag_spec)

# Individual mode contributions
print("Computing C_k(t) for each mode...")
C_k = np.zeros((max_lag_spec, N_modes))
for k in range(N_modes):
    if np.abs(kappa[k]) > 1e-15:
        C_k[:, k] = compute_autocorr(kappa[k] * mode_coords_centered[:, k], max_lag_spec)

# Sum of individual contributions (ignoring cross-correlations)
C_sum_individual = np.sum(C_k, axis=1)

print(f"  C_total(0) = {C_total[0]:.6e}")
print(f"  C_model(0) = {C_model[0]:.6e}")
print(f"  C_sum_individual(0) = {C_sum_individual[0]:.6e}")

# ============================================================================
# Spectral density
# ============================================================================
print("\n" + "="*70)
print("Computing Spectral Densities")
print("="*70)

# Frequency grid
freq_grid = np.linspace(420, 1800, 1381)
print(f"Frequency grid: {len(freq_grid)} points (420-1800 cm⁻¹)")

def compute_spectral_density(C, freq_grid, dt, T):
    """Compute spectral density from correlation function"""
    t = np.arange(len(C)) * dt
    
    J = np.zeros(len(freq_grid))
    for i, nu in enumerate(freq_grid):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, x=t)
    
    # Multiply by prefactor: ω/(kB*T)
    J *= (2 * np.pi * c_cm_fs * freq_grid) / (kB_cm * T)
    return J

print("\nComputing J_total(ω)...")
J_total = compute_spectral_density(C_total, freq_grid, dt, T)

print("Computing J_model(ω)...")
J_model = compute_spectral_density(C_model, freq_grid, dt, T)

print("Computing J_residual(ω)...")
J_residual = compute_spectral_density(C_residual, freq_grid, dt, T)

print("Computing J_k(ω) for each mode...")
J_k = np.zeros((len(freq_grid), N_modes))
for k in range(N_modes):
    if np.abs(kappa[k]) > 1e-15:
        J_k[:, k] = compute_spectral_density(C_k[:, k], freq_grid, dt, T)

J_sum_individual = np.sum(J_k, axis=1)

# Integrals
J_total_int = np.trapezoid(J_total, freq_grid)
J_model_int = np.trapezoid(J_model, freq_grid)
J_residual_int = np.trapezoid(J_residual, freq_grid)
J_sum_int = np.trapezoid(J_sum_individual, freq_grid)

print(f"\nSpectral density integrals:")
print(f"  ∫J_total = {J_total_int:.6e}")
print(f"  ∫J_model = {J_model_int:.6e} ({100*J_model_int/J_total_int:.1f}%)")
print(f"  ∫J_sum_individual = {J_sum_int:.6e} ({100*J_sum_int/J_total_int:.1f}%)")
print(f"  ∫J_residual = {J_residual_int:.6e} ({100*J_residual_int/J_total_int:.1f}%)")

# R² in frequency domain
ss_tot_freq = np.sum((J_total - np.mean(J_total))**2)
ss_res_freq = np.sum((J_total - J_model)**2)
R2_freq = 1 - ss_res_freq / ss_tot_freq

print(f"\nR² in frequency domain: {R2_freq:.4f}")

# ============================================================================
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Creating Plots")
print("="*70)

# Figure 1: Time domain analysis
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Time series (use full data)
ax = axes[0, 0]
t_plot = np.arange(min(2000, N_frames)) * dt / 1000  # ps
n_plot = len(t_plot)
ax.plot(t_plot, energy_diff_centered[:n_plot] * 27211.4, 'k-', alpha=0.7, lw=0.8, label='Actual ΔE(t)')
ax.plot(t_plot, y_pred_full[:n_plot] * 27211.4, 'r-', alpha=0.8, lw=1.0, label=f'Model (R²={R2:.4f})')
ax.set_xlabel('Time (ps)', fontsize=11)
ax.set_ylabel('ΔE (meV)', fontsize=11)
ax.set_title('Time Series: ΔE(t)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Residuals
ax = axes[0, 1]
ax.plot(t_plot, residuals_full[:n_plot] * 27211.4, 'b-', alpha=0.7, lw=0.8)
ax.axhline(0, color='k', linestyle='--', lw=0.5)
ax.set_xlabel('Time (ps)', fontsize=11)
ax.set_ylabel('Residual (meV)', fontsize=11)
ax.set_title('Residuals: ΔE - ΔE_model', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Panel 3: Correlation function
ax = axes[1, 0]
t_corr = np.arange(max_lag_spec) * dt / 1000  # ps
ax.plot(t_corr, C_total * 27211.4**2, 'k-', lw=1.5, label='C_total(t)', alpha=0.8)
ax.plot(t_corr, C_model * 27211.4**2, 'r-', lw=1.5, label='C_model(t)', alpha=0.8)
ax.plot(t_corr, C_sum_individual * 27211.4**2, 'b--', lw=1.5, label='Σ C_k(t)', alpha=0.8)
ax.set_xlabel('Time lag (ps)', fontsize=11)
ax.set_ylabel('C(t) (meV²)', fontsize=11)
ax.set_title(f'Correlation Functions ({max_lag_spec * dt / 1000:.0f} ps)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, t_corr[-1])

# Panel 4: Coupling coefficients
ax = axes[1, 1]
active_mask = np.abs(kappa) > 1e-15
colors = ['red' if k > 0 else 'blue' for k in kappa[active_mask]]
ax.scatter(all_mode_freqs[active_mask], kappa[active_mask] * 27211.4, s=30, alpha=0.6, c=colors)
ax.axhline(0, color='k', linestyle='--', lw=0.5)
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('κ (meV)', fontsize=11)
ax.set_title(f'Coupling Coefficients (pos:{n_positive}, neg:{n_negative})', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fit_time_domain_linear_time.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_time_domain_linear_time.png")
plt.close()

# Figure 2: Frequency domain (spectral density)
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Spectral density comparison
ax = axes[0, 0]
ax.plot(freq_grid, J_total, 'k-', lw=2, label='J_total(ω) (data)', alpha=0.8)
ax.plot(freq_grid, J_model, 'r-', lw=1.5, label=f'J_model(ω) (R²={R2_freq:.4f})', alpha=0.8)
ax.fill_between(freq_grid, 0, J_model, alpha=0.2, color='red')
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('J(ω) (cm⁻¹²)', fontsize=11)
ax.set_title('Spectral Density: J_total vs J_model', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Residual in frequency domain
ax = axes[0, 1]
J_diff = J_total - J_model
ax.plot(freq_grid, J_diff, 'b-', lw=1.5, alpha=0.8)
ax.axhline(0, color='k', linestyle='--', lw=0.5)
ax.fill_between(freq_grid, J_diff, 0, where=(J_diff > 0), alpha=0.3, color='blue')
ax.fill_between(freq_grid, J_diff, 0, where=(J_diff < 0), alpha=0.3, color='red')
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('J_total - J_model (cm⁻¹²)', fontsize=11)
ax.set_title('Spectral Density Residual', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Panel 3: Individual mode contributions (stacked)
ax = axes[1, 0]
# Get top contributing modes
J_k_integrals = np.array([np.trapezoid(J_k[:, k], freq_grid) for k in range(N_modes)])
top_modes = np.argsort(J_k_integrals)[-10:][::-1]

colors = plt.cm.tab10(np.linspace(0, 1, 10))
for i, k in enumerate(top_modes):
    if J_k_integrals[k] > 0:
        ax.plot(freq_grid, J_k[:, k], '-', lw=1.5, alpha=0.7, color=colors[i],
               label=f'Mode {all_mode_indices[k]} ({all_mode_freqs[k]:.0f} cm⁻¹)')

ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('J_k(ω) (cm⁻¹²)', fontsize=11)
ax.set_title('Top 10 Mode Contributions', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper right', ncol=2)
ax.grid(True, alpha=0.3)

# Panel 4: Decomposition
ax = axes[1, 1]
ax.plot(freq_grid, J_total, 'k-', lw=2, label='J_total', alpha=0.8)
ax.plot(freq_grid, J_sum_individual, 'g--', lw=1.5, label=f'Σ J_k ({100*J_sum_int/J_total_int:.0f}%)', alpha=0.8)
ax.plot(freq_grid, J_residual, 'b-', lw=1, label=f'J_residual ({100*J_residual_int/J_total_int:.0f}%)', alpha=0.6)
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('J(ω) (cm⁻¹²)', fontsize=11)
ax.set_title('Spectral Density Decomposition', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fit_time_domain_linear_freq.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_time_domain_linear_freq.png")
plt.close()

# Figure 3: Combined overview
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# Row 1: Time domain (use full data)
ax1 = fig.add_subplot(gs[0, :2])
ax1.plot(t_plot, energy_diff_centered[:n_plot] * 27211.4, 'k-', alpha=0.7, lw=0.8, label='Actual')
ax1.plot(t_plot, y_pred_full[:n_plot] * 27211.4, 'r-', alpha=0.8, lw=1.0, label='Model')
ax1.set_xlabel('Time (ps)', fontsize=11)
ax1.set_ylabel('ΔE (meV)', fontsize=11)
ax1.set_title(f'Time Domain: R² = {R2:.4f}', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

ax2 = fig.add_subplot(gs[0, 2])
top_n = min(10, n_active)
top_idx = np.argsort(np.abs(kappa))[-top_n:][::-1]
colors_bar = ['red' if kappa[i] > 0 else 'blue' for i in top_idx]
ax2.barh(range(top_n), kappa[top_idx] * 27211.4, color=colors_bar, alpha=0.7)
ax2.set_yticks(range(top_n))
ax2.set_yticklabels([f"Mode {all_mode_indices[i]} ({all_mode_freqs[i]:.0f})" for i in top_idx], fontsize=8)
ax2.axvline(0, color='k', linestyle='-', lw=0.5)
ax2.set_xlabel('κ (meV)', fontsize=11)
ax2.set_title('Top Modes', fontsize=11, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')

# Row 2: Correlation functions
ax3 = fig.add_subplot(gs[1, :2])
ax3.plot(t_corr, C_total * 27211.4**2, 'k-', lw=1.5, label='C_total(t)')
ax3.plot(t_corr, C_model * 27211.4**2, 'r-', lw=1.5, label='C_model(t)')
ax3.set_xlabel('Time lag (ps)', fontsize=11)
ax3.set_ylabel('C(t) (meV²)', fontsize=11)
ax3.set_title('Correlation Functions', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, t_corr[-1])

ax4 = fig.add_subplot(gs[1, 2])
active_mask_plot = np.abs(kappa) > 1e-15
colors_scatter = ['red' if k > 0 else 'blue' for k in kappa[active_mask_plot]]
ax4.scatter(all_mode_freqs[active_mask_plot], kappa[active_mask_plot] * 27211.4, s=25, alpha=0.6, c=colors_scatter)
ax4.axhline(0, color='k', linestyle='--', lw=0.5)
ax4.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax4.set_ylabel('κ (meV)', fontsize=11)
ax4.set_title('κ vs Frequency', fontsize=11, fontweight='bold')
ax4.grid(True, alpha=0.3)

# Row 3: Frequency domain
ax5 = fig.add_subplot(gs[2, :2])
ax5.plot(freq_grid, J_total, 'k-', lw=2, label='J_total(ω)')
ax5.plot(freq_grid, J_model, 'r-', lw=1.5, label=f'J_model(ω)')
ax5.fill_between(freq_grid, 0, J_model, alpha=0.2, color='red')
ax5.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax5.set_ylabel('J(ω) (cm⁻¹²)', fontsize=11)
ax5.set_title(f'Spectral Density: R² = {R2_freq:.4f}', fontsize=12, fontweight='bold')
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)

ax6 = fig.add_subplot(gs[2, 2])
ax6.bar(['Time\nDomain', 'Frequency\nDomain'], [R2, R2_freq], color=['blue', 'green'], alpha=0.7)
ax6.set_ylabel('R²', fontsize=11)
ax6.set_title('Model Performance', fontsize=11, fontweight='bold')
ax6.set_ylim(0, 1)
for i, r2 in enumerate([R2, R2_freq]):
    ax6.text(i, r2 + 0.02, f'{r2:.4f}', ha='center', va='bottom', fontsize=10)
ax6.grid(True, alpha=0.3, axis='y')

plt.suptitle('Time-Domain Linear Fitting: ΔE(t) = Σ_k κ_k Δq_k(t)', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('fit_time_domain_linear_overview.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_time_domain_linear_overview.png")
plt.close()

# ============================================================================
# Save results
# ============================================================================
print("\n" + "="*70)
print("Saving Results")
print("="*70)

# Coefficients
with open('kappa_time_domain_linear_ols.dat', 'w') as f:
    f.write(f"# Time-domain linear fitting (OLS, no constraint)\n")
    f.write(f"# Frames: {N_frames}\n")
    f.write(f"# R² (time domain) = {R2:.6f}\n")
    f.write(f"# R² (frequency domain) = {R2_freq:.6f}\n")
    f.write(f"# RMSE: {rmse:.6e} a.u. ({rmse * 27211.4:.3f} meV)\n")
    f.write(f"# Active modes: {n_active} (pos:{n_positive}, neg:{n_negative})\n")
    f.write(f"# Mode  Freq(cm⁻¹)  κ(a.u.)  κ(meV)  Sign\n")
    for i in range(N_modes):
        sign = '+' if kappa[i] > 0 else '-' if kappa[i] < 0 else '0'
        f.write(f"{all_mode_indices[i]:4d}  {all_mode_freqs[i]:10.2f}  {kappa[i]:16.8e}  {kappa[i]*27211.4:12.6f}  {sign}\n")
print("  Saved: kappa_time_domain_linear_ols.dat")

# Spectral densities
np.savetxt('spectral_density_time_domain_linear.dat',
           np.column_stack([freq_grid, J_total, J_model, J_total - J_model]),
           header=f'Time-domain linear fitting spectral density\nR² (time) = {R2:.6f}, R² (freq) = {R2_freq:.6f}\nFreq(cm⁻¹)  J_total  J_model  J_residual',
           fmt='%.6e')
print("  Saved: spectral_density_time_domain_linear.dat")

# All results in npz
np.savez('fit_time_domain_linear_results.npz',
         mode_indices=all_mode_indices,
         mode_freqs=all_mode_freqs,
         kappa=kappa,
         freq_grid=freq_grid,
         J_total=J_total,
         J_model=J_model,
         C_total=C_total,
         C_model=C_model,
         R2_time=R2,
         R2_freq=R2_freq,
         dt=dt)
print("  Saved: fit_time_domain_linear_results.npz")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nTime-Domain Linear Fitting: ΔE(t) = Σ_k κ_k Δq_k(t)")
print(f"  Constraint: NNLS (κ_k ≥ 0)")
print(f"  Active modes: {n_active}/{np.sum(freq_mask)}")
print(f"\nPerformance:")
print(f"  Time domain R²:      {R2:.4f}")
print(f"  Frequency domain R²: {R2_freq:.4f}")
print(f"  RMSE: {rmse * 27211.4:.3f} meV")
print(f"\nSpectral integrals:")
print(f"  ∫J_total: {J_total_int:.6e}")
print(f"  ∫J_model: {J_model_int:.6e} ({100*J_model_int/J_total_int:.1f}%)")

print("\n" + "="*70)
print("Done!")
print("="*70)
