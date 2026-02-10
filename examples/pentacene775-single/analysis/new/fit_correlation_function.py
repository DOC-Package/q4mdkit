#!/usr/bin/env python3
"""
Time-domain correlation function fitting:
  C_total(t) = Σ_k c_k * C_k(t)

where:
  C_total(t) = <ΔE(0)ΔE(t)>
  C_k(t) = <Δq_k(0)Δq_k(t)>

This is the time-domain equivalent of spectral_fit_with_cross_terms_v3.py
(since J(ω) is the Fourier transform of C(t)).
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
max_lag_fit = 1000  # frames for fitting (1000 * 4 fs = 4 ps)
max_lag_spec = 5000 # frames for spectral density (5000 * 4 fs = 20 ps)

print("="*70)
print("Time-Domain Correlation Function Fitting")
print("C_total(t) = Σ_k c_k * C_k(t)")
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
print(f"  Fitting correlation time: {max_lag_fit * dt / 1000:.2f} ps ({max_lag_fit} frames)")
print(f"  Spectral density correlation time: {max_lag_spec * dt / 1000:.2f} ps ({max_lag_spec} frames)")

# ============================================================================
# Filter modes: 450 ≤ freq ≤ 1800 cm⁻¹
# ============================================================================
freq_mask = (all_mode_freqs >= 450) & (all_mode_freqs <= 1800)
n_filtered = np.sum(freq_mask)
print(f"\nFiltered modes (450-1800 cm⁻¹): {n_filtered}")

# ============================================================================
# Subtract means
# ============================================================================
print("\nSubtracting means...")

energy_mean = np.mean(energy_diff)
energy_diff_centered = energy_diff - energy_mean

mode_means = np.mean(mode_coords, axis=0)
mode_coords_centered = mode_coords - mode_means

print(f"  Energy std: {np.std(energy_diff_centered):.6e} a.u. ({np.std(energy_diff_centered)*27211.4:.3f} meV)")

# ============================================================================
# Compute correlation functions
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

print("\nComputing C_total(t) = <ΔE(0)ΔE(t)> for fitting (1 ps)...")
C_total_fit = compute_autocorr(energy_diff_centered, max_lag_fit)

print("Computing C_k(t) = <Δq_k(0)Δq_k(t)> for fitting (1 ps)...")
C_k_fit = np.zeros((max_lag_fit, N_modes))
for k in range(N_modes):
    C_k_fit[:, k] = compute_autocorr(mode_coords_centered[:, k], max_lag_fit)
    if (k + 1) % 20 == 0:
        print(f"  Progress: {k+1}/{N_modes}")

print(f"\n  C_total(0) = {C_total_fit[0]:.6e} a.u.² ({C_total_fit[0] * 27211.4**2:.3f} meV²)")

# ============================================================================
# Fit correlation function with NNLS (c_k ≥ 0)
# ============================================================================
print("\n" + "="*70)
print("Fitting Correlation Function (NNLS, c_k ≥ 0)")
print("="*70)

# Use only filtered modes
C_k_filtered = C_k_fit[:, freq_mask]

# Normalize for better conditioning
C_k_integrals = np.trapezoid(C_k_filtered, axis=0) * dt
C_total_integral_fit = np.trapezoid(C_total_fit) * dt

# NNLS fit: C_total(t) ≈ Σ_k c_k * C_k(t)
print("\nFitting with NNLS...")
c_k_filtered, residual_norm = nnls(C_k_filtered, C_total_fit)

# Expand to full mode array
c_k = np.zeros(N_modes)
c_k[freq_mask] = c_k_filtered

# Reconstruct (for fitting range)
C_model_fit = C_k_fit @ c_k
C_residual_fit = C_total_fit - C_model_fit

# Metrics for fitting range
C_model_integral_fit = np.trapezoid(C_model_fit) * dt
C_residual_integral_fit = np.trapezoid(C_residual_fit) * dt

ss_tot = np.sum(C_total_fit**2)
ss_res = np.sum(C_residual_fit**2)
R2_fit = 1 - ss_res / ss_tot

n_active = np.sum(c_k > 1e-15)

print(f"\nFitting Results (0-1 ps):")
print(f"  R² = {R2_fit:.4f}")
print(f"  Active modes: {n_active}/{n_filtered}")
print(f"  C(0) capture: {100 * C_model_fit[0] / C_total_fit[0]:.1f}%")
print(f"  Integral capture: {100 * C_model_integral_fit / C_total_integral_fit:.1f}%")

# Top contributing modes
print("\nTop 10 modes by c_k:")
top_idx = np.argsort(c_k)[-10:][::-1]
for rank, idx in enumerate(top_idx, 1):
    if c_k[idx] > 1e-15:
        print(f"  {rank:2d}. Mode {all_mode_indices[idx]:2d} ({all_mode_freqs[idx]:6.1f} cm⁻¹): c_k = {c_k[idx]:.6e}")

# ============================================================================
# Compute long correlation functions for high-resolution spectral density
# ============================================================================
print("\n" + "="*70)
print(f"Computing Long Correlation Functions ({max_lag_spec * dt / 1000:.0f} ps)")
print("="*70)

print("\nComputing C_total(t) for spectral density...")
C_total_spec = compute_autocorr(energy_diff_centered, max_lag_spec)

print("Computing C_k(t) for spectral density...")
C_k_spec = np.zeros((max_lag_spec, N_modes))
for k in range(N_modes):
    C_k_spec[:, k] = compute_autocorr(mode_coords_centered[:, k], max_lag_spec)
    if (k + 1) % 20 == 0:
        print(f"  Progress: {k+1}/{N_modes}")

# Reconstruct model using fitted coefficients
C_model_spec = C_k_spec @ c_k
C_residual_spec = C_total_spec - C_model_spec

# Check fit quality over full range
ss_tot_spec = np.sum(C_total_spec**2)
ss_res_spec = np.sum(C_residual_spec**2)
R2_spec = 1 - ss_res_spec / ss_tot_spec

print(f"\nFull range reconstruction (0-{max_lag_spec * dt / 1000:.0f} ps):")
print(f"  R² = {R2_spec:.4f}")

# ============================================================================
# Compute spectral densities from long correlation functions
# ============================================================================
print("\n" + "="*70)
print("Computing High-Resolution Spectral Densities")
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

print(f"\nComputing J_total(ω) from {max_lag_spec * dt / 1000:.0f} ps correlation...")
J_total = compute_spectral_density(C_total_spec, freq_grid, dt, T)

print(f"Computing J_model(ω) from {max_lag_spec * dt / 1000:.0f} ps correlation...")
J_model = compute_spectral_density(C_model_spec, freq_grid, dt, T)

print(f"Computing J_residual(ω)...")
J_residual = compute_spectral_density(C_residual_spec, freq_grid, dt, T)

# Individual mode spectral densities (weighted by c_k)
print("Computing J_k(ω) for each active mode...")
J_k = np.zeros((len(freq_grid), N_modes))
for k in range(N_modes):
    if c_k[k] > 1e-15:
        J_k[:, k] = c_k[k] * compute_spectral_density(C_k_spec[:, k], freq_grid, dt, T)

J_sum_individual = np.sum(J_k, axis=1)

# Integrals
J_total_int = np.trapezoid(J_total, freq_grid)
J_model_int = np.trapezoid(J_model, freq_grid)
J_residual_int = np.trapezoid(J_residual, freq_grid)

print(f"\nSpectral density integrals (from {max_lag_spec * dt / 1000:.0f} ps correlation):")
print(f"  ∫J_total = {J_total_int:.6e}")
print(f"  ∫J_model = {J_model_int:.6e} ({100*J_model_int/J_total_int:.1f}%)")
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

# Figure 1: Time domain (correlation functions) - Fitting range
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

t_corr_fit = np.arange(max_lag_fit) * dt / 1000  # ps

# Panel 1: Correlation function comparison (fitting range)
ax = axes[0, 0]
ax.plot(t_corr_fit, C_total_fit * 27211.4**2, 'k-', lw=2, label='C_total(t)', alpha=0.8)
ax.plot(t_corr_fit, C_model_fit * 27211.4**2, 'r--', lw=1.5, label=f'C_model(t) (R²={R2_fit:.4f})', alpha=0.8)
ax.set_xlabel('Time lag (ps)', fontsize=11)
ax.set_ylabel('C(t) (meV²)', fontsize=11)
ax.set_title('Correlation Function Fit (0-1 ps)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Residual correlation (fitting range)
ax = axes[0, 1]
ax.plot(t_corr_fit, C_residual_fit * 27211.4**2, 'b-', lw=1.5, alpha=0.8)
ax.axhline(0, color='k', linestyle='--', lw=0.5)
ax.set_xlabel('Time lag (ps)', fontsize=11)
ax.set_ylabel('C_total - C_model (meV²)', fontsize=11)
ax.set_title(f'Residual (fitting range)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Panel 3: Full range correlation
ax = axes[1, 0]
t_corr_spec = np.arange(max_lag_spec) * dt / 1000  # ps
ax.plot(t_corr_spec, C_total_spec * 27211.4**2, 'k-', lw=1.5, label='C_total(t)', alpha=0.8)
ax.plot(t_corr_spec, C_model_spec * 27211.4**2, 'r-', lw=1, label=f'C_model(t) (R²={R2_spec:.4f})', alpha=0.8)
ax.axvline(max_lag_fit * dt / 1000, color='blue', linestyle='--', lw=1, label='Fit range (1 ps)')
ax.set_xlabel('Time lag (ps)', fontsize=11)
ax.set_ylabel('C(t) (meV²)', fontsize=11)
ax.set_title(f'Full Range Correlation (0-{max_lag_spec * dt / 1000:.0f} ps)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 4: Coefficients vs frequency
ax = axes[1, 1]
active_mask = c_k > 1e-15
ax.scatter(all_mode_freqs[active_mask], c_k[active_mask], s=40, alpha=0.6, color='red')
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('c_k', fontsize=11)
ax.set_title(f'Fitting Coefficients ({n_active} active)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fit_correlation_time.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_correlation_time.png")
plt.close()

# Figure 2: Frequency domain (spectral density)
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Spectral density comparison
ax = axes[0, 0]
ax.plot(freq_grid, J_total, 'k-', lw=2, label='J_total(ω)', alpha=0.8)
ax.plot(freq_grid, J_model, 'r--', lw=1.5, label=f'J_model(ω) (R²={R2_freq:.4f})', alpha=0.8)
ax.fill_between(freq_grid, 0, J_model, alpha=0.2, color='red')
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('J(ω) (cm⁻¹²)', fontsize=11)
ax.set_title(f'High-Res Spectral Density (from {max_lag_spec * dt / 1000:.0f} ps corr)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Residual spectral density
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

# Panel 3: Individual mode spectral densities
ax = axes[1, 0]
colors = plt.cm.tab10(np.linspace(0, 1, 10))
for i, idx in enumerate(top_idx[:10]):
    if c_k[idx] > 1e-15:
        ax.plot(freq_grid, J_k[:, idx], '-', lw=1.5, alpha=0.7, color=colors[i],
               label=f'Mode {all_mode_indices[idx]} ({all_mode_freqs[idx]:.0f})')
ax.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax.set_ylabel('c_k × J_k(ω) (cm⁻¹²)', fontsize=11)
ax.set_title('Top Mode Spectral Contributions', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

# Panel 4: Top modes bar chart
ax = axes[1, 1]
top_n = min(15, n_active)
top_idx_plot = np.argsort(c_k)[-top_n:][::-1]
ax.barh(range(top_n), c_k[top_idx_plot], color='red', alpha=0.7)
ax.set_yticks(range(top_n))
ax.set_yticklabels([f"Mode {all_mode_indices[i]} ({all_mode_freqs[i]:.0f})" for i in top_idx_plot], fontsize=8)
ax.set_xlabel('c_k', fontsize=11)
ax.set_title(f'Top {top_n} Modes', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('fit_correlation_freq.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_correlation_freq.png")
plt.close()

# Figure 3: Combined overview
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# Panel 1: Time domain fit (fitting range)
ax1 = fig.add_subplot(gs[0, :2])
ax1.plot(t_corr_fit, C_total_fit * 27211.4**2, 'k-', lw=2, label='C_total(t) (data)', alpha=0.8)
ax1.plot(t_corr_fit, C_model_fit * 27211.4**2, 'r-', lw=1.5, label=f'C_model(t) (fit)', alpha=0.8)
ax1.fill_between(t_corr_fit, 0, C_model_fit * 27211.4**2, alpha=0.2, color='red')
ax1.set_xlabel('Time lag (ps)', fontsize=11)
ax1.set_ylabel('C(t) (meV²)', fontsize=11)
ax1.set_title(f'Fitting Range (0-1 ps): R² = {R2_fit:.4f}', fontsize=12, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Panel 2: Top modes
ax2 = fig.add_subplot(gs[0, 2])
top_n_bar = min(10, n_active)
top_idx_bar = np.argsort(c_k)[-top_n_bar:][::-1]
ax2.barh(range(top_n_bar), c_k[top_idx_bar], color='red', alpha=0.7)
ax2.set_yticks(range(top_n_bar))
ax2.set_yticklabels([f"Mode {all_mode_indices[i]} ({all_mode_freqs[i]:.0f})" for i in top_idx_bar], fontsize=8)
ax2.set_xlabel('c_k', fontsize=11)
ax2.set_title('Top Modes', fontsize=11, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')

# Panel 3: Frequency domain
ax3 = fig.add_subplot(gs[1, :2])
ax3.plot(freq_grid, J_total, 'k-', lw=2, label='J_total(ω)', alpha=0.8)
ax3.plot(freq_grid, J_model, 'r-', lw=1.5, label=f'J_model(ω)', alpha=0.8)
ax3.fill_between(freq_grid, 0, J_model, alpha=0.2, color='red')
ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
ax3.set_ylabel('J(ω) (cm⁻¹²)', fontsize=11)
ax3.set_title(f'High-Res Spectral Density ({max_lag_spec * dt / 1000:.0f} ps): R² = {R2_freq:.4f}', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)

# Panel 4: Performance summary
ax4 = fig.add_subplot(gs[1, 2])
metrics = ['Fit\nR²', 'Full\nR²', 'Freq\nR²', 'C(0)\nCapture']
values = [R2_fit, R2_spec, R2_freq, C_model_fit[0]/C_total_fit[0]]
colors_bar = ['blue', 'cyan', 'green', 'orange']
bars = ax4.bar(metrics, values, color=colors_bar, alpha=0.7)
ax4.set_ylabel('Value', fontsize=11)
ax4.set_title('Performance Summary', fontsize=11, fontweight='bold')
ax4.set_ylim(0, 1.1)
for bar, val in zip(bars, values):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
             f'{val:.3f}', ha='center', va='bottom', fontsize=9)
ax4.grid(True, alpha=0.3, axis='y')

plt.suptitle(f'Correlation Fitting (1 ps) → High-Res Spectral Density ({max_lag_spec * dt / 1000:.0f} ps)', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('fit_correlation_overview.png', dpi=150, bbox_inches='tight')
print("  Saved: fit_correlation_overview.png")
plt.close()

# ============================================================================
# Save results
# ============================================================================
print("\n" + "="*70)
print("Saving Results")
print("="*70)

# Coefficients
with open('c_k_correlation_fit.dat', 'w') as f:
    f.write(f"# Correlation function fitting (NNLS, c_k ≥ 0)\n")
    f.write(f"# C_total(t) = Σ_k c_k × C_k(t)\n")
    f.write(f"# Fitting range: 0-{max_lag_fit * dt / 1000:.2f} ps ({max_lag_fit} frames)\n")
    f.write(f"# Spectral density from: 0-{max_lag_spec * dt / 1000:.0f} ps ({max_lag_spec} frames)\n")
    f.write(f"# R² (fit range) = {R2_fit:.6f}\n")
    f.write(f"# R² (full range) = {R2_spec:.6f}\n")
    f.write(f"# R² (frequency domain) = {R2_freq:.6f}\n")
    f.write(f"# Active modes: {n_active}\n")
    f.write(f"# Mode  Freq(cm⁻¹)  c_k  Active\n")
    for i in range(N_modes):
        active = 'Yes' if c_k[i] > 1e-15 else 'No'
        f.write(f"{all_mode_indices[i]:4d}  {all_mode_freqs[i]:10.2f}  {c_k[i]:16.8e}  {active}\n")
print("  Saved: c_k_correlation_fit.dat")

# Correlation functions (fitting range)
np.savetxt('correlation_functions_fit.dat',
           np.column_stack([t_corr_fit, C_total_fit * 27211.4**2, C_model_fit * 27211.4**2, C_residual_fit * 27211.4**2]),
           header=f'Correlation function fitting (0-1 ps)\nR² = {R2_fit:.6f}\nTime(ps)  C_total(meV²)  C_model(meV²)  C_residual(meV²)',
           fmt='%.6e')
print("  Saved: correlation_functions_fit.dat")

# Correlation functions (full range)
np.savetxt('correlation_functions_full.dat',
           np.column_stack([t_corr_spec, C_total_spec * 27211.4**2, C_model_spec * 27211.4**2, C_residual_spec * 27211.4**2]),
           header=f'Correlation function (0-{max_lag_spec * dt / 1000:.0f} ps)\nR² = {R2_spec:.6f}\nTime(ps)  C_total(meV²)  C_model(meV²)  C_residual(meV²)',
           fmt='%.6e')
print("  Saved: correlation_functions_full.dat")

# Spectral densities
np.savetxt('spectral_density_correlation_fit.dat',
           np.column_stack([freq_grid, J_total, J_model, J_total - J_model]),
           header=f'High-res spectral density from correlation fitting\nFit: 0-1ps, Spectral: 0-{max_lag_spec * dt / 1000:.0f}ps\nR² (freq) = {R2_freq:.6f}\nFreq(cm⁻¹)  J_total  J_model  J_residual',
           fmt='%.6e')
print("  Saved: spectral_density_correlation_fit.dat")

# All results in npz
np.savez('fit_correlation_results.npz',
         mode_indices=all_mode_indices,
         mode_freqs=all_mode_freqs,
         c_k=c_k,
         t_corr_fit=t_corr_fit,
         C_total_fit=C_total_fit,
         C_model_fit=C_model_fit,
         t_corr_spec=t_corr_spec,
         C_total_spec=C_total_spec,
         C_model_spec=C_model_spec,
         freq_grid=freq_grid,
         J_total=J_total,
         J_model=J_model,
         R2_fit=R2_fit,
         R2_spec=R2_spec,
         R2_freq=R2_freq,
         dt=dt,
         max_lag_fit=max_lag_fit,
         max_lag_spec=max_lag_spec)
print("  Saved: fit_correlation_results.npz")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nCorrelation Function Fitting: C_total(t) = Σ_k c_k × C_k(t)")
print(f"  Constraint: NNLS (c_k ≥ 0)")
print(f"  Fitting range: 0-{max_lag_fit * dt / 1000:.2f} ps")
print(f"  Spectral density from: 0-{max_lag_spec * dt / 1000:.0f} ps")
print(f"  Active modes: {n_active}/{n_filtered}")
print(f"\nPerformance:")
print(f"  Fit range R² (0-1 ps):     {R2_fit:.4f}")
print(f"  Full range R² (0-20 ps):   {R2_spec:.4f}")
print(f"  Frequency domain R²:       {R2_freq:.4f}")
print(f"  C(0) capture:              {100 * C_model_fit[0] / C_total_fit[0]:.1f}%")
print(f"  Spectral integral capture: {100 * J_model_int / J_total_int:.1f}%")

print("\n" + "="*70)
print("Done!")
print("="*70)
