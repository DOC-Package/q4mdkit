#!/usr/bin/env python3
"""
Spectral domain fitting with linear terms only.
Uses NNLS (Non-Negative Least Squares) to fit linear mode contributions
with physical constraint: c_k ≥ 0
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
print("Spectral Fitting: Linear Terms with NNLS")
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

# Filter modes for fitting: 450 ≤ freq ≤ 1800 cm⁻¹
stage1_mask = (mode_freqs >= 450) & (mode_freqs <= 1800)
print(f"\nTotal modes: {len(mode_indices)}")
print(f"Modes to fit (450-1800 cm⁻¹): {np.sum(stage1_mask)}")

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
# Linear terms with NNLS (only 450-1800 cm⁻¹ modes)
# ============================================================================
print("\n" + "="*70)
print("Linear Terms Fitting (NNLS, c_k ≥ 0)")
print(f"Using {np.sum(stage1_mask)} modes with 450 ≤ freq ≤ 1800 cm⁻¹")
print("="*70)

# Use only filtered modes for fitting
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
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Generating Visualization")
print("="*70)

# Create visualization
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Top panel: J_total vs J_linear
ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
ax1.plot(freq_grid, J_linear, 'b-', lw=1.5, label=f'Linear NNLS (R²={R2_linear:.4f})', alpha=0.7)
ax1.fill_between(freq_grid, 0, J_linear, alpha=0.2, color='blue')
ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
ax1.set_title(f'Linear Spectral Density Fitting ({100*J_linear_integral/J_total_integral:.1f}% capture)', 
              fontsize=13, fontweight='bold')
ax1.legend(fontsize=11, loc='upper right')
ax1.grid(True, alpha=0.3)

# Bottom panel: Residual
ax2.plot(freq_grid, J_residual_1, 'r-', lw=1.5, label=f'Residual ({100*J_residual_1_integral/J_total_integral:.1f}%)', alpha=0.7)
ax2.axhline(0, color='k', ls='--', lw=0.5)
ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
ax2.set_ylabel('Residual (cm⁻¹²)', fontsize=12)
ax2.set_title('Fitting Residual', fontsize=11, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectral_fit_linear.png', dpi=150, bbox_inches='tight')
print(f"Saved: spectral_fit_linear.png")
plt.close()

# ============================================================================
# Save results
# ============================================================================
print("\nSaving coefficients...")

# Linear coefficients
with open('spectral_coefficients_linear_nnls.dat', 'w') as f:
    f.write(f"# Linear coefficients from NNLS\n")
    f.write(f"# All c_k ≥ 0 (physical constraint)\n")
    f.write(f"# R² = {R2_linear:.6f}\n")
    f.write(f"# Spectral capture = {100*J_linear_integral/J_total_integral:.2f}%\n")
    f.write(f"# Active modes: {n_active}/{N_modes}\n")
    f.write(f"# Mode  Frequency(cm⁻¹)  c_k  Active\n")
    for i in range(N_modes):
        active = 'Yes' if c_k_physical[i] > 1e-10 else 'No'
        f.write(f"{mode_indices[i]:4d}  {mode_freqs[i]:10.2f}  {c_k_physical[i]:12.6e}  {active}\n")

print(f"Saved: spectral_coefficients_linear_nnls.dat")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\nLinear Fitting Summary:")
print(f"  R² = {R2_linear:.4f}")
print(f"  Spectral capture = {100*J_linear_integral/J_total_integral:.1f}%")
print(f"  Active modes: {n_active}/{N_modes}")

