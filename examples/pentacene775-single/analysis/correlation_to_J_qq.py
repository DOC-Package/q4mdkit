#!/usr/bin/env python3
"""
Compute J_qq matrix from pre-computed correlation functions.
This is much faster than compute_J_qq.py when experimenting with different
frequency grids, since correlation functions are already computed.
"""
import numpy as np
from scipy.integrate import simpson

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs

print("="*70)
print("Computing J_qq Matrix from Cached Correlation Functions")
print("="*70)

# Load correlation functions
print("\nLoading correlation functions...")
try:
    acf_matrix = np.load('acf_matrix.npy')
    ccf_matrix = np.load('ccf_matrix.npy')
    metadata = np.load('correlation_metadata.npz')
    
    N_modes = metadata['N_modes']
    max_lag = metadata['max_lag']
    mode_indices = metadata['mode_indices']
    mode_freqs = metadata['mode_freqs']
    
    print(f"  ACF matrix shape: {acf_matrix.shape}")
    print(f"  CCF matrix shape: {ccf_matrix.shape}")
    print(f"  N_modes: {N_modes}")
    print(f"  max_lag: {max_lag}")
    
except FileNotFoundError as e:
    print(f"  ERROR: {e}")
    print("  Please run compute_correlation_functions.py first.")
    exit(1)

# Frequency grid
freq_min, freq_max = 420, 1800
n_freqs = 1381
freq_grid = np.linspace(freq_min, freq_max, n_freqs)

print(f"\nFrequency grid: {freq_min}-{freq_max} cm⁻¹, {n_freqs} points")

def correlation_to_spectral_density(corr, dt, freq_grid, max_lag=5000):
    """
    Convert correlation function to spectral density.
    
    J(ν) = (2πcν/kBT) × ∫ C(τ) cos(2πcντ) dτ
    
    Use shorter integration time to avoid artificial peak sharpening.
    """
    # Reduce max_lag to avoid numerical artifacts from long integration
    max_lag_reduced = min(max_lag, 5000)  # 20ps at dt=4fs
    tau = np.arange(max_lag_reduced) * dt
    window = np.hanning(max_lag_reduced)
    corr_windowed = corr[:max_lag_reduced] #* window
    
    J = np.zeros(len(freq_grid))
    for i, nu in enumerate(freq_grid):
        omega = 2 * np.pi * c_cm_fs * nu
        integrand = corr_windowed * np.cos(omega * tau)
        integral = simpson(integrand, x=tau)
        #prefactor = 2 * np.pi * c_cm_fs * omega / (kB_cm * T)
        prefactor = 2 * np.pi * c_cm_fs
        J[i] = prefactor * integral
    
    return J

# Compute J_qq matrix
print(f"\nComputing J_qq matrix ({N_modes}×{N_modes}×{n_freqs})...")
J_qq = np.zeros((N_modes, N_modes, n_freqs))

for k in range(N_modes):
    if (k+1) % 10 == 0:
        print(f"  Progress: {k+1}/{N_modes}")
    
    # Diagonal
    J_qq[k, k, :] = correlation_to_spectral_density(
        acf_matrix[k, :], dt, freq_grid, max_lag
    )
    
    # Off-diagonal
    for l in range(k+1, N_modes):
        J_kl = correlation_to_spectral_density(
            ccf_matrix[k, l, :], dt, freq_grid, max_lag
        )
        J_qq[k, l, :] = J_kl
        J_qq[l, k, :] = J_kl

print("  Done.")

# Save J_qq matrix
print(f"\nSaving J_qq_matrix.npy...")
np.save('J_qq_matrix.npy', J_qq)
print(f"  Saved: {J_qq.nbytes / 1024**2:.1f} MB")

# Statistics
diag_mean = np.mean(np.abs(J_qq[np.arange(N_modes), np.arange(N_modes), :]))
mask = ~np.eye(N_modes, dtype=bool)
off_diag_indices = np.where(mask)
offdiag_mean = np.mean(np.abs(J_qq[off_diag_indices[0], off_diag_indices[1], :]))

print(f"\nCross-spectrum analysis:")
print(f"  Mean |J_qq,kk|: {diag_mean:.6e}")
print(f"  Mean |J_qq,kl| (k≠l): {offdiag_mean:.6e}")
print(f"  Ratio: {offdiag_mean/diag_mean:.4f}")

print("\n" + "="*70)
print("COMPLETE")
print("="*70)
print("\nYou can now run spectral_fit_quadratic_form_cached.py")
