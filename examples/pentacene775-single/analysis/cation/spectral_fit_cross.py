#!/usr/bin/env python3
"""
Spectral domain fitting using quadratic form optimization:
J(ω) = g^T J_qq(ω) g

Optimization: minimize ||J_total(ω) - g^T J_qq(ω) g||²

J_qq includes:
- Auto-correlations: J_qq,kk(ω) from q_k(t)
- Cross-correlations: J_qq,kl(ω) from q_k(t) and q_l(t) (asymmetric)

Caches computed spectral densities for faster re-runs.
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
print("Spectral Fitting: Quadratic Form Optimization")
print("Minimize ||J_total(ω) - g^T J_qq(ω) g||²")
print("="*70)

# Load data
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff = energy_diff_data[:, 4]

# Load mode information from mode_coords_stats.txt (cation modes)
mode_indices_all = []
mode_freqs_all = []
with open('mode_coords_stats.txt') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            mode_indices_all.append(int(parts[0]))
            mode_freqs_all.append(float(parts[1]))

mode_indices_all = np.array(mode_indices_all)
mode_freqs_all = np.array(mode_freqs_all)

# Filter modes for fitting: 400 ≤ freq ≤ 1800 cm⁻¹
fit_mask = (mode_freqs_all >= 400) & (mode_freqs_all <= 1800)
mode_indices = mode_indices_all[fit_mask]
mode_freqs = mode_freqs_all[fit_mask]
print(f"\nFiltered to {len(mode_indices)} modes with 400 ≤ freq ≤ 1800 cm⁻¹")

# Load mode coordinates from mode_coords_cation.txt
q_k_data = np.loadtxt('mode_coords_cation.txt')
q_k_all = q_k_data[:, 1:]
# Get column indices for selected modes (mode index 7 corresponds to column 0)
mode_col_indices = mode_indices - 7
mode_coords = q_k_all[:, mode_col_indices]

N_frames = min(len(energy_diff), len(mode_coords))
energy_diff = energy_diff[:N_frames]
mode_coords = mode_coords[:N_frames, :]
N_modes = len(mode_freqs)

print(f"  Frames: {N_frames}")
print(f"  Modes: {N_modes}")
print(f"  Timestep: {dt} fs")

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
    """Compute cross-correlation C_xy(τ) = <x(t)y(t+τ)> with Hanning window"""
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

# Frequency grid (420-1800 cm⁻¹)
freq_grid = np.linspace(420, 1800, 1381)
print(f"\nFrequency grid: {len(freq_grid)} points (420-1800 cm⁻¹)")

# ============================================================================
# Load or compute J_qq matrix and J_total
# ============================================================================
cache_file_total = 'spectral_cache_J_total.npy'
cache_file_J_qq = 'spectral_cache_J_qq.npz'

# Check if cache exists and has compatible size
compute_total = True
compute_J_qq = True

if os.path.exists(cache_file_total):
    cached_J_total = np.load(cache_file_total)
    if len(cached_J_total) == len(freq_grid):
        compute_total = False
        J_total = cached_J_total
        print(f"\nLoading J_total from {cache_file_total}...")
    else:
        print(f"\nCache size mismatch for J_total ({len(cached_J_total)} vs {len(freq_grid)}), recomputing...")

if os.path.exists(cache_file_J_qq):
    cache_data = np.load(cache_file_J_qq)
    cached_J_qq = cache_data['J_qq']
    if cached_J_qq.shape[2] == len(freq_grid):
        compute_J_qq = False
        J_qq = cached_J_qq
        print(f"\nLoading J_qq matrix from {cache_file_J_qq}...")
        print(f"  Shape: {J_qq.shape}")
        print(f"  File size: {os.path.getsize(cache_file_J_qq) / 1024**2:.1f} MB")
    else:
        print(f"\nCache size mismatch for J_qq ({cached_J_qq.shape[2]} vs {len(freq_grid)}), recomputing...")

if compute_total:
    print("\nComputing J_total...")
    J_total = compute_spectral_density(energy_diff, freq_grid)
    np.save(cache_file_total, J_total)
    print(f"  Saved to {cache_file_total}")

J_total_integral = np.trapezoid(J_total, freq_grid)
print(f"  ∫J_total = {J_total_integral:.6e} cm⁻¹²")

if compute_J_qq:
    print(f"\nComputing J_qq matrix ({N_modes} × {N_modes} × {len(freq_grid)})...")
    print(f"  This may take several minutes...")
    
    J_qq = np.zeros((N_modes, N_modes, len(freq_grid)))
    
    # Compute all elements of J_qq matrix
    for k in range(N_modes):
        # Auto-correlation: J_qq[k,k] from q_k
        J_qq[k, k, :] = compute_spectral_density(mode_coords[:, k], freq_grid)
        
        # Cross-correlations: J_qq[k,l] from q_k and q_l (k ≠ l)
        for l in range(N_modes):
            if k != l:
                J_qq[k, l, :] = compute_cross_spectral_density(mode_coords[:, k], mode_coords[:, l], freq_grid)
        
        if (k + 1) % 5 == 0:
            print(f"    Progress: {k+1}/{N_modes} modes")
    
    # Save to file
    np.savez_compressed(cache_file_J_qq, J_qq=J_qq)
    print(f"  Saved to {cache_file_J_qq}")
    print(f"  File size: {os.path.getsize(cache_file_J_qq) / 1024**2:.1f} MB")
else:
    print(f"\nLoading J_qq matrix from {cache_file_J_qq}...")
    cache_data = np.load(cache_file_J_qq)
    J_qq = cache_data['J_qq']
    print(f"  Shape: {J_qq.shape}")
    print(f"  File size: {os.path.getsize(cache_file_J_qq) / 1024**2:.1f} MB")

# Analyze J_qq matrix
diag_mean = np.mean(np.abs(J_qq[np.arange(N_modes), np.arange(N_modes), :]))
mask = ~np.eye(N_modes, dtype=bool)
off_diag_indices = np.where(mask)
offdiag_mean = np.mean(np.abs(J_qq[off_diag_indices[0], off_diag_indices[1], :]))

print(f"\nJ_qq matrix analysis:")
print(f"  Mean |J_qq,kk| (auto): {diag_mean:.6e}")
print(f"  Mean |J_qq,kl| (cross, k≠l): {offdiag_mean:.6e}")
print(f"  Ratio (cross/auto): {offdiag_mean/diag_mean:.4f}")

# ============================================================================
# Normalization for numerical stability
# ============================================================================
print("\n" + "="*70)
print("Applying Normalization")
print("="*70)

# Use J_total RMS as the common scale for both J_total and J_qq
J_total_rms = np.sqrt(np.mean(J_total**2))
J_qq_rms = np.sqrt(np.mean(J_qq**2))

print(f"\nOriginal scales (RMS):")
print(f"  J_total RMS: {J_total_rms:.6e}")
print(f"  J_qq RMS: {J_qq_rms:.6e}")
print(f"  Ratio (J_total/J_qq): {J_total_rms / J_qq_rms:.6e}")

# Normalize both by J_total_rms to keep them on the same scale
scale_factor = J_total_rms
J_total_normalized = J_total / scale_factor
J_qq_normalized = J_qq / scale_factor

print(f"\nNormalization:")
print(f"  Common scale factor: {scale_factor:.6e}")
print(f"  J_total_normalized RMS: {np.sqrt(np.mean(J_total_normalized**2)):.6f}")
print(f"  J_qq_normalized RMS: {np.sqrt(np.mean(J_qq_normalized**2)):.6f}")

# ============================================================================
# Optimization: minimize ||J_total(ω) - g^T J_qq(ω) g||² + λ||g||²
# ============================================================================
print("\n" + "="*70)
print("Optimizing g to minimize ||J_total - g^T J_qq g||² + λ||g||²")
print("="*70)

# Regularization parameter
lambda_reg = 1e-6

def objective(g):
    """Objective function with L2 regularization"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq_normalized, g)
    residual = J_total_normalized - J_model
    data_term = np.sum(residual**2)
    reg_term = lambda_reg * np.sum(g**2)
    return data_term + reg_term

def gradient(g):
    """Gradient with L2 regularization"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq_normalized, g)
    residual = J_total_normalized - J_model
    # Data term gradient
    grad_data = -4 * np.einsum('k,ijk,j->i', residual, J_qq_normalized, g)
    # Regularization gradient
    grad_reg = 2 * lambda_reg * g
    return grad_data + grad_reg

# Better initial guess: use linear approximation
# Assume g ≈ sqrt(c_k) where c_k would be from linear fit
g_init = np.ones(N_modes) * 0.1

print(f"\nOptimization settings:")
print(f"  Method: L-BFGS-B")
print(f"  Regularization: λ = {lambda_reg:.0e}")
print(f"  Bounds: -100 ≤ g_k ≤ 100")
print(f"  Using normalized J_total and J_qq")
print(f"  Initial objective: {objective(g_init):.6e}")

print("\nOptimizing...")
bounds = [(-100, 100) for _ in range(N_modes)]
result = minimize(
    objective,
    g_init,
    method='L-BFGS-B',
    jac=gradient,
    bounds=bounds,
    options={'maxiter': 500, 'disp': True, 'ftol': 1e-10, 'gtol': 1e-8}
)

g_opt = result.x

print(f"\nOptimization result:")
print(f"  Success: {result.success}")
print(f"  Message: {result.message}")
print(f"  Iterations: {result.nit}")
print(f"  Final objective: {result.fun:.6e}")
print(f"  Mean |g|: {np.mean(np.abs(g_opt)):.6e}")
print(f"  Max |g|: {np.max(np.abs(g_opt)):.6e}")

# Compute final J_model (in physical units)
J_model_normalized = np.einsum('i,ijk,j->k', g_opt, J_qq_normalized, g_opt)
# Convert back using the common scale factor
J_model = J_model_normalized * scale_factor

# Decompose into diagonal and off-diagonal contributions
# Diagonal: Σ_k g_k^2 J_qq[k,k,ω]
J_diag_normalized = np.einsum('i,i,iik->k', g_opt, g_opt, J_qq_normalized)
J_diag = J_diag_normalized * scale_factor

# Off-diagonal: Σ_k Σ_{l≠k} g_k g_l J_qq[k,l,ω]
J_offdiag = J_model - J_diag

# Integrals
J_diag_integral = np.trapezoid(J_diag, freq_grid)
J_offdiag_integral = np.trapezoid(J_offdiag, freq_grid)

print(f"\nContribution decomposition:")
print(f"  Diagonal (auto-correlation): {100*J_diag_integral/J_total_integral:.1f}%")
print(f"  Off-diagonal (cross-correlation): {100*J_offdiag_integral/J_total_integral:.1f}%")
print(f"  Total: {100*(J_diag_integral+J_offdiag_integral)/J_total_integral:.1f}%")

# Check for negative values
neg_mask = J_model < 0
n_neg = np.sum(neg_mask)
print(f"\nPhysical constraint check:")
print(f"  Negative J points: {n_neg}/{len(freq_grid)} ({100*n_neg/len(freq_grid):.1f}%)")
if n_neg > 0:
    print(f"  Min J: {J_model.min():.6e}")
    print(f"  WARNING: Spectral density is negative at some frequencies!")

# ============================================================================
# Fit Quality
# ============================================================================
print("\n" + "="*70)
print("Fit Quality")
print("="*70)

# R² and capture
residual = J_total - J_model
ss_res = np.sum(residual**2)
ss_tot = np.sum((J_total - np.mean(J_total))**2)
R2 = 1 - ss_res / ss_tot

J_model_integral = np.trapezoid(J_model, freq_grid)
capture = 100 * J_model_integral / J_total_integral

print(f"\n  R² = {R2:.4f}")
print(f"  ∫J_total = {J_total_integral:.6e} cm⁻¹²")
print(f"  ∫J_model = {J_model_integral:.6e} cm⁻¹²")
print(f"  Capture = {capture:.1f}%")

# Count active coefficients
active = np.abs(g_opt) > 1e-6
n_active = np.sum(active)
n_positive = np.sum(g_opt > 1e-6)
n_negative = np.sum(g_opt < -1e-6)

print(f"\n  Active coefficients: {n_active}/{N_modes}")
print(f"    Positive: {n_positive}")
print(f"    Negative: {n_negative}")

# Top modes by |g_k|
top_idx = np.argsort(np.abs(g_opt))[-15:][::-1]
print(f"\nTop 15 modes by |g_k|:")
for rank, idx in enumerate(top_idx, 1):
    sign = '+' if g_opt[idx] > 0 else '-'
    print(f"  {rank:2d}. Mode {mode_indices[idx]:2d} ({mode_freqs[idx]:6.1f} cm⁻¹): "
          f"{sign} {abs(g_opt[idx]):.6e}")



# ============================================================================
# Save results
# ============================================================================
print("\nSaving coefficients...")

# Coefficients
with open('g_k_quadratic_form.dat', 'w') as f:
    f.write(f"# Coefficients from Quadratic Form fit\n")
    f.write(f"# J(ω) = g^T J_qq(ω) g\n")
    f.write(f"# Fitted with normalized J_total and J_qq\n")
    f.write(f"# Normalization:\n")
    f.write(f"#   Common scale factor = {scale_factor:.6e}\n")
    f.write(f"#   (Both J_total and J_qq normalized by J_total_rms)\n")
    f.write(f"# R² = {R2:.6f}\n")
    f.write(f"# Capture = {capture:.2f}%\n")
    f.write(f"# Fitting range: {freq_grid[0]:.0f}-{freq_grid[-1]:.0f} cm⁻¹\n")
    f.write(f"# Mode range: 400-1800 cm⁻¹\n")
    f.write(f"#\n")
    f.write("# Mode  Frequency(cm⁻¹)  g_k  |g_k|  Sign\n")
    for i in range(N_modes):
        sign = '+' if g_opt[i] > 0 else '-'
        f.write(f"{mode_indices[i]:4d}  {mode_freqs[i]:10.2f}  {g_opt[i]:12.6e}  {abs(g_opt[i]):12.6e}  {sign}\n")

print(f"Saved: g_k_quadratic_form.dat")

# Spectral densities
with open('spectral_density_quadratic_form.dat', 'w') as f:
    f.write("# Frequency(cm⁻¹)  J_total(cm⁻¹²)  J_model(cm⁻¹²)  J_diag(cm⁻¹²)  J_offdiag(cm⁻¹²)  Residual(cm⁻¹²)\n")
    for nu, Jt, Jm, Jd, Jo, res in zip(freq_grid, J_total, J_model, J_diag, J_offdiag, residual):
        f.write(f"{nu:10.4f}  {Jt:14.6e}  {Jm:14.6e}  {Jd:14.6e}  {Jo:14.6e}  {res:14.6e}\n")

print(f"Saved: spectral_density_quadratic_form.dat")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\nQuadratic Form Optimization Summary:")
print(f"  Method: minimize ||J_total(ω) - g^T J_qq(ω) g||²")
print(f"  Normalization: Common scale factor for numerical stability")
print(f"    Scale factor = {scale_factor:.6e}")
print(f"    (Both J_total and J_qq normalized by J_total_rms)")
print(f"  R² = {R2:.4f}")
print(f"  Spectral capture = {capture:.1f}%")
print(f"  Active modes: {n_active}/{N_modes}")
print(f"    Positive: {n_positive}")
print(f"    Negative: {n_negative}")
print(f"  Negative J points: {n_neg}/{len(freq_grid)} ({100*n_neg/len(freq_grid):.1f}%)")
print(f"\n  Cache files created:")
print(f"    - {cache_file_total}")
print(f"    - {cache_file_J_qq}")
print(f"\n  Advantages:")
print(f"    ✓ Naturally includes mode coupling via J_qq[k,l]")
print(f"    ✓ Theoretically exact for linear model e(t) = g^T q(t)")
print(f"    ✓ No ad-hoc separation into stages")
print(f"    ✓ Normalized with common scale for consistent units")
