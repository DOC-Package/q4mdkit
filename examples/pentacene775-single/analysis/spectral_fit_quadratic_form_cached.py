#!/usr/bin/env python3
"""
Quadratic form fitting using pre-computed J_qq matrix:
J(ω) = g^T J_qq(ω) g

Loads J_qq_matrix.npy for fast execution.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs

print("="*70)
print("Quadratic Form Fitting (Using Cached J_qq Matrix)")
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

# Filter modes: 400 ≤ freq ≤ 1800 cm⁻¹ (must match cached J_qq)
freq_mask = (mode_freqs >= 400) & (mode_freqs <= 1800)
mode_indices = mode_indices[freq_mask]
mode_freqs = mode_freqs[freq_mask]
print(f"\nFiltered to {len(mode_indices)} modes with 400 ≤ freq ≤ 1800 cm⁻¹")

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
print(f"  Timestep: {dt} fs")

# ============================================================================
# Load cached J_qq matrix and compute J_total
# ============================================================================
print("\n" + "="*70)
print("Loading Cached Spectral Densities")
print("="*70)

def compute_spectral_density(x, freq_grid, max_lag=5000):
    """Compute spectral density J(ν) from time series"""
    x_mean = np.mean(x)
    x_centered = x - x_mean
    N = len(x)
    
    # Auto-correlation function
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(x_centered[:N-lag] * x_centered[lag:])
    
    # Hanning window
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    acf = acf #* window
    
    # Fourier transform
    tau_array = np.arange(max_lag) * dt
    J = np.zeros(len(freq_grid))
    for i, nu in enumerate(freq_grid):
        omega = 2 * np.pi * c_cm_fs * nu  # Angular frequency in rad/fs
        integrand = acf * np.cos(omega * tau_array)
        J[i] = simpson(integrand, x=tau_array)
    
    # Apply prefactor: (2πcν / kBT)
    #factor = 2 * np.pi * c_cm_fs * freq_grid / (kB_cm * T)
    factor = 2 * np.pi * c_cm_fs
    J = J * factor
    
    return J

# Frequency grid for fitting
freq_min, freq_max = 420, 1800
freq_grid = np.linspace(freq_min, freq_max, 1381)
n_freqs = len(freq_grid)

print(f"\nFrequency grid: {freq_min}-{freq_max} cm⁻¹, {n_freqs} points")

# Compute J_total(ω)
print("\nComputing J_total(ω)...")
J_total = compute_spectral_density(energy_diff, freq_grid)
integral_total = simpson(J_total, x=freq_grid)
print(f"  ∫J_total = {integral_total:.6e} cm⁻¹²")

# Load cached J_qq matrix
print("\nLoading J_qq matrix from cache...")
try:
    J_qq = np.load('J_qq_matrix.npy')
    print(f"  Loaded from J_qq_matrix.npy")
    print(f"  Shape: {J_qq.shape}")
    print(f"  Size: {J_qq.nbytes / 1024**2:.1f} MB")
    
    if J_qq.shape != (N_modes, N_modes, n_freqs):
        raise ValueError(f"Shape mismatch! Expected ({N_modes}, {N_modes}, {n_freqs}), got {J_qq.shape}")
    
except FileNotFoundError:
    print("  ERROR: J_qq_matrix.npy not found!")
    print("  Please run compute_J_qq.py first.")
    exit(1)
except Exception as e:
    print(f"  ERROR: {e}")
    exit(1)

# Check cross-spectrum magnitude
diag_mean = np.mean(np.abs(J_qq[np.arange(N_modes), np.arange(N_modes), :]))
mask = ~np.eye(N_modes, dtype=bool)
off_diag_indices = np.where(mask)
offdiag_mean = np.mean(np.abs(J_qq[off_diag_indices[0], off_diag_indices[1], :]))

print(f"\nCross-spectrum analysis:")
print(f"  Mean |J_qq,kk|: {diag_mean:.6e}")
print(f"  Mean |J_qq,kl| (k≠l): {offdiag_mean:.6e}")
print(f"  Ratio: {offdiag_mean/diag_mean:.4f}")

# ============================================================================
# Optimization: minimize ||J_total(ω) - g^T J_qq(ω) g||²
# ============================================================================
print("\n" + "="*70)
print("Optimizing g to minimize ||J_total - g^T J_qq g||²")
print("="*70)

def objective(g):
    """Objective function: sum of squared residuals"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq, g)
    residual = J_total - J_model
    return np.sum(residual**2)

def gradient(g):
    """Gradient of objective function"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq, g)
    residual = J_total - J_model
    # ∂/∂g_k = -2 Σ_ω residual(ω) * 2 Σ_l g_l J_qq[k,l,ω]
    grad = -4 * np.einsum('k,ijk,j->i', residual, J_qq, g)
    return grad

# Initial guess: use traditional linear fit results if available
g_init = np.ones(N_modes) * 0.001
print("\nUsing uniform initial guess")

print(f"\nOptimization settings:")
print(f"  Method: L-BFGS-B")
print(f"  Bounds: -1000 ≤ g_k ≤ 1000")
print(f"  Initial objective: {objective(g_init):.6e}")

print("\nOptimizing...")
bounds = [(-1000, 1000) for _ in range(N_modes)]
result = minimize(
    objective,
    g_init,
    method='L-BFGS-B',
    jac=gradient,
    bounds=bounds,
    options={'maxiter': 100, 'disp': True}
)

g_opt = result.x

print(f"\nOptimization result:")
print(f"  Success: {result.success}")
print(f"  Message: {result.message}")
print(f"  Iterations: {result.nit}")
print(f"  Final objective: {result.fun:.6e}")
print(f"  Mean |g|: {np.mean(np.abs(g_opt)):.6e}")
print(f"  Max |g|: {np.max(np.abs(g_opt)):.6e}")

# Compute final J_model
J_model = np.einsum('i,ijk,j->k', g_opt, J_qq, g_opt)

# Check for negative values
neg_mask = J_model < 0
n_neg = np.sum(neg_mask)
print(f"\nPhysical constraint check:")
print(f"  Negative J points: {n_neg}/{n_freqs} ({100*n_neg/n_freqs:.1f}%)")
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
r_squared = 1 - ss_res / ss_tot

integral_model = simpson(J_model, x=freq_grid)
capture = 100 * integral_model / integral_total

print(f"\n  R² = {r_squared:.4f}")
print(f"  ∫J_total = {integral_total:.6e} cm⁻¹²")
print(f"  ∫J_model = {integral_model:.6e} cm⁻¹²")
print(f"  Capture = {capture:.1f}%")

# Count active coefficients
active = np.abs(g_opt) > 1e-6
n_active = np.sum(active)
n_positive = np.sum(g_opt > 1e-6)
n_negative = np.sum(g_opt < -1e-6)

print(f"\n  Active coefficients: {n_active}/{N_modes}")
print(f"    Positive: {n_positive}")
print(f"    Negative: {n_negative}")

# ============================================================================
# Visualization
# ============================================================================
print("\n" + "="*70)
print("Generating Visualization")
print("="*70)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot 1: J(ω) comparison
ax1.plot(freq_grid, J_total, 'k-', linewidth=1.5, label='$J_{\\rm total}(\\omega)$ (data)', alpha=0.7)
ax1.plot(freq_grid, J_model, 'r-', linewidth=1.0, label='$J_{\\rm model}(\\omega)$ (fit)')
ax1.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax1.set_xlabel('Frequency (cm$^{-1}$)', fontsize=11)
ax1.set_ylabel('$J(\\omega)$ (cm$^{-1}$)', fontsize=11)
ax1.set_xlim(freq_min, freq_max)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_title(f'Quadratic Form Fit: $R^2$ = {r_squared:.4f}, Capture = {capture:.1f}%', fontsize=12)

# Plot 2: Residual
ax2.plot(freq_grid, residual, 'b-', linewidth=0.8)
ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=11)
ax2.set_ylabel('Residual (cm$^{-1}$)', fontsize=11)
ax2.set_xlim(freq_min, freq_max)
ax2.grid(True, alpha=0.3)
ax2.set_title('Residual: $J_{\\rm total} - J_{\\rm model}$', fontsize=12)

plt.tight_layout()
plt.savefig('quadratic_form_fit_cached.png', dpi=150, bbox_inches='tight')
print("Saved: quadratic_form_fit_cached.png")

# Save coefficients
print("\nSaving coefficients...")
with open('g_k_quadratic_form_cached.dat', 'w') as f:
    f.write("# Coefficients from Quadratic Form fit (using cached J_qq)\n")
    f.write(f"# R² = {r_squared:.6f}\n")
    f.write(f"# Capture = {capture:.2f}%\n")
    f.write(f"# Fitting range: {freq_min}-{freq_max} cm⁻¹\n")
    f.write(f"# Mode range: 400-1800 cm⁻¹\n")
    f.write("#\n")
    f.write("# Mode_Index  Frequency(cm⁻¹)  g_k\n")
    for idx, freq, g in zip(mode_indices, mode_freqs, g_opt):
        f.write(f"{idx:4d}  {freq:10.4f}  {g:12.6e}\n")

print("Saved: g_k_quadratic_form_cached.dat")

# Save spectral densities
print("\nSaving spectral densities...")
with open('spectral_density_quadratic_form_cached.dat', 'w') as f:
    f.write("# Frequency(cm⁻¹)  J_total(cm⁻¹)  J_model(cm⁻¹)  Residual(cm⁻¹)\n")
    for nu, Jt, Jm, res in zip(freq_grid, J_total, J_model, residual):
        f.write(f"{nu:10.4f}  {Jt:14.6e}  {Jm:14.6e}  {res:14.6e}\n")

print("Saved: spectral_density_quadratic_form_cached.dat")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)

print(f"\nQuadratic Form Method Summary:")
print(f"  Fitting range: {freq_min}-{freq_max} cm⁻¹")
print(f"  Mode range: 400-1800 cm⁻¹")
print(f"  R² = {r_squared:.4f}")
print(f"  Spectral capture = {capture:.1f}%")
print(f"  Active modes: {n_active}/{N_modes}")
print(f"  Negative J points: {n_neg}/{n_freqs} ({100*n_neg/n_freqs:.1f}%)")

print("\nAdvantages:")
print("  ✓ Naturally includes mode coupling via J_qq,kl")
print("  ✓ Theoretically exact for linear model e(t) = g^T q(t)")
print("  ✓ Fast execution using cached J_qq matrix")
