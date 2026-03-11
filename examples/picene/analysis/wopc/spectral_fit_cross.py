#!/usr/bin/env python3
"""
Spectral domain fitting using quadratic form optimization (Internal only):
J_internal(ω) = g^T J_qq(ω) g

Optimization: minimize ||J_internal(ω) - g^T J_qq(ω) g||²

J_qq includes:
- Auto-correlations: J_qq,kk(ω) from q_k(t)
- Cross-correlations: J_qq,kl(ω) from q_k(t) and q_l(t) (asymmetric)

Caches computed spectral densities for faster re-runs.
Mode range: 500-1900 cm⁻¹
"""
import numpy as np
from scipy.integrate import simpson
from scipy.fft import fft, ifft
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os.path

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
two_pi_c = 2 * np.pi * c_cm_fs  # rad·cm/fs
T = 300.0  # K
dt = 4.0  # fs

# Segmentation parameters (matching plot_local_v2.py)
segment_ps = 100.0  # ps
corr_ps = 8.0  # ps

print("="*70)
print("Spectral Fitting: Quadratic Form Optimization (Internal only)")
print("Minimize ||J_internal(ω) - g^T J_qq(ω) g||²")
print("="*70)

# Unit conversion
cm_to_meV = 0.123984  # meV/cm⁻¹

# Load J_internal from spectral_density_comparison.dat
print("\nLoading J_internal from spectral_density_comparison.dat...")
spectral_data = np.loadtxt('spectral_density.dat')
freq_data = spectral_data[:, 0]  # nu (cm⁻¹)
J_internal_meV = spectral_data[:, 1]  # J_internal (meV)
J_internal_data = J_internal_meV / cm_to_meV  # Convert back to cm⁻¹

print(f"  Loaded {len(freq_data)} frequency points")
print(f"  Frequency range: {freq_data[0]:.1f} - {freq_data[-1]:.1f} cm⁻¹")

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

# Filter modes for fitting: 500 ≤ freq ≤ 1900 cm⁻¹
fit_mask = (mode_freqs_all >= 500) & (mode_freqs_all <= 1900)
mode_indices = mode_indices_all[fit_mask]
mode_freqs = mode_freqs_all[fit_mask]
print(f"\nFiltered to {len(mode_indices)} modes with 500 ≤ freq ≤ 1900 cm⁻¹")

# Load mode coordinates from mode_coords_normalized.txt
q_k_data = np.loadtxt('mode_coords_normalized.txt')
q_k_all = q_k_data[:, 1:]
# Get column indices for selected modes (mode index 7 corresponds to column 0)
mode_col_indices = mode_indices - 7
mode_coords = q_k_all[:, mode_col_indices]

N_frames = len(mode_coords)
N_modes = len(mode_freqs)

# Compute segment and correlation lengths
segment_length = int(segment_ps * 1000 / dt)  # frames
corr_length = int(corr_ps * 1000 / dt)  # frames

print(f"  Frames: {N_frames}")
print(f"  Modes: {N_modes}")
print(f"  Timestep: {dt} fs")
print(f"  Segment length: {segment_ps} ps ({segment_length} frames)")
print(f"  Correlation length: {corr_ps} ps ({corr_length} frames)")

# =============================================================================
# Correlation function computation (matching plot_local_v2.py)
# =============================================================================
def compute_correlation_fft(delta: np.ndarray, dt: float) -> tuple:
    """Compute autocorrelation function using FFT."""
    N = len(delta)
    delta_padded = np.concatenate([delta, np.zeros(N)])
    fft_delta = fft(delta_padded)
    power_spectrum = np.abs(fft_delta)**2
    autocorr_full = np.real(ifft(power_spectrum))[:N]
    norm = np.arange(N, 0, -1)
    C = autocorr_full / norm
    t_corr = np.arange(N) * dt
    return t_corr, C


def compute_crosscorr_fft(x: np.ndarray, y: np.ndarray, dt: float) -> tuple:
    """Compute cross-correlation function using FFT.
    C_xy(τ) = <x(t) y(t+τ)>
    """
    N = len(x)
    x_padded = np.concatenate([x, np.zeros(N)])
    y_padded = np.concatenate([y, np.zeros(N)])
    fft_x = fft(x_padded)
    fft_y = fft(y_padded)
    cross_spectrum = np.conj(fft_x) * fft_y
    crosscorr_full = np.real(ifft(cross_spectrum))[:N]
    norm = np.arange(N, 0, -1)
    C = crosscorr_full / norm
    t_corr = np.arange(N) * dt
    return t_corr, C


def compute_correlation_segmented(delta: np.ndarray, dt: float, 
                                   segment_length: int, corr_length: int) -> tuple:
    """Compute autocorrelation function with segment averaging."""
    N = len(delta)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        segment = delta[start:end]
        segment = segment - np.mean(segment)
        _, C_seg = compute_correlation_fft(segment, dt)
        C_sum += C_seg[:corr_length]
    
    C_avg = C_sum / n_segments
    t_corr = np.arange(corr_length) * dt
    
    return t_corr, C_avg, n_segments


def compute_crosscorr_segmented(x: np.ndarray, y: np.ndarray, dt: float,
                                 segment_length: int, corr_length: int) -> tuple:
    """Compute cross-correlation function with segment averaging.
    C_xy(τ) = <x(t) y(t+τ)>
    """
    N = len(x)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        seg_x = x[start:end] - np.mean(x[start:end])
        seg_y = y[start:end] - np.mean(y[start:end])
        _, C_seg = compute_crosscorr_fft(seg_x, seg_y, dt)
        C_sum += C_seg[:corr_length]
    
    C_avg = C_sum / n_segments
    t_corr = np.arange(corr_length) * dt
    
    return t_corr, C_avg, n_segments


# =============================================================================
# Spectral density computation (matching plot_local_v2.py)
# =============================================================================
def correlation_to_spectral_density(t_corr: np.ndarray, C: np.ndarray, 
                                     nu_out: np.ndarray, T: float,
                                     use_window: bool = False) -> np.ndarray:
    """
    Compute spectral density from correlation function.
    J(ν̃) = (2πc ν̃ / 2k_B T) ∫_0^∞ dt C(t) cos(2πc ν̃ t)
    """
    dt = t_corr[1] - t_corr[0]
    
    if use_window:
        window = 0.5 * (1 + np.cos(np.pi * t_corr / t_corr[-1]))
        C_windowed = C * window
    else:
        C_windowed = C
    
    J = np.zeros(len(nu_out))
    
    for i, nu in enumerate(nu_out):
        omega = two_pi_c * nu
        integrand = C_windowed * np.cos(omega * t_corr)
        integral = simpson(integrand, dx=dt)
        J[i] = two_pi_c * nu / (2.0 * kB_cm * T) * integral
    
    return J


def compute_spectral_density(x, freq):
    """Compute spectral density using segmented correlation (matching plot_local_v2.py)."""
    t_corr, C, _ = compute_correlation_segmented(x, dt, segment_length, corr_length)
    return correlation_to_spectral_density(t_corr, C, freq, T)


def compute_cross_spectral_density(x, y, freq):
    """Compute cross-spectral density using segmented correlation (matching plot_local_v2.py)."""
    t_corr, C, _ = compute_crosscorr_segmented(x, y, dt, segment_length, corr_length)
    return correlation_to_spectral_density(t_corr, C, freq, T)

# Frequency grid (500-1900 cm⁻¹)
freq_grid = np.linspace(500, 1900, 701)
print(f"\nFrequency grid: {len(freq_grid)} points (500-1900 cm⁻¹)")

# Interpolate J_internal to the frequency grid
from scipy.interpolate import interp1d
J_internal_interp = interp1d(freq_data, J_internal_data, kind='linear', fill_value='extrapolate')
J_target = J_internal_interp(freq_grid)
print(f"  Interpolated J_internal to fitting grid")

# ============================================================================
# Load or compute J_qq matrix (J_internal already loaded above)
# ============================================================================
# Cache file name includes segmentation parameters to invalidate old caches
cache_file_J_qq = f'spectral_cache_J_qq_seg{int(segment_ps)}ps_corr{int(corr_ps)}ps.npz'

# Check if cache exists and has compatible size
compute_J_qq = True

J_target_integral = np.trapezoid(J_target, freq_grid)
print(f"  ∫J_internal = {J_target_integral:.6e} cm⁻¹²")

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

# Use J_target RMS as the common scale for both J_target and J_qq
J_target_rms = np.sqrt(np.mean(J_target**2))
J_qq_rms = np.sqrt(np.mean(J_qq**2))

print(f"\nOriginal scales (RMS):")
print(f"  J_internal RMS: {J_target_rms:.6e}")
print(f"  J_qq RMS: {J_qq_rms:.6e}")
print(f"  Ratio (J_internal/J_qq): {J_target_rms / J_qq_rms:.6e}")

# Normalize both by J_target_rms to keep them on the same scale
scale_factor = J_target_rms
J_target_normalized = J_target / scale_factor
J_qq_normalized = J_qq / scale_factor

print(f"\nNormalization:")
print(f"  Common scale factor: {scale_factor:.6e}")
print(f"  J_internal_normalized RMS: {np.sqrt(np.mean(J_target_normalized**2)):.6f}")
print(f"  J_qq_normalized RMS: {np.sqrt(np.mean(J_qq_normalized**2)):.6f}")

# ============================================================================
# Optimization: minimize ||J_internal(ω) - κ^T J_qq(ω) κ||² + λ||κ||²
# ============================================================================
print("\n" + "="*70)
print("Optimizing κ to minimize ||J_internal - κ^T J_qq κ||² + λ||κ||²")
print("="*70)

# Regularization parameter
lambda_reg = 1e-6

def objective(g):
    """Objective function with L2 regularization"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq_normalized, g)
    residual = J_target_normalized - J_model
    data_term = np.sum(residual**2)
    reg_term = lambda_reg * np.sum(g**2)
    return data_term + reg_term

def gradient(g):
    """Gradient with L2 regularization"""
    J_model = np.einsum('i,ijk,j->k', g, J_qq_normalized, g)
    residual = J_target_normalized - J_model
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
print(f"  Bounds: -1e5 ≤ κ_k ≤ 1e5")
print(f"  Using normalized J_internal and J_qq")
print(f"  Initial objective: {objective(g_init):.6e}")

print("\nOptimizing...")
bounds = [(-1e5, 1e5) for _ in range(N_modes)]
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
print(f"  Mean |κ|: {np.mean(np.abs(g_opt)):.6e}")
print(f"  Max |κ|: {np.max(np.abs(g_opt)):.6e}")

# Compute final J_model (in physical units)
J_model_normalized = np.einsum('i,ijk,j->k', g_opt, J_qq_normalized, g_opt)
# Convert back using the common scale factor
J_model = J_model_normalized * scale_factor

# Decompose into diagonal and off-diagonal contributions
# Diagonal: Σ_k κ_k² J_qq[k,k,ω]
J_diag_normalized = np.einsum('i,i,iik->k', g_opt, g_opt, J_qq_normalized)
J_diag = J_diag_normalized * scale_factor

# Off-diagonal: Σ_k Σ_{l≠k} κ_k κ_l J_qq[k,l,ω]
J_offdiag = J_model - J_diag

# Integrals
J_diag_integral = np.trapezoid(J_diag, freq_grid)
J_offdiag_integral = np.trapezoid(J_offdiag, freq_grid)

print(f"\nContribution decomposition:")
print(f"  Diagonal (auto-correlation): {100*J_diag_integral/J_target_integral:.1f}%")
print(f"  Off-diagonal (cross-correlation): {100*J_offdiag_integral/J_target_integral:.1f}%")
print(f"  Total: {100*(J_diag_integral+J_offdiag_integral)/J_target_integral:.1f}%")

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
residual = J_target - J_model
ss_res = np.sum(residual**2)
ss_tot = np.sum((J_target - np.mean(J_target))**2)
R2 = 1 - ss_res / ss_tot

J_model_integral = np.trapezoid(J_model, freq_grid)
capture = 100 * J_model_integral / J_target_integral

print(f"\n  R² = {R2:.4f}")
print(f"  ∫J_internal = {J_target_integral:.6e} cm⁻¹²")
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

# Top modes by |κ_k|
top_idx = np.argsort(np.abs(g_opt))[-15:][::-1]
print(f"\nTop 15 modes by |κ_k|:")
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
    f.write(f"# Coefficients from Quadratic Form fit (Internal only)\n")
    f.write(f"# J(ω) = g^T J_qq(ω) g\n")
    f.write(f"# Fitted with normalized J_internal and J_qq\n")
    f.write(f"# Normalization:\n")
    f.write(f"#   Common scale factor = {scale_factor:.6e}\n")
    f.write(f"#   (Both J_internal and J_qq normalized by J_target_rms)\n")
    f.write(f"# R² = {R2:.6f}\n")
    f.write(f"# Capture = {capture:.2f}%\n")
    f.write(f"# Fitting range: {freq_grid[0]:.0f}-{freq_grid[-1]:.0f} cm⁻¹\n")
    f.write(f"# Mode range: 500-1900 cm⁻¹\n")
    f.write(f"#\n")
    f.write("# Mode  Frequency(cm⁻¹)  g_k  |g_k|  Sign\n")
    for i in range(N_modes):
        sign = '+' if g_opt[i] > 0 else '-'
        f.write(f"{mode_indices[i]:4d}  {mode_freqs[i]:10.2f}  {g_opt[i]:12.6e}  {abs(g_opt[i]):12.6e}  {sign}\n")

print(f"Saved: g_k_quadratic_form.dat")

# Spectral densities
with open('spectral_density_quadratic_form.dat', 'w') as f:
    f.write("# Frequency(cm⁻¹)  J_internal(cm⁻¹²)  J_model(cm⁻¹²)  J_diag(cm⁻¹²)  J_offdiag(cm⁻¹²)  Residual(cm⁻¹²)\n")
    for nu, Jt, Jm, Jd, Jo, res in zip(freq_grid, J_target, J_model, J_diag, J_offdiag, residual):
        f.write(f"{nu:10.4f}  {Jt:14.6e}  {Jm:14.6e}  {Jd:14.6e}  {Jo:14.6e}  {res:14.6e}\n")

print(f"Saved: spectral_density_quadratic_form.dat")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\nQuadratic Form Optimization Summary:")
print(f"  Method: minimize ||J_internal(ω) - g^T J_qq(ω) g||²")
print(f"  Normalization: Common scale factor for numerical stability")
print(f"    Scale factor = {scale_factor:.6e}")
print(f"    (Both J_internal and J_qq normalized by J_target_rms)")
print(f"  R² = {R2:.4f}")
print(f"  Spectral capture = {capture:.1f}%")
print(f"  Active modes: {n_active}/{N_modes}")
print(f"    Positive: {n_positive}")
print(f"    Negative: {n_negative}")
print(f"  Negative J points: {n_neg}/{len(freq_grid)} ({100*n_neg/len(freq_grid):.1f}%)")
print(f"\n  Cache files created:")
print(f"    - {cache_file_J_qq}")
print(f"\n  Advantages:")
print(f"    ✓ Naturally includes mode coupling via J_qq[k,l]")
print(f"    ✓ Theoretically exact for linear model e(t) = g^T q(t)")
print(f"    ✓ No ad-hoc separation into stages")
print(f"    ✓ Normalized with common scale for consistent units")

# ============================================================================
# Plot results
# ============================================================================
print("\n" + "="*70)
print("Generating Plots")
print("="*70)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left panel: Spectral densities
ax = axes[0]
ax.plot(freq_grid, J_target * cm_to_meV, 'k-', linewidth=2, label='Internal')
ax.plot(freq_grid, J_model * cm_to_meV, 'r--', linewidth=1.5, label='Fit')
ax.plot(freq_grid, J_diag * cm_to_meV, 'b:', linewidth=1.5, label='Auto')
ax.plot(freq_grid, J_offdiag * cm_to_meV, 'g-.', linewidth=1.5, label='Cross')
ax.axhline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
ax.set_ylabel(r'Spectral Density (meV)', fontsize=16)
ax.set_xlim(freq_grid[0], freq_grid[-1])
# Allow negative values for cross-correlation
y_max = max(np.max(J_target), np.max(J_model)) * cm_to_meV * 1.1
y_min = min(0, np.min(J_offdiag) * cm_to_meV * 1.1)
ax.set_ylim(y_min, y_max)
ax.legend(fontsize=14, loc='upper left')
ax.grid(True, alpha=0.3)
ax.tick_params(axis='both', labelsize=16)
#ax.text(0.02, 0.98, f'R² = {R2:.4f}\nCapture = {capture:.1f}%', 
#ax.text(0.02, 0.98, f'R² = {R2:.4f}', 
#        transform=ax.transAxes, fontsize=14, verticalalignment='top',
#        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Right panel: Mode contributions (|κ_k|²)
ax = axes[1]
g_sq = g_opt**2
# Sort by contribution
sort_idx = np.argsort(g_sq)[::-1]

# Bar plot with colors based on sign of κ
colors = ['C0' if g_opt[i] > 0 else 'C3' for i in sort_idx[:20]]
ax.bar(range(min(20, len(sort_idx))), g_sq[sort_idx[:20]], color=colors, alpha=0.7)
ax.set_xlabel(r'Mode rank by $|\kappa|^2$', fontsize=16)
ax.set_ylabel(r'$|\kappa_k|^2$', fontsize=16)
ax.set_xlim(-0.5, min(20, len(sort_idx)) - 0.5)
ax.set_ylim(0, np.max(g_sq[sort_idx[:20]]) * 1.1)
ax.grid(True, alpha=0.3, axis='y')
ax.tick_params(axis='both', labelsize=16)

# Add frequency labels on top
for i, idx in enumerate(sort_idx[:min(10, len(sort_idx))]):
    ax.text(i, g_sq[idx], f'{mode_freqs[idx]:.0f}', 
            ha='center', va='bottom', fontsize=12, rotation=45)

ax.legend([plt.Rectangle((0,0),1,1,fc='C0',alpha=0.7), 
           plt.Rectangle((0,0),1,1,fc='C3',alpha=0.7)],
          [r'$\kappa > 0$', r'$\kappa < 0$'], fontsize=16, loc='upper right')

plt.tight_layout()
plt.savefig('spectral_fit_internal.png', dpi=300, bbox_inches='tight')
plt.savefig('spectral_fit_internal.pdf', bbox_inches='tight')
print("  Saved: spectral_fit_internal.png, .pdf")
plt.close()

# Additional plot: Residual
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(freq_grid, residual * cm_to_meV, 'k-', linewidth=1)
ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax.fill_between(freq_grid, 0, residual * cm_to_meV, alpha=0.3, 
                where=residual > 0, color='C0', label='Underfit')
ax.fill_between(freq_grid, 0, residual * cm_to_meV, alpha=0.3, 
                where=residual < 0, color='C3', label='Overfit')
ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
ax.set_ylabel(r'Residual (meV)', fontsize=16)
ax.set_xlim(freq_grid[0], freq_grid[-1])
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.tick_params(axis='both', labelsize=12)
plt.tight_layout()
plt.savefig('spectral_fit_residual.png', dpi=300, bbox_inches='tight')
plt.savefig('spectral_fit_residual.pdf', bbox_inches='tight')
print("  Saved: spectral_fit_residual.png, .pdf")
plt.close()

print("\nAll plots saved.")
