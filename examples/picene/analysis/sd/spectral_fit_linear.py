#!/usr/bin/env python3
"""
Spectral density fitting with linear terms only (NNLS).

Stage 1: Linear terms with NNLS (c_k ≥ 0)
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import nnls
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 3000  # frames

print("="*70)
print("Spectral Density Fitting: Linear Terms Only (NNLS)")
print("="*70)


def interpolate_nan(data, name="data"):
    """Interpolate NaN values using linear interpolation."""
    nan_mask = np.isnan(data)
    n_nan = np.sum(nan_mask)
    
    if n_nan == 0:
        return data
    
    print(f"  Found {n_nan} NaN values in {name} ({100*n_nan/len(data):.2f}%), interpolating...")
    
    valid_mask = ~nan_mask
    if np.sum(valid_mask) < 2:
        raise ValueError(f"Not enough valid data points for interpolation in {name}")
    
    indices = np.arange(len(data))
    f_interp = interp1d(indices[valid_mask], data[valid_mask], 
                        kind='linear', fill_value='extrapolate')
    
    data_interp = data.copy()
    data_interp[nan_mask] = f_interp(indices[nan_mask])
    
    return data_interp


# Load data
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff.dat')
energy_diff = energy_diff_data[:, 2]  # 3列目：原子単位のエネルギー差

# Interpolate NaN values in energy_diff
energy_diff = interpolate_nan(energy_diff, "energy_diff")

# Load mode information
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

mode_indices = np.array(mode_indices_all)
mode_freqs = np.array(mode_freqs_all)

# Load normalized displacement
normalized_data = np.loadtxt('normalized_displacement.txt')
# 最初の列はフレームインデックスなのでスキップ
normalized_coords = normalized_data[:, 1:]  # モード7から始まる

# Interpolate NaN values in normalized_coords (each mode separately)
for k in range(normalized_coords.shape[1]):
    normalized_coords[:, k] = interpolate_nan(normalized_coords[:, k], f"mode_{k+7}")

N_frames = min(len(energy_diff), len(normalized_coords))
energy_diff = energy_diff[:N_frames]
normalized_coords = normalized_coords[:N_frames, :]
N_modes = len(mode_freqs)

print(f"  Frames: {N_frames}")
print(f"  Modes: {N_modes}")

# ============================================================================
# 相関関数の計算
# ============================================================================

def compute_autocorr(x):
    """Compute autocorrelation with Hanning window"""
    x_mean = np.mean(x)
    x_centered = x - x_mean
    N = len(x)
    
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(x_centered[:N-lag] * x_centered[lag:])
    
    t = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t / max_lag))
    return acf #* window

def compute_spectral_density(x, freq):
    """Compute spectral density using Simpson integration"""
    C = compute_autocorr(x)
    t = np.arange(max_lag) * dt
    
    J = np.zeros(len(freq))
    for i, nu in enumerate(freq):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        J[i] = simpson(integrand, t)
    
    J *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return J

# ============================================================================
# Frequency grid (0-1800 cm⁻¹)
# ============================================================================
freq_grid = np.linspace(0, 1800, 1801)
print(f"\nFrequency grid: {len(freq_grid)} points (0-1800 cm⁻¹)")

# ============================================================================
# Stage 1: Linear terms with NNLS (全モード 0-1800 cm⁻¹)
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, c_k ≥ 0)")
print("="*70)

# Cache for Stage 1
cache_file_linear = 'spectral_cache_linear_normalized.npz'

if os.path.exists(cache_file_linear):
    print(f"\nLoading cached linear spectral densities from {cache_file_linear}...")
    cache_data = np.load(cache_file_linear)
    J_total = cache_data['J_total']
    J_k_matrix = cache_data['J_k_matrix']
    J_k_integrals = cache_data['J_k_integrals']
    # 周波数グリッドの互換性を確認
    if len(cache_data['freq_grid']) != len(freq_grid):
        print("  Cache freq_grid mismatch, recomputing...")
        os.remove(cache_file_linear)
        cache_data = None
else:
    cache_data = None

if cache_data is None:
    print("\nComputing J_total...")
    J_total = compute_spectral_density(energy_diff, freq_grid)
    
    print("\nComputing J_k for all modes...")
    J_k_matrix = np.zeros((len(freq_grid), N_modes))
    J_k_integrals = np.zeros(N_modes)
    for k in range(N_modes):
        if (k+1) % 10 == 0:
            print(f"  Mode {k+1}/{N_modes}")
        J_k_matrix[:, k] = compute_spectral_density(normalized_coords[:, k], freq_grid)
        J_k_integrals[k] = np.trapezoid(J_k_matrix[:, k], freq_grid)
    
    np.savez_compressed(cache_file_linear,
                        J_total=J_total,
                        J_k_matrix=J_k_matrix,
                        J_k_integrals=J_k_integrals,
                        freq_grid=freq_grid)
    print(f"  Saved to {cache_file_linear}")

J_total_integral = np.trapezoid(J_total, freq_grid)
print(f"  ∫J_total = {J_total_integral:.6e} cm⁻¹²")

# Stage 1: 全モード(0-1800 cm⁻¹)でフィット
stage1_mask = (mode_freqs >= 0) & (mode_freqs <= 1800)
print(f"\nUsing {np.sum(stage1_mask)} modes with 0 ≤ freq ≤ 1800 cm⁻¹")

J_k_stage1 = J_k_matrix[:, stage1_mask]
J_k_integrals_stage1 = J_k_integrals[stage1_mask]

# 積分がゼロのモードを除外
valid_integrals = J_k_integrals_stage1 > 1e-15
J_k_stage1 = J_k_stage1[:, valid_integrals]
J_k_integrals_stage1 = J_k_integrals_stage1[valid_integrals]
stage1_indices = np.where(stage1_mask)[0][valid_integrals]

print(f"  Valid modes (non-zero integral): {len(stage1_indices)}")

# Normalize and fit
J_k_normalized = J_k_stage1 / J_k_integrals_stage1[np.newaxis, :]
J_total_normalized = J_total / J_total_integral

c_k_nnls_stage1, _ = nnls(J_k_normalized, J_total_normalized)
c_k_physical_stage1 = c_k_nnls_stage1 * J_total_integral / J_k_integrals_stage1

c_k_physical = np.zeros(N_modes)
c_k_physical[stage1_indices] = c_k_physical_stage1

J_linear = J_k_matrix @ c_k_physical
J_linear_integral = np.trapezoid(J_linear, freq_grid)

J_residual = J_total - J_linear
J_residual_integral = np.trapezoid(J_residual, freq_grid)

ss_res = np.sum(J_residual**2)
ss_tot = np.sum((J_total - np.mean(J_total))**2)
R2_linear = 1 - ss_res / ss_tot

n_active = np.sum(c_k_physical > 1e-10)

print(f"\nResults:")
print(f"  R² = {R2_linear:.4f}")
print(f"  ∫J_linear = {J_linear_integral:.6e} cm⁻¹² ({100*J_linear_integral/J_total_integral:.1f}%)")
print(f"  ∫J_residual = {J_residual_integral:.6e} cm⁻¹² ({100*J_residual_integral/J_total_integral:.1f}%)")
print(f"  Active modes: {n_active}/{N_modes}")

# ============================================================================
# 可視化
# ============================================================================
print("\n可視化を生成中...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: J_total vs J_linear
ax1 = axes[0, 0]
ax1.plot(freq_grid, J_total, 'b-', label='J_total (MD)', linewidth=1.5)
ax1.plot(freq_grid, J_linear, 'r--', label='J_linear (fit)', linewidth=1.5)
ax1.set_xlabel('Frequency (cm⁻¹)')
ax1.set_ylabel('J(ω) (cm⁻¹²)')
ax1.set_title(f'Linear Fit (R² = {R2_linear:.4f})')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: 残差
ax2 = axes[0, 1]
ax2.plot(freq_grid, J_residual, 'g-', linewidth=1)
ax2.axhline(0, color='k', linestyle='--', linewidth=0.5)
ax2.set_xlabel('Frequency (cm⁻¹)')
ax2.set_ylabel('Residual (cm⁻¹²)')
ax2.set_title('Residual (J_total - J_linear)')
ax2.grid(True, alpha=0.3)

# Plot 3: 係数の棒グラフ
ax3 = axes[1, 0]
active_mask = c_k_physical > 1e-10
ax3.bar(mode_freqs[active_mask], c_k_physical[active_mask], width=10, alpha=0.7)
ax3.set_xlabel('Mode Frequency (cm⁻¹)')
ax3.set_ylabel('Coefficient')
ax3.set_title(f'Linear Coefficients ({n_active} active modes)')
ax3.grid(True, alpha=0.3)

# Plot 4: 各モードの寄与（上位10モード）
ax4 = axes[1, 1]
sorted_idx = np.argsort(c_k_physical)[::-1]
colors = plt.cm.tab10(np.linspace(0, 1, 10))
for rank, i in enumerate(sorted_idx[:10]):
    if c_k_physical[i] > 1e-10:
        ax4.plot(freq_grid, c_k_physical[i] * J_k_matrix[:, i], 
                 color=colors[rank], label=f'Mode {i+7} ({mode_freqs[i]:.0f} cm⁻¹)', 
                 linewidth=1, alpha=0.7)
ax4.set_xlabel('Frequency (cm⁻¹)')
ax4.set_ylabel('J_k contribution')
ax4.set_title('Top 10 mode contributions')
ax4.legend(fontsize=7, loc='upper right')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectral_fit_linear.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_fit_linear.png")
plt.close()

# ============================================================================
# 係数をファイル出力
# ============================================================================
print("\n" + "="*70)
print("係数をファイル出力")
print("="*70)

coef_file = 'spectral_fit_linear_coefficients.txt'
with open(coef_file, 'w') as f:
    f.write("# Linear NNLS coefficients\n")
    f.write(f"# R² = {R2_linear:.6f}\n")
    f.write(f"# ∫J_linear = {J_linear_integral:.6e} cm⁻¹²\n")
    f.write(f"# Active modes: {n_active}/{N_modes}\n")
    f.write("#\n")
    f.write("# Mode  Freq(cm⁻¹)  Coefficient\n")
    for i in range(N_modes):
        mode_num = i + 7
        f.write(f"{mode_num:5d}  {mode_freqs[i]:10.4f}  {c_k_physical[i]:15.6e}\n")
print(f"Saved: {coef_file}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Linear fit R² = {R2_linear:.4f}")
print(f"Active modes: {n_active}/{N_modes}")

print(f"\n上位線形係数:")
sorted_idx_linear = np.argsort(c_k_physical)[::-1]
for rank, i in enumerate(sorted_idx_linear[:10], 1):
    if c_k_physical[i] > 1e-10:
        print(f"  {rank}. Mode {i+7}: c = {c_k_physical[i]:.6e}, freq = {mode_freqs[i]:.1f} cm⁻¹")

print("\n" + "="*70)
print("Complete!")
print("="*70)
