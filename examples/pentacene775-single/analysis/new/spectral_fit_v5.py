#!/usr/bin/env python3
"""
v5: 周波数グリッドを2010 cm⁻¹まで拡張して37,71の和周波数も確認

Stage 1: Linear terms with NNLS (c_k ≥ 0)
Stage 2: d_37,71 の4次相関でフィット (d²)
  - 差周波数 524 cm⁻¹ と 和周波数 2006 cm⁻¹ の両方を確認
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import nnls, least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("v5: 周波数グリッド拡張 (0-2010 cm⁻¹)")
print("="*70)

# Load data
print("\nLoading data...")
energy_diff_data = np.loadtxt('energy_diff_all.dat')
energy_diff = energy_diff_data[:, 4]

# Load mode information
mode_indices_all = []
mode_freqs_all = []
with open('mode_coords_new_stats.txt') as f:
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

# Load mode coordinates
q_k_data = np.loadtxt('mode_coords_new.txt')
q_k_all = q_k_data[:, 1:]
mode_coords = q_k_all[:, mode_indices - 7]

N_frames = min(len(energy_diff), len(mode_coords))
energy_diff = energy_diff[:N_frames]
mode_coords = mode_coords[:N_frames, :]
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
    return acf * window

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

def compute_correlation_4th(qi, qj, max_lag_corr=5000):
    """4次相関関数 C4(t) = <δ(q_i(t)q_j(t)) * δ(q_i(0)q_j(0))>"""
    n_frames = len(qi)
    
    qi_qj = qi * qj
    qi_qj_fluct = qi_qj - np.mean(qi_qj)
    
    C4 = np.zeros(max_lag_corr)
    for t in range(max_lag_corr):
        if t < n_frames:
            C4[t] = np.mean(qi_qj_fluct[t:] * qi_qj_fluct[:n_frames-t])
    
    t_arr = np.arange(max_lag_corr)
    window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag_corr))
    return C4 * window

def compute_spectrum_from_corr(C, freq):
    """相関関数からスペクトルを計算"""
    t = np.arange(len(C)) * dt
    
    spectrum = np.zeros(len(freq))
    for idx, nu in enumerate(freq):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * nu * t)
        spectrum[idx] = simpson(integrand, x=t)
    
    spectrum *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return spectrum

# ============================================================================
# Frequency grid (0-2010 cm⁻¹に拡張)
# ============================================================================
freq_grid = np.linspace(0, 2010, 2011)
print(f"\nFrequency grid: {len(freq_grid)} points (0-2010 cm⁻¹)")

# ============================================================================
# Stage 1: Linear terms with NNLS (全モード)
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, c_k ≥ 0)")
print("="*70)

cache_file_linear = 'spectral_cache_linear_v5.npz'

if os.path.exists(cache_file_linear):
    print(f"\nLoading cached linear spectral densities from {cache_file_linear}...")
    cache_data = np.load(cache_file_linear)
    if len(cache_data['freq_grid']) == len(freq_grid):
        J_total = cache_data['J_total']
        J_k_matrix = cache_data['J_k_matrix']
        J_k_integrals = cache_data['J_k_integrals']
    else:
        print("  Cache freq_grid mismatch, recomputing...")
        cache_data = None
else:
    cache_data = None

if cache_data is None or 'J_total' not in dir():
    print("\nComputing J_total...")
    J_total = compute_spectral_density(energy_diff, freq_grid)
    
    print("\nComputing J_k for all modes...")
    J_k_matrix = np.zeros((len(freq_grid), N_modes))
    J_k_integrals = np.zeros(N_modes)
    for k in range(N_modes):
        if (k+1) % 10 == 0:
            print(f"  Mode {k+1}/{N_modes}")
        J_k_matrix[:, k] = compute_spectral_density(mode_coords[:, k], freq_grid)
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

J_residual_1 = J_total - J_linear

ss_res_1 = np.sum(J_residual_1**2)
ss_tot = np.sum((J_total - np.mean(J_total))**2)
R2_linear = 1 - ss_res_1 / ss_tot

n_active = np.sum(c_k_physical > 1e-10)

print(f"\nStage 1 Results:")
print(f"  R² = {R2_linear:.4f}")
print(f"  ∫J_linear = {J_linear_integral:.6e} cm⁻¹²")
print(f"  Active modes: {n_active}/{N_modes}")

# ============================================================================
# Stage 2: d_37,71 をフィット（4次相関のみ）
# ============================================================================
print("\n" + "="*70)
print("Stage 2: d_37,71 をフィット（4次相関）")
print("="*70)

idx_37 = np.where(mode_indices == 37)[0][0]
idx_71 = np.where(mode_indices == 71)[0][0]

q37 = mode_coords[:, idx_37]
q71 = mode_coords[:, idx_71]

print(f"\nモード周波数:")
print(f"  Mode 37: {mode_freqs[idx_37]:.1f} cm⁻¹")
print(f"  Mode 71: {mode_freqs[idx_71]:.1f} cm⁻¹")
print(f"  和周波数 37+71: {mode_freqs[idx_37] + mode_freqs[idx_71]:.1f} cm⁻¹")
print(f"  差周波数 71-37: {mode_freqs[idx_71] - mode_freqs[idx_37]:.1f} cm⁻¹")

# 4次相関のスペクトルを計算
print("\n4次相関のスペクトルを計算中...")

cache_file_v5 = 'spectral_cache_v5.npz'

if os.path.exists(cache_file_v5):
    print(f"Loading cached from {cache_file_v5}...")
    cache_v5 = np.load(cache_file_v5)
    if len(cache_v5['freq_grid']) == len(freq_grid):
        J_4th_37_71 = cache_v5['J_4th_37_71']
    else:
        print("  Cache freq_grid mismatch, recomputing...")
        cache_v5 = None
else:
    cache_v5 = None

if cache_v5 is None or 'J_4th_37_71' not in dir():
    print("  Computing 4th(37,71)...")
    C4_37_71 = compute_correlation_4th(q37, q71)
    J_4th_37_71 = compute_spectrum_from_corr(C4_37_71, freq_grid)
    
    np.savez_compressed(cache_file_v5,
                        J_4th_37_71=J_4th_37_71,
                        freq_grid=freq_grid)
    print(f"  Saved to {cache_file_v5}")

# ピーク位置の確認
idx_diff = np.argmin(np.abs(freq_grid - 524))
idx_sum = np.argmin(np.abs(freq_grid - 2006))

print(f"\nJ_4th(37,71) のスペクトル特性:")
print(f"  差周波数 {freq_grid[idx_diff]:.0f} cm⁻¹: {J_4th_37_71[idx_diff]:.6e}")
print(f"  和周波数 {freq_grid[idx_sum]:.0f} cm⁻¹: {J_4th_37_71[idx_sum]:.6e}")
print(f"  和/差 比: {J_4th_37_71[idx_sum]/J_4th_37_71[idx_diff]:.4f}")

# フィッティング: d² * J_4th_37_71
def residual_func(d):
    J_model = d[0]**2 * J_4th_37_71
    return J_residual_1 - J_model

result = least_squares(residual_func, [0.1], method='lm')
d_37_71 = result.x[0]

print(f"\n非線形最小二乗法の結果:")
print(f"  d_37,71 = {d_37_71:+.6e}")
print(f"  d_37,71² = {d_37_71**2:+.6e}")

J_4th_contrib = d_37_71**2 * J_4th_37_71
J_model_best = J_linear + J_4th_contrib

ss_res_best = np.sum((J_total - J_model_best)**2)
R2_best = 1 - ss_res_best / ss_tot

print(f"\nR² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

# ============================================================================
# 各周波数での分析
# ============================================================================
print("\n" + "="*70)
print("各周波数での分析")
print("="*70)

for freq_target in [524, 2006]:
    idx = np.argmin(np.abs(freq_grid - freq_target))
    print(f"\n{freq_grid[idx]:.0f} cm⁻¹:")
    print(f"  J_total:      {J_total[idx]:.6e}")
    print(f"  J_linear:     {J_linear[idx]:.6e}")
    print(f"  J_4th(37,71): {J_4th_contrib[idx]:.6e}")
    print(f"  J_model:      {J_model_best[idx]:.6e}")
    print(f"  Residual:     {J_total[idx] - J_model_best[idx]:.6e}")

# ============================================================================
# 可視化
# ============================================================================
print("\n可視化を生成中...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: 全体のフィット
ax1 = axes[0, 0]
ax1.plot(freq_grid, J_total, 'b-', label='J_total (MD)', linewidth=1.5)
ax1.plot(freq_grid, J_linear, 'g--', label='J_linear (Stage 1)', linewidth=1.5)
ax1.plot(freq_grid, J_model_best, 'r-', label='J_model (Stage 1+2)', linewidth=1.5)
ax1.axvline(524, color='gray', linestyle=':', alpha=0.5, label='524 cm⁻¹')
ax1.axvline(2006, color='orange', linestyle=':', alpha=0.5, label='2006 cm⁻¹')
ax1.set_xlabel('Frequency (cm⁻¹)')
ax1.set_ylabel('J(ω)')
ax1.set_title(f'Overall Fit (R² = {R2_best:.4f})')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# Plot 2: 差周波数 524 cm⁻¹ 付近
ax2 = axes[0, 1]
mask_524 = (freq_grid >= 480) & (freq_grid <= 580)
ax2.plot(freq_grid[mask_524], J_total[mask_524], 'b-', label='J_total', linewidth=2)
ax2.plot(freq_grid[mask_524], J_linear[mask_524], 'g--', label='Stage 1', linewidth=2)
ax2.plot(freq_grid[mask_524], J_model_best[mask_524], 'r-', label='Stage 1+2', linewidth=2)
ax2.plot(freq_grid[mask_524], J_4th_contrib[mask_524], 'orange', label='4th(37,71)', linewidth=1.5, alpha=0.7)
ax2.axvline(524, color='gray', linestyle=':', linewidth=1)
ax2.set_xlabel('Frequency (cm⁻¹)')
ax2.set_ylabel('J(ω)')
ax2.set_title('差周波数 524 cm⁻¹ (71-37)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: 和周波数 2006 cm⁻¹ 付近
ax3 = axes[1, 0]
mask_2006 = (freq_grid >= 1950) & (freq_grid <= 2010)
ax3.plot(freq_grid[mask_2006], J_total[mask_2006], 'b-', label='J_total', linewidth=2)
ax3.plot(freq_grid[mask_2006], J_linear[mask_2006], 'g--', label='Stage 1', linewidth=2)
ax3.plot(freq_grid[mask_2006], J_model_best[mask_2006], 'r-', label='Stage 1+2', linewidth=2)
ax3.plot(freq_grid[mask_2006], J_4th_contrib[mask_2006], 'orange', label='4th(37,71)', linewidth=1.5, alpha=0.7)
ax3.axvline(2006, color='gray', linestyle=':', linewidth=1)
ax3.set_xlabel('Frequency (cm⁻¹)')
ax3.set_ylabel('J(ω)')
ax3.set_title('和周波数 2006 cm⁻¹ (37+71)')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: J_4th(37,71) の全スペクトル
ax4 = axes[1, 1]
ax4.plot(freq_grid, J_4th_37_71, 'b-', label='J_4th(37,71)', linewidth=1.5)
ax4.axvline(524, color='orange', linestyle='--', alpha=0.7, label='差 524 cm⁻¹')
ax4.axvline(2006, color='red', linestyle='--', alpha=0.7, label='和 2006 cm⁻¹')
ax4.set_xlabel('Frequency (cm⁻¹)')
ax4.set_ylabel('J_4th(37,71)')
ax4.set_title('4th order correlation spectrum (37,71)')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectral_fit_v5.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_fit_v5.png")
plt.close()

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Stage 1 (Linear NNLS): R² = {R2_linear:.4f}")
print(f"Stage 2 (+d_37,71²):   R² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")
print(f"\nd_37,71 = {d_37_71:+.6e}")
print(f"d_37,71² = {d_37_71**2:+.6e}")

print(f"\nJ_4th(37,71) のピーク:")
print(f"  差周波数 524 cm⁻¹:  {J_4th_37_71[idx_diff]:.6e}")
print(f"  和周波数 2006 cm⁻¹: {J_4th_37_71[idx_sum]:.6e}")
print(f"  和/差 比: {J_4th_37_71[idx_sum]/J_4th_37_71[idx_diff]:.4f}")

print(f"\nフィット後の寄与 (d²×J_4th):")
print(f"  差周波数 524 cm⁻¹:  {J_4th_contrib[idx_diff]:.6e}")
print(f"  和周波数 2006 cm⁻¹: {J_4th_contrib[idx_sum]:.6e}")

print("\n" + "="*70)
print("Complete!")
print("="*70)
