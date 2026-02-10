#!/usr/bin/env python3
"""
v6: Stage 1の残差をモード25, 26, 36の4次相関でフィット

Stage 1: Linear terms with NNLS (c_k ≥ 0)
Stage 2: モード25, 26, 36の全4次相関ペア(6組)で非線形最小二乗フィット
  - 4次相関: <δ(q_i(t)q_j(t)) * δ(q_i(0)q_j(0))> (i,j ∈ {25,26,36})
  - 係数は d² (非線形)
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import nnls
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
print("v6: Stage 1残差をモード25,26,36の4次相関でフィット")
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

def compute_correlation_3rd(qi, qj, qk, max_lag_corr=5000):
    """3次相関関数 C3(t) = <δ(q_i(t)q_j(t)) * δq_k(0)> (揺らぎの相関)"""
    n_frames = len(qi)
    
    # 揺らぎを計算
    qi_qj = qi * qj
    qi_qj_fluct = qi_qj - np.mean(qi_qj)
    qk_fluct = qk - np.mean(qk)
    
    C3 = np.zeros(max_lag_corr)
    for t in range(max_lag_corr):
        if t < n_frames:
            C3[t] = np.mean(qi_qj_fluct[t:] * qk_fluct[:n_frames-t])
    
    # Hanning窓を適用
    t_arr = np.arange(max_lag_corr)
    window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag_corr))
    return C3 * window

def compute_correlation_4th(qi, qj, max_lag_corr=5000):
    """4次相関関数 C4(t) = <δ(q_i(t)q_j(t)) * δ(q_i(0)q_j(0))> (揺らぎの相関)"""
    n_frames = len(qi)
    
    # 揺らぎを計算
    qi_qj = qi * qj
    qi_qj_fluct = qi_qj - np.mean(qi_qj)
    
    C4 = np.zeros(max_lag_corr)
    for t in range(max_lag_corr):
        if t < n_frames:
            C4[t] = np.mean(qi_qj_fluct[t:] * qi_qj_fluct[:n_frames-t])
    
    # Hanning窓を適用
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
    
    # 温度ファクター
    spectrum *= (2 * np.pi * c_cm_fs * freq) / (kB_cm * T)
    return spectrum

# ============================================================================
# Frequency grid (0-1800 cm⁻¹)
# ============================================================================
freq_grid = np.linspace(0, 1800, 1801)
print(f"\nFrequency grid: {len(freq_grid)} points (0-1800 cm⁻¹)")

# ============================================================================
# Stage 1: Linear terms with NNLS (全モード 0-1800 cm⁻¹)
# ============================================================================
print("\n" + "="*70)
print("Stage 1: Linear Terms (NNLS, c_k ≥ 0) - 全モード使用")
print("="*70)

# Cache for Stage 1 (v4用の新しいキャッシュ)
cache_file_linear = 'spectral_cache_linear_v4.npz'

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

J_residual_1 = J_total - J_linear
J_residual_1_integral = np.trapezoid(J_residual_1, freq_grid)

ss_res_1 = np.sum(J_residual_1**2)
ss_tot = np.sum((J_total - np.mean(J_total))**2)
R2_linear = 1 - ss_res_1 / ss_tot

n_active = np.sum(c_k_physical > 1e-10)

print(f"\nStage 1 Results:")
print(f"  R² = {R2_linear:.4f}")
print(f"  ∫J_linear = {J_linear_integral:.6e} cm⁻¹² ({100*J_linear_integral/J_total_integral:.1f}%)")
print(f"  ∫J_residual = {J_residual_1_integral:.6e} cm⁻¹² ({100*J_residual_1_integral/J_total_integral:.1f}%)")
print(f"  Active modes: {n_active}/{N_modes}")

# ============================================================================
# Stage 2: モード25,26,36の4次相関をフィット
# ============================================================================
print("\n" + "="*70)
print("Stage 2: モード25,26,36の4次相関をフィット")
print("="*70)

# モードのインデックス（配列インデックスはモード番号-7）
idx_25 = 25 - 7  # = 18
idx_26 = 26 - 7  # = 19
idx_36 = 36 - 7  # = 29

q25 = mode_coords[:, idx_25]
q26 = mode_coords[:, idx_26]
q36 = mode_coords[:, idx_36]

print(f"\nモード周波数:")
print(f"  Mode 25: {mode_freqs[idx_25]:.1f} cm⁻¹")
print(f"  Mode 26: {mode_freqs[idx_26]:.1f} cm⁻¹")
print(f"  Mode 36: {mode_freqs[idx_36]:.1f} cm⁻¹")

print(f"\n期待されるピーク位置:")
print(f"  (25,25): 2×{mode_freqs[idx_25]:.1f} = {2*mode_freqs[idx_25]:.1f} cm⁻¹")
print(f"  (25,26): {mode_freqs[idx_25]:.1f}+{mode_freqs[idx_26]:.1f} = {mode_freqs[idx_25]+mode_freqs[idx_26]:.1f} cm⁻¹")
print(f"  (25,36): {mode_freqs[idx_25]:.1f}+{mode_freqs[idx_36]:.1f} = {mode_freqs[idx_25]+mode_freqs[idx_36]:.1f} cm⁻¹")
print(f"  (26,26): 2×{mode_freqs[idx_26]:.1f} = {2*mode_freqs[idx_26]:.1f} cm⁻¹")
print(f"  (26,36): {mode_freqs[idx_26]:.1f}+{mode_freqs[idx_36]:.1f} = {mode_freqs[idx_26]+mode_freqs[idx_36]:.1f} cm⁻¹")
print(f"  (36,36): 2×{mode_freqs[idx_36]:.1f} = {2*mode_freqs[idx_36]:.1f} cm⁻¹")

# ============================================================================
# 4次相関のスペクトルを計算（6ペア）
# ============================================================================
print("\n4次相関のスペクトルを計算中...")

cache_file_v6 = 'spectral_cache_v6_25_26_36.npz'

# ペアのリスト
pairs = [(25,25), (25,26), (25,36), (26,26), (26,36), (36,36)]
pair_names = ['(25,25)', '(25,26)', '(25,36)', '(26,26)', '(26,36)', '(36,36)']
q_dict = {25: q25, 26: q26, 36: q36}

if os.path.exists(cache_file_v6):
    print(f"Loading cached from {cache_file_v6}...")
    cache_v6 = np.load(cache_file_v6)
    if len(cache_v6['freq_grid']) != len(freq_grid):
        print("  Cache freq_grid mismatch, recomputing...")
        os.remove(cache_file_v6)
        cache_v6 = None
    else:
        J_4th = {}
        for i, j in pairs:
            key = f'J_4th_{i}_{j}'
            J_4th[(i,j)] = cache_v6[key]
else:
    cache_v6 = None

if cache_v6 is None:
    J_4th = {}
    for i, j in pairs:
        print(f"  Computing 4th({i},{j})...")
        C4 = compute_correlation_4th(q_dict[i], q_dict[j])
        J_4th[(i,j)] = compute_spectrum_from_corr(C4, freq_grid)
    
    save_dict = {'freq_grid': freq_grid}
    for i, j in pairs:
        save_dict[f'J_4th_{i}_{j}'] = J_4th[(i,j)]
    np.savez_compressed(cache_file_v6, **save_dict)
    print(f"  Saved to {cache_file_v6}")

print(f"\n各4次相関スペクトルの特性:")
for (i,j), name in zip(pairs, pair_names):
    J = J_4th[(i,j)]
    idx_peak = np.argmax(np.abs(J))
    print(f"  4th{name}: max={np.max(np.abs(J)):.6e}, peak={freq_grid[idx_peak]:.0f} cm⁻¹")

# ============================================================================
# Stage 2 フィッティング（非線形: d² の係数）
# J_cross = Σ d_ij² * J_4th(i,j)
# ============================================================================
print("\n" + "="*70)
print("Stage 2 フィッティング (6つの d 係数, 非線形)")
print("="*70)

from scipy.optimize import least_squares

# J_4th をリストに変換
J_4th_list = [J_4th[p] for p in pairs]

def residual_func(d_params):
    J_model = np.zeros_like(freq_grid)
    for k, d in enumerate(d_params):
        J_model += d**2 * J_4th_list[k]
    return J_residual_1 - J_model

# 初期値
d0 = [0.01] * 6
result = least_squares(residual_func, d0, method='lm')
d_best = result.x

print(f"\n非線形フィット結果 (d, d²):")
for k, name in enumerate(pair_names):
    print(f"  d{name} = {d_best[k]:+.6e}, d² = {d_best[k]**2:.6e}")

# 各成分の寄与
J_4th_contrib = {}
J_cross_best = np.zeros_like(freq_grid)
for k, ((i,j), name) in enumerate(zip(pairs, pair_names)):
    J_4th_contrib[(i,j)] = d_best[k]**2 * J_4th[(i,j)]
    J_cross_best += J_4th_contrib[(i,j)]

J_model_best = J_linear + J_cross_best
ss_res_best = np.sum((J_total - J_model_best)**2)
R2_best = 1 - ss_res_best / ss_tot

print(f"\nR² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

# ============================================================================
# 可視化
# ============================================================================
print("\n可視化を生成中...")

fig, axes = plt.subplots(3, 2, figsize=(14, 15))

# Plot 1: J_total vs J_linear vs J_model
ax1 = axes[0, 0]
ax1.plot(freq_grid, J_total, 'b-', label='J_total (MD)', linewidth=1.5)
ax1.plot(freq_grid, J_linear, 'g--', label='J_linear (Stage 1)', linewidth=1.5)
ax1.plot(freq_grid, J_model_best, 'r-', label='J_model (Stage 1+2)', linewidth=1.5)
ax1.set_xlabel('Frequency (cm⁻¹)')
ax1.set_ylabel('J(ω) (cm⁻¹²)')
ax1.set_title(f'Overall Fit (R² = {R2_best:.4f})')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: 残差の比較
ax2 = axes[0, 1]
ax2.plot(freq_grid, J_residual_1, 'g-', label='Stage 1 residual', linewidth=1)
ax2.plot(freq_grid, J_total - J_model_best, 'r-', label='Stage 2 residual', linewidth=1)
ax2.axhline(0, color='k', linestyle='--', linewidth=0.5)
ax2.set_xlabel('Frequency (cm⁻¹)')
ax2.set_ylabel('Residual (cm⁻¹²)')
ax2.set_title('Residual comparison')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: 1230 cm⁻¹付近の拡大
ax3 = axes[1, 0]
mask_1230 = (freq_grid >= 1180) & (freq_grid <= 1280)
ax3.plot(freq_grid[mask_1230], J_total[mask_1230], 'b-', label='J_total', linewidth=2)
ax3.plot(freq_grid[mask_1230], J_linear[mask_1230], 'g--', label='Stage 1', linewidth=2)
ax3.plot(freq_grid[mask_1230], J_model_best[mask_1230], 'r-', label='Stage 1+2', linewidth=2)
ax3.axvline(1230, color='gray', linestyle=':', linewidth=1)
ax3.set_xlabel('Frequency (cm⁻¹)')
ax3.set_ylabel('J(ω) (cm⁻¹²)')
ax3.set_title('Zoom: 1230 cm⁻¹ region')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: 各ペアの4次相関の寄与
ax4 = axes[1, 1]
colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
for k, ((i,j), name) in enumerate(zip(pairs, pair_names)):
    ax4.plot(freq_grid, J_4th_contrib[(i,j)], color=colors[k], label=f'4th{name}', linewidth=1, alpha=0.7)
ax4.plot(freq_grid, J_cross_best, 'k--', linewidth=2, label='Total')
ax4.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax4.set_xlabel('Frequency (cm⁻¹)')
ax4.set_ylabel('J(ω)')
ax4.set_title('4th order contributions by pair')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.3)

# Plot 5: 1200-1280 cm⁻¹付近の拡大（和周波数領域）
ax5 = axes[2, 0]
mask_1220 = (freq_grid >= 1180) & (freq_grid <= 1280)
ax5.plot(freq_grid[mask_1220], J_total[mask_1220], 'b-', label='J_total', linewidth=2)
ax5.plot(freq_grid[mask_1220], J_linear[mask_1220], 'g--', label='Stage 1', linewidth=2)
ax5.plot(freq_grid[mask_1220], J_model_best[mask_1220], 'r-', label='Stage 1+2', linewidth=2)
for k, ((i,j), name) in enumerate(zip(pairs, pair_names)):
    ax5.plot(freq_grid[mask_1220], J_4th_contrib[(i,j)][mask_1220], color=colors[k], 
             label=f'4th{name}', linewidth=1, alpha=0.5)
ax5.axvline(1230, color='gray', linestyle=':', linewidth=1)
ax5.set_xlabel('Frequency (cm⁻¹)')
ax5.set_ylabel('J(ω) (cm⁻¹²)')
ax5.set_title('Zoom: 1220 cm⁻¹ region (sum frequencies 25/26+36)')
ax5.legend(fontsize=6)
ax5.grid(True, alpha=0.3)

# Plot 6: Stage 1 残差と cross term の比較
ax6 = axes[2, 1]
ax6.plot(freq_grid, J_residual_1, 'g-', label='Stage 1 residual', linewidth=1.5)
ax6.plot(freq_grid, J_cross_best, 'r-', label='Cross term (3rd+4th)', linewidth=1.5)
ax6.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax6.set_xlabel('Frequency (cm$^{-1}$)')
ax6.set_ylabel('J($\\omega$)')
ax6.set_title('Stage 1 residual vs Cross term')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectral_fit_v6.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_fit_v6.png")
plt.close()

# ============================================================================
# 1230 cm⁻¹ における分析
# ============================================================================
print("\n" + "="*70)
print("1230 cm⁻¹ における分析（4次相関の寄与）")
print("="*70)

idx_1230 = np.argmin(np.abs(freq_grid - 1230))
print(f"\n周波数: {freq_grid[idx_1230]:.1f} cm⁻¹")
print(f"  J_total:      {J_total[idx_1230]:.6e}")
print(f"  J_linear:     {J_linear[idx_1230]:.6e}")
for (i,j), name in zip(pairs, pair_names):
    print(f"  J_4th{name}: {J_4th_contrib[(i,j)][idx_1230]:.6e}")
print(f"  J_model:      {J_model_best[idx_1230]:.6e}")
print(f"  Residual:     {J_total[idx_1230] - J_model_best[idx_1230]:.6e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Stage 1 (Linear NNLS, all modes): R² = {R2_linear:.4f}")
print(f"Stage 2 (4th order, 6 pairs):     R² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

print(f"\n非線形係数 d (d²):")
for k, name in enumerate(pair_names):
    if np.abs(d_best[k]) > 1e-6:
        print(f"  d{name} = {d_best[k]:+.6e} ({d_best[k]**2:.6e})")

print(f"\n上位線形係数:")
sorted_idx_linear = np.argsort(c_k_physical)[::-1]
for rank, i in enumerate(sorted_idx_linear[:10], 1):
    if c_k_physical[i] > 1e-10:
        print(f"  {rank}. Mode {i+7}: c = {c_k_physical[i]:.6e}, freq = {mode_freqs[i]:.1f} cm⁻¹")

print("\n" + "="*70)
print("Complete!")
print("="*70)
