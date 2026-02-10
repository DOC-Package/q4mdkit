#!/usr/bin/env python3
"""
v7: Stage 1の残差を17モードの4次相関でフィット (1100 cm⁻¹領域)

4次相関上位30に登場する17モード:
  低周波: 13, 14, 15, 16, 17, 28
  中周波: 35, 36, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51

Stage 1: Linear terms with NNLS (c_k ≥ 0)
Stage 2: 17モードの全4次相関ペアで非線形最小二乗フィット
  - 4次相関: <δ(q_i(t)q_j(t)) * δ(q_i(0)q_j(0))>
  - 係数は d² (非線形)
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from itertools import combinations_with_replacement

# Constants
kB_cm = 0.695034800  # cm⁻¹/K
c_cm_fs = 2.99792458e-5  # cm/fs
T = 300.0  # K
dt = 4.0  # fs
max_lag = 5000  # frames

print("="*70)
print("v7: Stage 1残差を17モードの4次相関でフィット (1100 cm⁻¹領域)")
print("="*70)

# 17モードのリスト
target_modes = [13, 14, 15, 16, 17, 28, 35, 36, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51]
print(f"\nターゲットモード ({len(target_modes)}個): {target_modes}")

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
print(f"  Total modes: {N_modes}")

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
# Stage 2: 17モードの4次相関をフィット
# ============================================================================
print("\n" + "="*70)
print("Stage 2: 17モードの4次相関をフィット (1100 cm⁻¹領域)")
print("="*70)

# モードのインデックス（配列インデックスはモード番号-7）
mode_indices_map = {m: m - 7 for m in target_modes}

# モード座標を辞書に格納
q_dict = {}
for m in target_modes:
    q_dict[m] = mode_coords[:, mode_indices_map[m]]

print(f"\nモード周波数:")
for m in target_modes:
    idx = mode_indices_map[m]
    print(f"  Mode {m}: {mode_freqs[idx]:.1f} cm⁻¹")

# 全ペアを生成 (i <= j)
pairs = list(combinations_with_replacement(target_modes, 2))
n_pairs = len(pairs)
print(f"\n4次相関ペア数: {n_pairs}")

# 期待されるピーク位置（1050-1150 cm⁻¹ 付近のもの）
print(f"\n1050-1150 cm⁻¹付近の和周波数を持つペア:")
relevant_pairs = []
for i, j in pairs:
    freq_i = mode_freqs[mode_indices_map[i]]
    freq_j = mode_freqs[mode_indices_map[j]]
    sum_freq = freq_i + freq_j
    if 1050 <= sum_freq <= 1150:
        print(f"  ({i},{j}): {freq_i:.1f}+{freq_j:.1f} = {sum_freq:.1f} cm⁻¹")
        relevant_pairs.append((i, j, sum_freq))

# ============================================================================
# 4次相関のスペクトルを計算
# ============================================================================
print("\n4次相関のスペクトルを計算中...")

cache_file_v7 = 'spectral_cache_v7_with39.npz'

if os.path.exists(cache_file_v7):
    print(f"Loading cached from {cache_file_v7}...")
    cache_v7 = np.load(cache_file_v7)
    if len(cache_v7['freq_grid']) != len(freq_grid):
        print("  Cache freq_grid mismatch, recomputing...")
        os.remove(cache_file_v7)
        cache_v7 = None
    else:
        J_4th = {}
        for i, j in pairs:
            key = f'J_4th_{i}_{j}'
            if key in cache_v7:
                J_4th[(i,j)] = cache_v7[key]
            else:
                print(f"  Key {key} not found, recomputing all...")
                os.remove(cache_file_v7)
                cache_v7 = None
                break
else:
    cache_v7 = None

if cache_v7 is None:
    J_4th = {}
    for idx, (i, j) in enumerate(pairs):
        if (idx + 1) % 20 == 0:
            print(f"  Computing pair {idx+1}/{n_pairs}...")
        C4 = compute_correlation_4th(q_dict[i], q_dict[j])
        J_4th[(i,j)] = compute_spectrum_from_corr(C4, freq_grid)
    
    save_dict = {'freq_grid': freq_grid}
    for i, j in pairs:
        save_dict[f'J_4th_{i}_{j}'] = J_4th[(i,j)]
    np.savez_compressed(cache_file_v7, **save_dict)
    print(f"  Saved to {cache_file_v7}")

print(f"\n4次相関スペクトルの計算完了: {len(J_4th)} ペア")

# ============================================================================
# Stage 2 フィッティング（非線形: d² の係数）
# J_cross = Σ d_ij² * J_4th(i,j)
# ============================================================================
print("\n" + "="*70)
print(f"Stage 2 フィッティング ({n_pairs}個の d 係数, 非線形)")
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
d0 = [0.01] * n_pairs
result = least_squares(residual_func, d0, method='lm')
d_best = result.x

# 各成分の寄与
J_4th_contrib = {}
J_cross_best = np.zeros_like(freq_grid)
for k, (i, j) in enumerate(pairs):
    J_4th_contrib[(i,j)] = d_best[k]**2 * J_4th[(i,j)]
    J_cross_best += J_4th_contrib[(i,j)]

J_model_best = J_linear + J_cross_best
ss_res_best = np.sum((J_total - J_model_best)**2)
R2_best = 1 - ss_res_best / ss_tot

print(f"\nR² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

# 上位の非線形係数を表示
print(f"\n非線形係数 d² (上位20):")
d_squared = [(k, pairs[k], d_best[k], d_best[k]**2) for k in range(n_pairs)]
d_squared_sorted = sorted(d_squared, key=lambda x: abs(x[3]), reverse=True)
for rank, (k, (i, j), d, d2) in enumerate(d_squared_sorted[:20], 1):
    freq_i = mode_freqs[mode_indices_map[i]]
    freq_j = mode_freqs[mode_indices_map[j]]
    sum_freq = freq_i + freq_j
    print(f"  {rank:2d}. ({i},{j}): d={d:+.4e}, d²={d2:.4e}, sum={sum_freq:.0f}cm⁻¹")

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

# Plot 3: 1100 cm⁻¹付近の拡大
ax3 = axes[1, 0]
mask_1100 = (freq_grid >= 1000) & (freq_grid <= 1200)
ax3.plot(freq_grid[mask_1100], J_total[mask_1100], 'b-', label='J_total', linewidth=2)
ax3.plot(freq_grid[mask_1100], J_linear[mask_1100], 'g--', label='Stage 1', linewidth=2)
ax3.plot(freq_grid[mask_1100], J_model_best[mask_1100], 'r-', label='Stage 1+2', linewidth=2)
ax3.axvline(1100, color='gray', linestyle=':', linewidth=1)
ax3.set_xlabel('Frequency (cm⁻¹)')
ax3.set_ylabel('J(ω) (cm⁻¹²)')
ax3.set_title('Zoom: 1100 cm⁻¹ region')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: 上位の4次相関の寄与（1050-1150領域に寄与するペア）
ax4 = axes[1, 1]
# 1100付近で寄与の大きいペアを選択
idx_1100 = np.argmin(np.abs(freq_grid - 1100))
pair_contrib_1100 = [(p, J_4th_contrib[p][idx_1100]) for p in pairs]
pair_contrib_1100_sorted = sorted(pair_contrib_1100, key=lambda x: abs(x[1]), reverse=True)
colors = plt.cm.tab20(np.linspace(0, 1, 10))
for rank, (p, _) in enumerate(pair_contrib_1100_sorted[:10]):
    ax4.plot(freq_grid, J_4th_contrib[p], color=colors[rank], 
             label=f'4th({p[0]},{p[1]})', linewidth=1, alpha=0.8)
ax4.plot(freq_grid, J_cross_best, 'k--', linewidth=2, label='Total')
ax4.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax4.set_xlabel('Frequency (cm⁻¹)')
ax4.set_ylabel('J(ω)')
ax4.set_title('Top 10 4th order contributions (by 1100 cm⁻¹ value)')
ax4.legend(fontsize=6)
ax4.grid(True, alpha=0.3)

# Plot 5: 1050-1150 cm⁻¹付近の拡大（和周波数領域）
ax5 = axes[2, 0]
mask_detail = (freq_grid >= 1050) & (freq_grid <= 1150)
ax5.plot(freq_grid[mask_detail], J_total[mask_detail], 'b-', label='J_total', linewidth=2)
ax5.plot(freq_grid[mask_detail], J_linear[mask_detail], 'g--', label='Stage 1', linewidth=2)
ax5.plot(freq_grid[mask_detail], J_model_best[mask_detail], 'r-', label='Stage 1+2', linewidth=2)
# 上位5ペアを個別表示
for rank, (p, _) in enumerate(pair_contrib_1100_sorted[:5]):
    ax5.plot(freq_grid[mask_detail], J_4th_contrib[p][mask_detail], 
             color=colors[rank], label=f'4th({p[0]},{p[1]})', linewidth=1, alpha=0.7)
ax5.axvline(1100, color='gray', linestyle=':', linewidth=1)
ax5.set_xlabel('Frequency (cm⁻¹)')
ax5.set_ylabel('J(ω) (cm⁻¹²)')
ax5.set_title('Zoom: 1050-1150 cm⁻¹ region')
ax5.legend(fontsize=6)
ax5.grid(True, alpha=0.3)

# Plot 6: Stage 1 残差と cross term の比較
ax6 = axes[2, 1]
ax6.plot(freq_grid, J_residual_1, 'g-', label='Stage 1 residual', linewidth=1.5)
ax6.plot(freq_grid, J_cross_best, 'r-', label='Cross term (4th)', linewidth=1.5)
ax6.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax6.set_xlabel('Frequency (cm$^{-1}$)')
ax6.set_ylabel('J($\\omega$)')
ax6.set_title('Stage 1 residual vs Cross term')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectral_fit_v7.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_fit_v7.png")
plt.close()

# ============================================================================
# 1100 cm⁻¹ における詳細分析
# ============================================================================
print("\n" + "="*70)
print("1100 cm⁻¹ における分析（4次相関の寄与）")
print("="*70)

idx_1100 = np.argmin(np.abs(freq_grid - 1100))
print(f"\n周波数: {freq_grid[idx_1100]:.1f} cm⁻¹")
print(f"  J_total:      {J_total[idx_1100]:.6e}")
print(f"  J_linear:     {J_linear[idx_1100]:.6e}")
print(f"  J_4th_total:  {J_cross_best[idx_1100]:.6e}")
print(f"  J_model:      {J_model_best[idx_1100]:.6e}")
print(f"  Residual:     {J_total[idx_1100] - J_model_best[idx_1100]:.6e}")

# 1100 cm⁻¹での4次相関寄与（上位10）
print(f"\n1100 cm⁻¹での4次相関寄与 (上位10):")
for rank, (p, val) in enumerate(pair_contrib_1100_sorted[:10], 1):
    freq_i = mode_freqs[mode_indices_map[p[0]]]
    freq_j = mode_freqs[mode_indices_map[p[1]]]
    print(f"  {rank:2d}. ({p[0]},{p[1]}): {val:+.6e}  (sum={freq_i+freq_j:.0f}cm⁻¹)")

# 1060, 1100, 1140 cm⁻¹でも確認
print("\n" + "-"*50)
for target_freq in [1060, 1080, 1100, 1120, 1140]:
    idx = np.argmin(np.abs(freq_grid - target_freq))
    print(f"\n{freq_grid[idx]:.0f} cm⁻¹:")
    print(f"  J_total = {J_total[idx]:.6e}")
    print(f"  J_model = {J_model_best[idx]:.6e}")
    print(f"  Explained = {100*J_model_best[idx]/J_total[idx]:.1f}%")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Stage 1 (Linear NNLS, all modes): R² = {R2_linear:.4f}")
print(f"Stage 2 (4th order, {n_pairs} pairs): R² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

print(f"\n上位非線形係数 d (d²):")
for rank, (k, (i, j), d, d2) in enumerate(d_squared_sorted[:15], 1):
    if d2 > 1e-8:
        freq_i = mode_freqs[mode_indices_map[i]]
        freq_j = mode_freqs[mode_indices_map[j]]
        print(f"  {rank:2d}. ({i},{j}): d={d:+.4e}, d²={d2:.4e}, νi+νj={freq_i+freq_j:.0f}cm⁻¹")

print(f"\n上位線形係数:")
sorted_idx_linear = np.argsort(c_k_physical)[::-1]
for rank, i in enumerate(sorted_idx_linear[:10], 1):
    if c_k_physical[i] > 1e-10:
        print(f"  {rank}. Mode {i+7}: c = {c_k_physical[i]:.6e}, freq = {mode_freqs[i]:.1f} cm⁻¹")

# 1100 cm⁻¹ 領域の説明力
print("\n" + "-"*50)
print("1050-1150 cm⁻¹ 領域の説明力:")
mask_region = (freq_grid >= 1050) & (freq_grid <= 1150)
J_total_region = np.sum(J_total[mask_region]**2)
J_residual_region = np.sum((J_total[mask_region] - J_model_best[mask_region])**2)
R2_region = 1 - J_residual_region / J_total_region
print(f"  領域R² = {R2_region:.4f}")
print(f"  J_total (region sum) = {np.sum(J_total[mask_region]):.6e}")
print(f"  J_model (region sum) = {np.sum(J_model_best[mask_region]):.6e}")

print("\n" + "="*70)
print("Complete!")
print("="*70)
