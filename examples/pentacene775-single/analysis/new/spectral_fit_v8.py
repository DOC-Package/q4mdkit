#!/usr/bin/env python3
"""
v8: Stage 1の残差をモード14, 15, 16, 40の3次・4次相関でフィット

Stage 1: Linear terms with NNLS (c_k ≥ 0) → g_k
Stage 2: モード14, 15, 16, 40の3次・4次相関で非線形最小二乗フィット
  - 3次相関: g_k * d_ij * <δ(q_i(t)q_j(t)) * δq_k(0)>
  - 4次相関: d_ij² * <δ(q_i(t)q_j(t)) * δ(q_i(0)q_j(0))>
  - g_k はStage 1で決定済み、フィットパラメータは d_ij のみ（10個）
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
print("v8: Stage 1残差をモード14,15,16,40の3次・4次相関でフィット")
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

# Cache for Stage 1 (v4と共通のキャッシュを使用)
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
# Stage 2: モード14,15,16,40の3次・4次相関をフィット
# ============================================================================
print("\n" + "="*70)
print("Stage 2: モード14,15,16,40の3次・4次相関をフィット")
print("="*70)

# モードのインデックス（配列インデックスはモード番号-7）
idx_14 = 14 - 7  # = 7
idx_15 = 15 - 7  # = 8
idx_16 = 16 - 7  # = 9
idx_40 = 40 - 7  # = 33

q14 = mode_coords[:, idx_14]
q15 = mode_coords[:, idx_15]
q16 = mode_coords[:, idx_16]
q40 = mode_coords[:, idx_40]

print(f"\nモード周波数:")
print(f"  Mode 14: {mode_freqs[idx_14]:.1f} cm⁻¹")
print(f"  Mode 15: {mode_freqs[idx_15]:.1f} cm⁻¹")
print(f"  Mode 16: {mode_freqs[idx_16]:.1f} cm⁻¹")
print(f"  Mode 40: {mode_freqs[idx_40]:.1f} cm⁻¹")

print(f"\n期待されるピーク位置:")
print(f"  (14,14): 2×{mode_freqs[idx_14]:.1f} = {2*mode_freqs[idx_14]:.1f} cm⁻¹")
print(f"  (14,15): {mode_freqs[idx_14]:.1f}+{mode_freqs[idx_15]:.1f} = {mode_freqs[idx_14]+mode_freqs[idx_15]:.1f} cm⁻¹")
print(f"  (14,16): {mode_freqs[idx_14]:.1f}+{mode_freqs[idx_16]:.1f} = {mode_freqs[idx_14]+mode_freqs[idx_16]:.1f} cm⁻¹")
print(f"  (14,40): {mode_freqs[idx_14]:.1f}+{mode_freqs[idx_40]:.1f} = {mode_freqs[idx_14]+mode_freqs[idx_40]:.1f} cm⁻¹")
print(f"  (15,15): 2×{mode_freqs[idx_15]:.1f} = {2*mode_freqs[idx_15]:.1f} cm⁻¹")
print(f"  (15,16): {mode_freqs[idx_15]:.1f}+{mode_freqs[idx_16]:.1f} = {mode_freqs[idx_15]+mode_freqs[idx_16]:.1f} cm⁻¹")
print(f"  (15,40): {mode_freqs[idx_15]:.1f}+{mode_freqs[idx_40]:.1f} = {mode_freqs[idx_15]+mode_freqs[idx_40]:.1f} cm⁻¹")
print(f"  (16,16): 2×{mode_freqs[idx_16]:.1f} = {2*mode_freqs[idx_16]:.1f} cm⁻¹")
print(f"  (16,40): {mode_freqs[idx_16]:.1f}+{mode_freqs[idx_40]:.1f} = {mode_freqs[idx_16]+mode_freqs[idx_40]:.1f} cm⁻¹")
print(f"  (40,40): 2×{mode_freqs[idx_40]:.1f} = {2*mode_freqs[idx_40]:.1f} cm⁻¹")

# ============================================================================
# 3次・4次相関のスペクトルを計算
# ============================================================================
print("\n3次・4次相関のスペクトルを計算中...")

cache_file_v8 = 'spectral_cache_v8_3rd_4th.npz'

# 4次相関ペアのリスト (4モードの場合は10ペア)
pairs_4th = [(14,14), (14,15), (14,16), (14,40), (15,15), (15,16), (15,40), (16,16), (16,40), (40,40)]
pair_names_4th = ['(14,14)', '(14,15)', '(14,16)', '(14,40)', '(15,15)', '(15,16)', '(15,40)', '(16,16)', '(16,40)', '(40,40)']

# 3次相関トリプレットのリスト (i<=j, k は任意)
# (i,j) ペアは i<=j で10通り、k は4通り → 40組
modes_list = [14, 15, 16, 40]
triples_3rd = []
triple_names_3rd = []
for i in modes_list:
    for j in modes_list:
        if j >= i:  # i <= j
            for k in modes_list:
                triples_3rd.append((i, j, k))
                triple_names_3rd.append(f'({i},{j},{k})')

print(f"  4次相関ペア: {len(pairs_4th)}組")
print(f"  3次相関トリプレット: {len(triples_3rd)}組")

q_dict = {14: q14, 15: q15, 16: q16, 40: q40}

if os.path.exists(cache_file_v8):
    print(f"Loading cached from {cache_file_v8}...")
    cache_v8 = np.load(cache_file_v8)
    if len(cache_v8['freq_grid']) != len(freq_grid):
        print("  Cache freq_grid mismatch, recomputing...")
        os.remove(cache_file_v8)
        cache_v8 = None
    else:
        J_4th = {}
        for i, j in pairs_4th:
            key = f'J_4th_{i}_{j}'
            J_4th[(i,j)] = cache_v8[key]
        J_3rd = {}
        for i, j, k in triples_3rd:
            key = f'J_3rd_{i}_{j}_{k}'
            J_3rd[(i,j,k)] = cache_v8[key]
else:
    cache_v8 = None

if cache_v8 is None:
    # 4次相関の計算
    J_4th = {}
    for i, j in pairs_4th:
        print(f"  Computing 4th({i},{j})...")
        C4 = compute_correlation_4th(q_dict[i], q_dict[j])
        J_4th[(i,j)] = compute_spectrum_from_corr(C4, freq_grid)
    
    # 3次相関の計算
    J_3rd = {}
    for idx, (i, j, k) in enumerate(triples_3rd):
        if (idx+1) % 10 == 0:
            print(f"  Computing 3rd: {idx+1}/{len(triples_3rd)}...")
        C3 = compute_correlation_3rd(q_dict[i], q_dict[j], q_dict[k])
        J_3rd[(i,j,k)] = compute_spectrum_from_corr(C3, freq_grid)
    
    save_dict = {'freq_grid': freq_grid}
    for i, j in pairs_4th:
        save_dict[f'J_4th_{i}_{j}'] = J_4th[(i,j)]
    for i, j, k in triples_3rd:
        save_dict[f'J_3rd_{i}_{j}_{k}'] = J_3rd[(i,j,k)]
    np.savez_compressed(cache_file_v8, **save_dict)
    print(f"  Saved to {cache_file_v8}")

print(f"\n各4次相関スペクトルの特性:")
for (i,j), name in zip(pairs_4th, pair_names_4th):
    J = J_4th[(i,j)]
    idx_peak = np.argmax(np.abs(J))
    print(f"  4th{name}: max={np.max(np.abs(J)):.6e}, peak={freq_grid[idx_peak]:.0f} cm⁻¹")

print(f"\n各3次相関スペクトルの特性 (上位10):")
J_3rd_max = [(t, np.max(np.abs(J_3rd[t]))) for t in triples_3rd]
J_3rd_max_sorted = sorted(J_3rd_max, key=lambda x: x[1], reverse=True)
for t, maxval in J_3rd_max_sorted[:10]:
    idx_peak = np.argmax(np.abs(J_3rd[t]))
    print(f"  3rd{t}: max={maxval:.6e}, peak={freq_grid[idx_peak]:.0f} cm⁻¹")

# ============================================================================
# Stage 2 フィッティング（3次: g_k * d_ij, 4次: d_ij²）
# J_cross = Σ g_k * d_ij * J_3rd(i,j,k) + Σ d_ij² * J_4th(i,j)
# g_k はStage 1で決定済み、フィットパラメータは d_ij のみ
# ============================================================================
print("\n" + "="*70)
print(f"Stage 2 フィッティング (d_ij: {len(pairs_4th)}個, g_k はStage 1から)")
print("="*70)

from scipy.optimize import least_squares

# g_k (Stage 1の線形係数) を取得
# c_k_physical[mode_idx] where mode_idx = mode_number - 7
g_k = {}
for m in modes_list:
    g_k[m] = c_k_physical[m - 7]
    print(f"  g_{m} = {g_k[m]:.6e}")

# J_4th をリストに変換
J_4th_list = [J_4th[p] for p in pairs_4th]

# 各(i,j)ペアに対して、Σ_k g_k * J_3rd(i,j,k) を事前計算
# これにより J_3rd_sum[(i,j)] = Σ_k g_k * J_3rd(i,j,k)
J_3rd_sum = {}
for i, j in pairs_4th:
    J_sum = np.zeros_like(freq_grid)
    for k in modes_list:
        J_sum += g_k[k] * J_3rd[(i,j,k)]
    J_3rd_sum[(i,j)] = J_sum

J_3rd_sum_list = [J_3rd_sum[p] for p in pairs_4th]

n_params = len(pairs_4th)  # d_ij のみ (10個)

def residual_func(d_params):
    J_model = np.zeros_like(freq_grid)
    for k, d in enumerate(d_params):
        # 3次相関の寄与: d_ij * Σ_k g_k * J_3rd(i,j,k)
        J_model += d * J_3rd_sum_list[k]
        # 4次相関の寄与: d_ij² * J_4th(i,j)
        J_model += d**2 * J_4th_list[k]
    return J_residual_1 - J_model

# 初期値
d0 = [0.01] * n_params
result = least_squares(residual_func, d0, method='lm')
d_best = result.x

print(f"\n非線形フィット結果:")
print(f"\n非線形係数 d_ij (d, d²):")
for k, name in enumerate(pair_names_4th):
    print(f"  d{name} = {d_best[k]:+.6e}, d² = {d_best[k]**2:.6e}")

# 各成分の寄与
J_3rd_contrib = {}  # (i,j) -> d_ij * Σ_k g_k * J_3rd(i,j,k)
J_3rd_total = np.zeros_like(freq_grid)
for k, (i,j) in enumerate(pairs_4th):
    J_3rd_contrib[(i,j)] = d_best[k] * J_3rd_sum[(i,j)]
    J_3rd_total += J_3rd_contrib[(i,j)]

J_4th_contrib = {}  # (i,j) -> d_ij² * J_4th(i,j)
J_4th_total = np.zeros_like(freq_grid)
for k, (i,j) in enumerate(pairs_4th):
    J_4th_contrib[(i,j)] = d_best[k]**2 * J_4th[(i,j)]
    J_4th_total += J_4th_contrib[(i,j)]

J_cross_best = J_3rd_total + J_4th_total

J_model_best = J_linear + J_cross_best
ss_res_best = np.sum((J_total - J_model_best)**2)
R2_best = 1 - ss_res_best / ss_tot

print(f"\nR² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")
print(f"  3次相関の寄与 ∫J_3rd = {np.trapezoid(J_3rd_total, freq_grid):.6e}")
print(f"  4次相関の寄与 ∫J_4th = {np.trapezoid(J_4th_total, freq_grid):.6e}")

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

# Plot 3: 524 cm⁻¹付近の拡大
ax3 = axes[1, 0]
mask_524 = (freq_grid >= 480) & (freq_grid <= 580)
ax3.plot(freq_grid[mask_524], J_total[mask_524], 'b-', label='J_total', linewidth=2)
ax3.plot(freq_grid[mask_524], J_linear[mask_524], 'g--', label='Stage 1', linewidth=2)
ax3.plot(freq_grid[mask_524], J_model_best[mask_524], 'r-', label='Stage 1+2', linewidth=2)
ax3.axvline(524, color='gray', linestyle=':', linewidth=1)
ax3.set_xlabel('Frequency (cm⁻¹)')
ax3.set_ylabel('J(ω) (cm⁻¹²)')
ax3.set_title('Zoom: 524 cm⁻¹ region')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: 3次・4次相関の寄与
ax4 = axes[1, 1]
ax4.plot(freq_grid, J_3rd_total, 'b-', label='3rd order total', linewidth=1.5, alpha=0.8)
ax4.plot(freq_grid, J_4th_total, 'r-', label='4th order total', linewidth=1.5, alpha=0.8)
ax4.plot(freq_grid, J_cross_best, 'k--', linewidth=2, label='3rd+4th Total')
ax4.axhline(0, color='gray', linestyle=':', linewidth=0.5)
ax4.set_xlabel('Frequency (cm⁻¹)')
ax4.set_ylabel('J(ω)')
ax4.set_title('3rd and 4th order contributions')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# Plot 5: 500 cm⁻¹付近の拡大（和周波数領域）
ax5 = axes[2, 0]
mask_500 = (freq_grid >= 450) & (freq_grid <= 600)
ax5.plot(freq_grid[mask_500], J_total[mask_500], 'b-', label='J_total', linewidth=2)
ax5.plot(freq_grid[mask_500], J_linear[mask_500], 'g--', label='Stage 1', linewidth=2)
ax5.plot(freq_grid[mask_500], J_model_best[mask_500], 'r-', label='Stage 1+2', linewidth=2)
ax5.plot(freq_grid[mask_500], J_3rd_total[mask_500], 'c-', label='3rd total', linewidth=1.5, alpha=0.7)
ax5.plot(freq_grid[mask_500], J_4th_total[mask_500], 'm-', label='4th total', linewidth=1.5, alpha=0.7)
ax5.axvline(524, color='gray', linestyle=':', linewidth=1)
ax5.set_xlabel('Frequency (cm⁻¹)')
ax5.set_ylabel('J(ω) (cm⁻¹²)')
ax5.set_title('Zoom: 500 cm⁻¹ region (sum frequencies)')
ax5.legend(fontsize=7)
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
plt.savefig('spectral_fit_v8.png', dpi=150, bbox_inches='tight')
print("Saved: spectral_fit_v8.png")
plt.close()

# ============================================================================
# 524 cm⁻¹ における分析
# ============================================================================
print("\n" + "="*70)
print("524 cm⁻¹ における分析（3次・4次相関の寄与）")
print("="*70)

idx_524 = np.argmin(np.abs(freq_grid - 524))
print(f"\n周波数: {freq_grid[idx_524]:.1f} cm⁻¹")
print(f"  J_total:      {J_total[idx_524]:.6e}")
print(f"  J_linear:     {J_linear[idx_524]:.6e}")
print(f"  J_3rd_total:  {J_3rd_total[idx_524]:.6e}")
print(f"  J_4th_total:  {J_4th_total[idx_524]:.6e}")
print(f"\n  各ペアの寄与 (3次: d*Σg_k*J_3rd, 4次: d²*J_4th):")
for (i,j), name in zip(pairs_4th, pair_names_4th):
    val_3rd = J_3rd_contrib[(i,j)][idx_524]
    val_4th = J_4th_contrib[(i,j)][idx_524]
    if np.abs(val_3rd) > 1e-10 or np.abs(val_4th) > 1e-10:
        print(f"    {name}: 3rd={val_3rd:+.6e}, 4th={val_4th:+.6e}")
print(f"\n  J_model:      {J_model_best[idx_524]:.6e}")
print(f"  Residual:     {J_total[idx_524] - J_model_best[idx_524]:.6e}")

# ============================================================================
# Summary
# ============================================================================
# ============================================================================
# 係数をファイル出力
# ============================================================================
print("\n" + "="*70)
print("係数をファイル出力")
print("="*70)

# 線形係数 g_k の出力
coeff_linear_file = 'coefficients_linear_v8.dat'
with open(coeff_linear_file, 'w') as f:
    f.write("# Linear coefficients g_k from Stage 1 (NNLS)\n")
    f.write("# Mode  Frequency(cm-1)  g_k\n")
    for i in range(N_modes):
        f.write(f"{i+7:4d}  {mode_freqs[i]:10.2f}  {c_k_physical[i]:16.8e}\n")
print(f"Saved: {coeff_linear_file}")

# 非線形係数 d_ij の出力
coeff_d_file = 'coefficients_d_v8.dat'
with open(coeff_d_file, 'w') as f:
    f.write("# Nonlinear coefficients d_ij from Stage 2\n")
    f.write("# 3rd order: g_k * d_ij * J_3rd(i,j,k)\n")
    f.write("# 4th order: d_ij^2 * J_4th(i,j)\n")
    f.write("# Pair(i,j)  d_ij  d_ij^2\n")
    for k, ((i,j), name) in enumerate(zip(pairs_4th, pair_names_4th)):
        f.write(f"{name:10s}  {d_best[k]:16.8e}  {d_best[k]**2:16.8e}\n")
print(f"Saved: {coeff_d_file}")

# サマリーファイルの出力
summary_file = 'fit_summary_v8.dat'
with open(summary_file, 'w') as f:
    f.write("# Spectral Fit v8 Summary\n")
    f.write(f"# Date: 2026-02-04\n")
    f.write(f"# Modes used for 3rd/4th order: 14, 15, 16, 40\n")
    f.write("#\n")
    f.write("# Model: J = sum_k g_k^2 J_k + sum_{i,j} [g_k d_ij J_3rd(i,j,k) + d_ij^2 J_4th(i,j)]\n")
    f.write("#\n")
    f.write(f"# Stage 1 (Linear NNLS -> g_k):\n")
    f.write(f"#   R² = {R2_linear:.6f}\n")
    f.write(f"#   ∫J_linear = {J_linear_integral:.6e} cm⁻¹²\n")
    f.write(f"#   Active modes: {n_active}/{N_modes}\n")
    f.write("#\n")
    f.write(f"# Stage 2 (Nonlinear fit -> d_ij, {len(pairs_4th)} parameters):\n")
    f.write(f"#   R² = {R2_best:.6f}\n")
    f.write(f"#   ΔR² = {R2_best - R2_linear:+.6f}\n")
    f.write(f"#   ∫J_3rd = {np.trapezoid(J_3rd_total, freq_grid):.6e}\n")
    f.write(f"#   ∫J_4th = {np.trapezoid(J_4th_total, freq_grid):.6e}\n")
    f.write("#\n")
    f.write("# g_k values for selected modes:\n")
    for m in modes_list:
        f.write(f"#   g_{m} = {g_k[m]:.6e}\n")
    f.write("#\n")
    f.write("# Top 10 linear coefficients (g_k):\n")
    sorted_idx_linear = np.argsort(c_k_physical)[::-1]
    for rank, i in enumerate(sorted_idx_linear[:10], 1):
        if c_k_physical[i] > 1e-10:
            f.write(f"#   {rank:2d}. Mode {i+7:3d}: g = {c_k_physical[i]:.6e}, freq = {mode_freqs[i]:.1f} cm⁻¹\n")
    f.write("#\n")
    f.write("# d_ij coefficients:\n")
    for k, name in enumerate(pair_names_4th):
        f.write(f"#   d{name} = {d_best[k]:+.6e} (d² = {d_best[k]**2:.6e})\n")
print(f"Saved: {summary_file}")

print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Stage 1 (Linear NNLS -> g_k): R² = {R2_linear:.4f}")
print(f"Stage 2 (d_ij: {len(pairs_4th)} params):  R² = {R2_best:.4f} (ΔR² = {R2_best - R2_linear:+.4f})")

print(f"\ng_k (選択モード):")
for m in modes_list:
    print(f"  g_{m} = {g_k[m]:.6e}")

print(f"\nd_ij 係数:")
for k, name in enumerate(pair_names_4th):
    print(f"  d{name} = {d_best[k]:+.6e} (d² = {d_best[k]**2:.6e})")

print(f"\n上位線形係数 g_k:")
sorted_idx_linear = np.argsort(c_k_physical)[::-1]
for rank, i in enumerate(sorted_idx_linear[:10], 1):
    if c_k_physical[i] > 1e-10:
        print(f"  {rank}. Mode {i+7}: g = {c_k_physical[i]:.6e}, freq = {mode_freqs[i]:.1f} cm⁻¹")

print("\n" + "="*70)
print("Complete!")
print("="*70)
