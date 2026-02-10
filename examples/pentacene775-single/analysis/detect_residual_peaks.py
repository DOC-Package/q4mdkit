#!/usr/bin/env python3
"""
残差（Residual）のピーク位置を検出する

J_total(ω) - J_fit(ω) の残差からピークを検出し、
そのピーク位置（周波数）と強度を出力する
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*70)
print("Residual Peak Detection")
print("="*70)

# データ読み込み
print("\nLoading spectral density data...")

# spectral_fit_cross.pyの出力ファイルを読み込み
fit_data = np.loadtxt('spectral_density_quadratic_form.dat')
omega = fit_data[:, 0]
J_total = fit_data[:, 1]
J_fit = fit_data[:, 2]  # J_model
J_diag = fit_data[:, 3]
J_offdiag = fit_data[:, 4]
residual_saved = fit_data[:, 5]

print(f"  Loaded from: spectral_density_quadratic_form.dat")
print(f"  Columns: omega, J_total, J_model, J_diag, J_offdiag, residual")
N_omega = len(omega)
print(f"  Omega points: {N_omega}")

# 残差を計算
residual = J_total - J_fit
rel_residual = residual / (np.abs(J_total) + 1e-20)

print(f"\nSpectral range: {omega.min():.1f} - {omega.max():.1f} cm⁻¹")
print(f"Number of points: {len(omega)}")
print(f"\nJ_total statistics:")
print(f"  Mean: {np.mean(J_total):.6e}")
print(f"  Max:  {np.max(J_total):.6e} at {omega[np.argmax(J_total)]:.1f} cm⁻¹")
print(f"\nResidual statistics (absolute):")
print(f"  Mean: {np.mean(residual):.6e}")
print(f"  Std:  {np.std(residual):.6e}")
print(f"  Max:  {np.max(residual):.6e} at {omega[np.argmax(residual)]:.1f} cm⁻¹")
print(f"  Min:  {np.min(residual):.6e} at {omega[np.argmin(residual)]:.1f} cm⁻¹")
print(f"\nRelative residual statistics:")
print(f"  Mean: {np.mean(np.abs(rel_residual)):.4f}")
print(f"  Max:  {np.max(np.abs(rel_residual)):.4f}")

# ピーク検出パラメータ（相対残差ベース）
prominence_threshold_rel = 0.05  # 5%以上の変化
min_height_rel = 0.05  # 5%以上の残差

print(f"\nPeak detection parameters:")
print(f"  Relative prominence threshold: {prominence_threshold_rel*100:.1f}%")
print(f"  Minimum relative height: {min_height_rel*100:.1f}%")

# 正のピークを検出（相対残差で）
peaks_pos, properties_pos = find_peaks(rel_residual, 
                                       prominence=prominence_threshold_rel,
                                       height=min_height_rel,
                                       distance=5)  # 5ポイント以上離れたピーク

# 負のピークを検出（反転して検出）
peaks_neg, properties_neg = find_peaks(-rel_residual, 
                                       prominence=prominence_threshold_rel,
                                       height=min_height_rel,
                                       distance=5)

print("\n" + "="*70)
print("Positive Peaks (Underestimated by fit)")
print("="*70)
print(f"{'Peak #':<8} {'Frequency':<12} {'Residual':<15} {'Rel.Res.(%)':<15} {'Prominence':<15}")
print("-" * 80)

if len(peaks_pos) > 0:
    # 強度順にソート
    sorted_idx = np.argsort(residual[peaks_pos])[::-1]
    peaks_pos_sorted = peaks_pos[sorted_idx]
    
    for i, peak_idx in enumerate(peaks_pos_sorted):
        freq = omega[peak_idx]
        height = residual[peak_idx]
        rel_height = rel_residual[peak_idx] * 100
        prom = properties_pos['prominences'][sorted_idx[i]]
        print(f"{i+1:<8} {freq:<12.2f} {height:<15.6e} {rel_height:<15.2f} {prom:<15.4f}")
else:
    print("  No significant positive peaks detected")

print("\n" + "="*70)
print("Negative Peaks (Overestimated by fit)")
print("="*70)
print(f"{'Peak #':<8} {'Frequency':<12} {'Residual':<15} {'Rel.Res.(%)':<15} {'Prominence':<15}")
print("-" * 80)

if len(peaks_neg) > 0:
    # 強度順にソート（絶対値）
    sorted_idx = np.argsort(np.abs(residual[peaks_neg]))[::-1]
    peaks_neg_sorted = peaks_neg[sorted_idx]
    
    for i, peak_idx in enumerate(peaks_neg_sorted):
        freq = omega[peak_idx]
        height = residual[peak_idx]
        rel_height = rel_residual[peak_idx] * 100
        prom = properties_neg['prominences'][sorted_idx[i]]
        print(f"{i+1:<8} {freq:<12.2f} {height:<15.6e} {rel_height:<15.2f} {prom:<15.4f}")
else:
    print("  No significant negative peaks detected")

# ピーク周波数範囲の統計
print("\n" + "="*70)
print("Peak Frequency Statistics")
print("="*70)

if len(peaks_pos) > 0:
    peak_freqs_pos = omega[peaks_pos]
    print(f"\nPositive peaks:")
    print(f"  Count: {len(peaks_pos)}")
    print(f"  Frequency range: {peak_freqs_pos.min():.1f} - {peak_freqs_pos.max():.1f} cm⁻¹")
    print(f"  Mean frequency: {np.mean(peak_freqs_pos):.1f} cm⁻¹")

if len(peaks_neg) > 0:
    peak_freqs_neg = omega[peaks_neg]
    print(f"\nNegative peaks:")
    print(f"  Count: {len(peaks_neg)}")
    print(f"  Frequency range: {peak_freqs_neg.min():.1f} - {peak_freqs_neg.max():.1f} cm⁻¹")
    print(f"  Mean frequency: {np.mean(peak_freqs_neg):.1f} cm⁻¹")

# 可視化
print("\n可視化を作成中...")
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# 1. スペクトル密度の比較
ax = axes[0]
ax.plot(omega, J_total, 'k-', label='J_total (target)', linewidth=1.5, alpha=0.8)
ax.plot(omega, J_fit, 'r-', label='J_fit (quadratic form)', linewidth=1.5, alpha=0.8)
ax.set_xlabel('Frequency (cm⁻¹)')
ax.set_ylabel('J(ω)')
ax.set_title('Spectral Density Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. 残差
ax = axes[1]
ax.plot(omega, residual, 'b-', linewidth=1, alpha=0.7, label='Residual')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax.axhline(y=np.std(residual), color='gray', linestyle='--', alpha=0.3, label='±1σ')
ax.axhline(y=-np.std(residual), color='gray', linestyle='--', alpha=0.3)

# ピークをマーク
if len(peaks_pos) > 0:
    ax.plot(omega[peaks_pos], residual[peaks_pos], 'ro', markersize=6, 
            label=f'Positive peaks ({len(peaks_pos)})')
if len(peaks_neg) > 0:
    ax.plot(omega[peaks_neg], residual[peaks_neg], 'go', markersize=6, 
            label=f'Negative peaks ({len(peaks_neg)})')

ax.set_xlabel('Frequency (cm⁻¹)')
ax.set_ylabel('Residual')
ax.set_title('Residual = J_total - J_fit')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. 残差の拡大図（ピーク領域）
ax = axes[2]
ax.plot(omega, residual, 'b-', linewidth=1.5, alpha=0.7)
ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)

# 正のピークにラベル付け
if len(peaks_pos) > 0:
    top_pos_peaks = peaks_pos_sorted[:min(10, len(peaks_pos))]
    for peak_idx in top_pos_peaks:
        freq = omega[peak_idx]
        height = residual[peak_idx]
        ax.plot(freq, height, 'ro', markersize=8)
        ax.annotate(f'{freq:.0f}', xy=(freq, height), 
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, color='red')

# 負のピークにラベル付け
if len(peaks_neg) > 0:
    top_neg_peaks = peaks_neg_sorted[:min(10, len(peaks_neg))]
    for peak_idx in top_neg_peaks:
        freq = omega[peak_idx]
        height = residual[peak_idx]
        ax.plot(freq, height, 'go', markersize=8)
        ax.annotate(f'{freq:.0f}', xy=(freq, height), 
                   xytext=(5, -15), textcoords='offset points',
                   fontsize=8, color='green')

ax.set_xlabel('Frequency (cm⁻¹)')
ax.set_ylabel('Residual')
ax.set_title('Residual with Peak Labels')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('residual_peaks.png', dpi=150, bbox_inches='tight')
print(f"  Saved: residual_peaks.png")

# ピーク情報をファイルに保存
with open('residual_peaks.dat', 'w') as f:
    f.write("# Residual Peak Detection Results\n")
    f.write(f"# Relative prominence threshold: {prominence_threshold_rel*100:.1f}%\n")
    f.write(f"# Minimum relative height: {min_height_rel*100:.1f}%\n")
    f.write("#\n")
    f.write("# Positive Peaks (Underestimated by fit)\n")
    f.write("# Peak_No  Frequency(cm-1)  Residual  Prominence\n")
    
    if len(peaks_pos) > 0:
        sorted_idx_pos = np.argsort(residual[peaks_pos])[::-1]
        peaks_pos_sorted_save = peaks_pos[sorted_idx_pos]
        for i, peak_idx in enumerate(peaks_pos_sorted_save):
            freq = omega[peak_idx]
            height = residual[peak_idx]
            prom = properties_pos['prominences'][sorted_idx_pos[i]]
            f.write(f"{i+1:8d}  {freq:15.6f}  {height:15.6e}  {prom:15.6e}\n")
    
    f.write("#\n")
    f.write("# Negative Peaks (Overestimated by fit)\n")
    f.write("# Peak_No  Frequency(cm-1)  Residual  Prominence\n")
    
    if len(peaks_neg) > 0:
        sorted_idx_neg = np.argsort(np.abs(residual[peaks_neg]))[::-1]
        peaks_neg_sorted_save = peaks_neg[sorted_idx_neg]
        for i, peak_idx in enumerate(peaks_neg_sorted_save):
            freq = omega[peak_idx]
            height = residual[peak_idx]
            prom = properties_neg['prominences'][sorted_idx_neg[i]]
            f.write(f"{i+1:8d}  {freq:15.6f}  {height:15.6e}  {prom:15.6e}\n")

print(f"  Saved: residual_peaks.dat")

print("\n" + "="*70)
print("Complete!")
print("="*70)
