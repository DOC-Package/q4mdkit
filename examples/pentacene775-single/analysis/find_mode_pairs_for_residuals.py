#!/usr/bin/env python3
"""
残差ピークの周波数を生成できるモードペアを探す

3次・4次相関のスペクトル密度は、モード周波数の和・差の位置にピークを持つ：
- 3次相関: ω_k ± ω_l
- 4次相関: ω_k ± ω_l

ターゲット周波数: 524, 800, 1092, 1231, 1562 cm⁻¹
"""
import numpy as np

print("="*80)
print("Finding Mode Pairs for Residual Peak Frequencies")
print("="*80)

# ターゲット周波数（残差ピーク）
target_freqs = [524, 800, 1092, 1231, 1562]
tolerance = 30  # cm⁻¹の許容範囲

print(f"\nTarget frequencies: {target_freqs} cm⁻¹")
print(f"Tolerance: ±{tolerance} cm⁻¹")

# モード情報を読み込み
with open('g_k_coefficients_below_2000cm.dat') as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith('#') and 'R²' not in l]
    data = [l.split() for l in lines if len(l.split()) >= 3]

mode_indices = np.array([int(d[0]) for d in data])
mode_freqs = np.array([float(d[1]) for d in data])

# フィルタリング: 420 ≤ freq ≤ 1300 cm⁻¹
freq_mask = (mode_freqs >= 420) & (mode_freqs <= 1300)
mode_indices_filtered = mode_indices[freq_mask]
mode_freqs_filtered = mode_freqs[freq_mask]

N_modes = len(mode_freqs_filtered)
print(f"\nTotal modes in range (420-1300 cm⁻¹): {N_modes}")
print(f"Frequency range: {mode_freqs_filtered.min():.1f} - {mode_freqs_filtered.max():.1f} cm⁻¹")

# 各ターゲット周波数について、候補となるモードペアを探す
all_pairs = {}

for target_freq in target_freqs:
    print("\n" + "="*80)
    print(f"Target frequency: {target_freq} cm⁻¹ (±{tolerance} cm⁻¹)")
    print("="*80)
    
    pairs_sum = []  # ω_k + ω_l ≈ target
    pairs_diff = []  # |ω_k - ω_l| ≈ target
    
    for i in range(N_modes):
        for j in range(i+1, N_modes):
            freq_k = mode_freqs_filtered[i]
            freq_l = mode_freqs_filtered[j]
            idx_k = mode_indices_filtered[i]
            idx_l = mode_indices_filtered[j]
            
            # 和をチェック
            freq_sum = freq_k + freq_l
            if abs(freq_sum - target_freq) <= tolerance:
                error = freq_sum - target_freq
                pairs_sum.append((idx_k, idx_l, freq_k, freq_l, freq_sum, error))
            
            # 差をチェック
            freq_diff = abs(freq_k - freq_l)
            if abs(freq_diff - target_freq) <= tolerance:
                error = freq_diff - target_freq
                pairs_diff.append((idx_k, idx_l, freq_k, freq_l, freq_diff, error))
    
    all_pairs[target_freq] = {
        'sum': pairs_sum,
        'diff': pairs_diff
    }
    
    # 和のペア
    if pairs_sum:
        print(f"\n【Sum combinations: ω_k + ω_l ≈ {target_freq}】")
        print(f"Found {len(pairs_sum)} pairs")
        print(f"{'Mode k':<8} {'ω_k(cm⁻¹)':<12} {'Mode l':<8} {'ω_l(cm⁻¹)':<12} {'Sum':<10} {'Error':<10}")
        print("-" * 70)
        
        # エラーの小さい順にソート
        pairs_sum_sorted = sorted(pairs_sum, key=lambda x: abs(x[5]))
        for idx_k, idx_l, freq_k, freq_l, freq_sum, error in pairs_sum_sorted[:15]:
            print(f"{idx_k:<8} {freq_k:<12.2f} {idx_l:<8} {freq_l:<12.2f} {freq_sum:<10.2f} {error:>+9.2f}")
    else:
        print(f"\n【Sum combinations: ω_k + ω_l ≈ {target_freq}】")
        print("  No pairs found")
    
    # 差のペア
    if pairs_diff:
        print(f"\n【Difference combinations: |ω_k - ω_l| ≈ {target_freq}】")
        print(f"Found {len(pairs_diff)} pairs")
        print(f"{'Mode k':<8} {'ω_k(cm⁻¹)':<12} {'Mode l':<8} {'ω_l(cm⁻¹)':<12} {'Diff':<10} {'Error':<10}")
        print("-" * 70)
        
        # エラーの小さい順にソート
        pairs_diff_sorted = sorted(pairs_diff, key=lambda x: abs(x[5]))
        for idx_k, idx_l, freq_k, freq_l, freq_diff, error in pairs_diff_sorted[:15]:
            print(f"{idx_k:<8} {freq_k:<12.2f} {idx_l:<8} {freq_l:<12.2f} {freq_diff:<10.2f} {error:>+9.2f}")
    else:
        print(f"\n【Difference combinations: |ω_k - ω_l| ≈ {target_freq}】")
        print("  No pairs found")

# 統計情報
print("\n" + "="*80)
print("Summary Statistics")
print("="*80)

total_sum_pairs = sum(len(all_pairs[f]['sum']) for f in target_freqs)
total_diff_pairs = sum(len(all_pairs[f]['diff']) for f in target_freqs)

print(f"\nTotal candidate pairs:")
print(f"  Sum combinations (ω_k + ω_l): {total_sum_pairs}")
print(f"  Diff combinations (|ω_k - ω_l|): {total_diff_pairs}")
print(f"  Total unique pairs needed: {total_sum_pairs + total_diff_pairs}")

# 最も頻繁に現れるモードを特定
mode_count = {}
for target_freq in target_freqs:
    for pair_type in ['sum', 'diff']:
        for idx_k, idx_l, *_ in all_pairs[target_freq][pair_type]:
            mode_count[idx_k] = mode_count.get(idx_k, 0) + 1
            mode_count[idx_l] = mode_count.get(idx_l, 0) + 1

if mode_count:
    print("\n" + "="*80)
    print("Most Frequently Appearing Modes (Top 20)")
    print("="*80)
    print(f"{'Mode':<8} {'Frequency':<15} {'Appearances':<12}")
    print("-" * 40)
    
    sorted_modes = sorted(mode_count.items(), key=lambda x: x[1], reverse=True)
    for mode_idx, count in sorted_modes[:20]:
        # 周波数を取得
        mask = mode_indices_filtered == mode_idx
        if np.any(mask):
            freq = mode_freqs_filtered[mask][0]
            print(f"{mode_idx:<8} {freq:<15.2f} {count:<12}")

# ファイルに保存
print("\n" + "="*80)
print("Saving Results")
print("="*80)

with open('mode_pairs_for_residuals.dat', 'w') as f:
    f.write("# Mode pairs for residual peak frequencies\n")
    f.write(f"# Tolerance: ±{tolerance} cm⁻¹\n")
    f.write("# Mode range: 420-1300 cm⁻¹\n")
    f.write("#\n")
    
    for target_freq in target_freqs:
        f.write(f"\n# Target: {target_freq} cm⁻¹\n")
        f.write("#\n")
        
        # Sum pairs
        f.write(f"# Sum combinations (ω_k + ω_l ≈ {target_freq})\n")
        if all_pairs[target_freq]['sum']:
            f.write("# Mode_k  Freq_k  Mode_l  Freq_l  Sum  Error  Type\n")
            for idx_k, idx_l, freq_k, freq_l, freq_sum, error in all_pairs[target_freq]['sum']:
                f.write(f"{idx_k:6d}  {freq_k:8.2f}  {idx_l:6d}  {freq_l:8.2f}  {freq_sum:8.2f}  {error:+7.2f}  sum\n")
        else:
            f.write("# No pairs found\n")
        
        f.write("#\n")
        
        # Diff pairs
        f.write(f"# Difference combinations (|ω_k - ω_l| ≈ {target_freq})\n")
        if all_pairs[target_freq]['diff']:
            f.write("# Mode_k  Freq_k  Mode_l  Freq_l  Diff  Error  Type\n")
            for idx_k, idx_l, freq_k, freq_l, freq_diff, error in all_pairs[target_freq]['diff']:
                f.write(f"{idx_k:6d}  {freq_k:8.2f}  {idx_l:6d}  {freq_l:8.2f}  {freq_diff:8.2f}  {error:+7.2f}  diff\n")
        else:
            f.write("# No pairs found\n")
        
        f.write("#\n")

print("  Saved: mode_pairs_for_residuals.dat")

print("\n" + "="*80)
print("Complete!")
print("="*80)
print(f"\nNext steps:")
print(f"  1. Use these mode pairs to compute 3rd/4th order correlation functions")
print(f"  2. Fit the residual with these higher-order terms")
print(f"  3. Expected parameters: ~{total_sum_pairs + total_diff_pairs} pairs")
