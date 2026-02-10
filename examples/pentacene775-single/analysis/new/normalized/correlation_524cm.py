#!/usr/bin/env python3
"""
524 cm⁻¹ピークに寄与する3次・4次相関の解析（正規化変位版）

正規化変位を使用:
    Δq̃_k(t) = Δq_k(t) / <Δq_k²>_t, where Δq_k(t) = q_k(t) - <q_k>_t

3次相関: C3 = <Δq̃_i(t)Δq̃_j(t)Δq̃_k(0)>
4次相関: C4 = <Δq̃_i(t)Δq̃_j(t)Δq̃_i(0)Δq̃_j(0)>  (同一モード対の自己相関のみ)

モード周波数の和・差のすべての組み合わせを考慮
上位50の組み合わせのみ解析
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import sys
import pickle
import os

# ターゲット周波数と許容幅
TARGET_FREQ = 524.0  # cm⁻¹
TOLERANCE = 40.0  # ±40 cm⁻¹ (486-564 cm⁻¹の範囲、モード14,14も含む)

# キャッシュファイル
CACHE_FILE = 'correlation_524cm_normalized_cache.pkl'

class Logger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w')
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = Logger('correlation_524cm_normalized.log')

def load_mode_data():
    """モードデータの読み込み（正規化変位を使用）"""
    # 周波数情報はmode_coords_new_stats.txtから
    stats_file = 'mode_coords_new_stats.txt'
    mode_info = {}
    with open(stats_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            mode_idx = int(parts[0])
            freq = float(parts[1])
            mode_info[mode_idx] = freq
    
    # 正規化変位を読み込む
    coords_file = 'normalized_displacement.txt'
    print(f"Loading normalized displacement from {coords_file}...")
    data = np.loadtxt(coords_file, comments='#')
    normalized_coords = data[:, 1:]  # 最初の列はフレームインデックス
    
    return normalized_coords, mode_info

def find_target_combinations(normalized_coords, mode_info, freq_min=200, freq_max=1600):
    """
    ターゲット周波数(524 cm⁻¹)にピークを生じる可能性のあるモード組み合わせを探索
    t=0の値でソートして上位を返す
    """
    valid_modes = [(idx, freq) for idx, freq in mode_info.items() 
                   if freq_min <= freq <= freq_max]
    valid_modes.sort(key=lambda x: x[1])
    
    print(f"対象モード数: {len(valid_modes)}")
    print(f"周波数範囲: {valid_modes[0][1]:.1f} - {valid_modes[-1][1]:.1f} cm⁻¹\n")
    
    # 3次相関の候補を収集（t=0の値で評価）
    print("3次相関の候補を収集中...")
    candidates_3rd = []
    
    for i, freq_i in valid_modes:
        qi = normalized_coords[:, i-7]
        for j, freq_j in valid_modes:
            if j < i:
                continue
            qj = normalized_coords[:, j-7]
            for k, freq_k in valid_modes:
                # 524 cm⁻¹に寄与する可能性があるか判定（和のみ、差は排除）
                mechanisms = []
                
                if abs(freq_k - TARGET_FREQ) <= TOLERANCE:
                    mechanisms.append(('ν_k', freq_k))
                sum_ij = freq_i + freq_j
                if abs(sum_ij - TARGET_FREQ) <= TOLERANCE:
                    mechanisms.append(('ν_i+ν_j', sum_ij))
                
                if mechanisms:
                    qk = normalized_coords[:, k-7]
                    C3_0 = np.mean(qi * qj * qk)
                    candidates_3rd.append({
                        'modes': (i, j, k),
                        'freqs': (freq_i, freq_j, freq_k),
                        'mechanisms': mechanisms,
                        'C_0': C3_0
                    })
    
    # |C(0)|でソート
    candidates_3rd.sort(key=lambda x: abs(x['C_0']), reverse=True)
    print(f"3次相関の候補数: {len(candidates_3rd)}, 上位50を選択")
    
    # 4次相関の候補を収集 (i,j)-(i,j)型のみ
    print("\n4次相関の候補を収集中 (<Δq̃_i(t)Δq̃_j(t)Δq̃_i(0)Δq̃_j(0)>型のみ)...")
    candidates_4th = []
    
    for i, freq_i in valid_modes:
        qi = normalized_coords[:, i-7]
        for j, freq_j in valid_modes:
            if j < i:
                continue
            qj = normalized_coords[:, j-7]
            
            # 524 cm⁻¹に寄与する可能性があるか判定（和のみ、差は排除）
            mechanisms = []
            
            sum_ij = freq_i + freq_j
            if abs(sum_ij - TARGET_FREQ) <= TOLERANCE:
                mechanisms.append(('ν_i+ν_j', sum_ij))
            
            if mechanisms:
                C4_0 = np.mean(qi * qj * qi * qj)
                candidates_4th.append({
                    'modes': (i, j),
                    'freqs': (freq_i, freq_j),
                    'mechanisms': mechanisms,
                    'C_0': C4_0
                })
    
    candidates_4th.sort(key=lambda x: abs(x['C_0']), reverse=True)
    print(f"4次相関の候補数: {len(candidates_4th)}, 上位50を選択")
    
    return candidates_3rd[:50], candidates_4th[:50]

def compute_correlation_3rd(normalized_coords, i, j, k, max_lag=3000):
    """
    3次相関関数 C3(t) = <Δq̃_i(t)Δq̃_j(t)Δq̃_k(0)>
    正規化変位はすでに平均が引かれているので、そのまま使用
    """
    n_frames = normalized_coords.shape[0]
    qi = normalized_coords[:, i-7]
    qj = normalized_coords[:, j-7]
    qk = normalized_coords[:, k-7]
    
    # 正規化変位の積
    qi_qj = qi * qj
    # 正規化変位の積の揺らぎ
    qi_qj_fluct = qi_qj - np.mean(qi_qj)
    
    C3 = np.zeros(max_lag)
    for t in range(max_lag):
        if t < n_frames:
            C3[t] = np.mean(qi_qj_fluct[t:] * qk[:n_frames-t])
    
    # Hanning窓を適用
    t_arr = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag))
    return C3 * window

def compute_correlation_4th(normalized_coords, i, j, max_lag=3000):
    """
    4次相関関数 C4(t) = <Δq̃_i(t)Δq̃_j(t)Δq̃_i(0)Δq̃_j(0)>
    正規化変位はすでに平均が引かれているので、そのまま使用
    """
    n_frames = normalized_coords.shape[0]
    qi = normalized_coords[:, i-7]
    qj = normalized_coords[:, j-7]
    
    # 正規化変位の積
    qi_qj = qi * qj
    # 揺らぎを計算
    qi_qj_fluct = qi_qj - np.mean(qi_qj)
    
    C4 = np.zeros(max_lag)
    for t in range(max_lag):
        if t < n_frames:
            C4[t] = np.mean(qi_qj_fluct[t:] * qi_qj_fluct[:n_frames-t])
    
    # Hanning窓を適用
    t_arr = np.arange(max_lag)
    window = 0.5 * (1 + np.cos(np.pi * t_arr / max_lag))
    return C4 * window

def compute_spectrum(C, dt=4.0, freq_grid=None):
    """相関関数のフーリエ変換"""
    if freq_grid is None:
        freq_grid = np.linspace(400, 650, 251)
    
    n_points = len(C)
    t = np.arange(n_points) * dt
    
    # 定数
    c_cm_fs = 2.99792458e-5  # cm/fs
    kB_cm = 0.695034800  # cm⁻¹/K
    T = 300.0  # K
    
    spectrum = np.zeros(len(freq_grid))
    for idx, freq in enumerate(freq_grid):
        integrand = C * np.cos(2 * np.pi * c_cm_fs * freq * t)
        spectrum[idx] = simpson(integrand, x=t)
    
    # prefactorを適用
    spectrum *= (2 * np.pi * c_cm_fs * freq_grid) / (kB_cm * T)
    
    return freq_grid, spectrum

def get_peak_at_target(freq_grid, spectrum, target=TARGET_FREQ, window=TOLERANCE):
    """ターゲット周波数付近のピーク高さを取得"""
    mask = (freq_grid >= target - window) & (freq_grid <= target + window)
    if not np.any(mask):
        return 0.0, target
    
    sub_spec = spectrum[mask]
    sub_freq = freq_grid[mask]
    max_idx = np.argmax(np.abs(sub_spec))
    return sub_spec[max_idx], sub_freq[max_idx]

def load_cache():
    """キャッシュを読み込む"""
    if os.path.exists(CACHE_FILE):
        print(f"キャッシュを読み込み中: {CACHE_FILE}")
        with open(CACHE_FILE, 'rb') as f:
            return pickle.load(f)
    return None

def save_cache(results_3rd, results_4th, freq_grid):
    """キャッシュを保存する"""
    cache = {
        'results_3rd': results_3rd,
        'results_4th': results_4th,
        'freq_grid': freq_grid,
        'target_freq': TARGET_FREQ,
        'tolerance': TOLERANCE
    }
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(cache, f)
    print(f"キャッシュを保存: {CACHE_FILE}")

def main():
    print("=" * 70)
    print(f"524 cm⁻¹ ピークに寄与する相関関数の解析（正規化変位版）")
    print(f"Using: Δq̃_k(t) = Δq_k(t) / <Δq_k²>_t")
    print(f"ターゲット: {TARGET_FREQ:.0f} ± {TOLERANCE:.0f} cm⁻¹")
    print("=" * 70)
    
    # キャッシュを確認
    cache = load_cache()
    
    if cache is not None:
        print("キャッシュから結果を復元しました")
        results_3rd = cache['results_3rd']
        results_4th = cache['results_4th']
        freq_grid = cache['freq_grid']
    else:
        normalized_coords, mode_info = load_mode_data()
        print(f"\n正規化変位: {normalized_coords.shape[0]} frames × {normalized_coords.shape[1]} modes")
        
        # 候補の選択
        candidates_3rd, candidates_4th = find_target_combinations(normalized_coords, mode_info)
        
        dt = 4.0
        max_lag = 3000
        freq_grid = np.linspace(400, 650, 251)
        
        # 3次相関の計算
        print("\n" + "=" * 70)
        print("3次相関関数の計算 (上位50)")
        print("=" * 70)
        
        results_3rd = []
        for idx, cand in enumerate(candidates_3rd):
            i, j, k = cand['modes']
            if (idx + 1) % 10 == 0:
                print(f"  Computing {idx+1}/50...")
            
            C3 = compute_correlation_3rd(normalized_coords, i, j, k, max_lag)
            freq, spec = compute_spectrum(C3, dt, freq_grid)
            peak_height, peak_freq = get_peak_at_target(freq, spec)
            
            results_3rd.append({
                **cand,
                'correlation': C3.copy(),
                'spectrum': spec.copy(),
                'peak_height': peak_height,
                'peak_freq': peak_freq
            })
        
        # 4次相関の計算
        print("\n" + "=" * 70)
        print("4次相関関数の計算 (上位50)")
        print("=" * 70)
        
        results_4th = []
        for idx, cand in enumerate(candidates_4th):
            i, j = cand['modes']
            if (idx + 1) % 10 == 0:
                print(f"  Computing {idx+1}/50...")
            
            C4 = compute_correlation_4th(normalized_coords, i, j, max_lag)
            freq, spec = compute_spectrum(C4, dt, freq_grid)
            peak_height, peak_freq = get_peak_at_target(freq, spec)
            
            results_4th.append({
                **cand,
                'correlation': C4.copy(),
                'spectrum': spec.copy(),
                'peak_height': peak_height,
                'peak_freq': peak_freq
            })
        
        # キャッシュを保存
        save_cache(results_3rd, results_4th, freq_grid)
    
    # ピーク高さでソート
    results_3rd.sort(key=lambda x: abs(x['peak_height']), reverse=True)
    results_4th.sort(key=lambda x: abs(x['peak_height']), reverse=True)
    
    # 結果表示
    print("\n" + "=" * 70)
    print(f"3次相関: {TARGET_FREQ:.0f} cm⁻¹ 付近のピーク (正規化変位)")
    print("=" * 70)
    print(f"{'Rank':>4} {'Modes':>15} {'Freqs(cm⁻¹)':>25} {'Mechanism':>15} {'|C(0)|':>12} {'Peak Height':>12}")
    print("-" * 90)
    
    for rank, r in enumerate(results_3rd, 1):
        i, j, k = r['modes']
        fi, fj, fk = r['freqs']
        mech_str = ','.join([m[0] for m in r['mechanisms']])
        print(f"{rank:>4} {f'{i}×{j}×{k}':>15} {f'{fi:.0f},{fj:.0f},{fk:.0f}':>25} {mech_str:>15} {abs(r['C_0']):>12.4e} {r['peak_height']:>+12.4e}")
    
    print("\n" + "=" * 70)
    print(f"4次相関: {TARGET_FREQ:.0f} cm⁻¹ 付近のピーク (正規化変位)")
    print("=" * 70)
    print(f"{'Rank':>4} {'Modes':>12} {'Freqs(cm⁻¹)':>20} {'Mechanism':>15} {'|C(0)|':>12} {'Peak Height':>12}")
    print("-" * 80)
    
    for rank, r in enumerate(results_4th, 1):
        i, j = r['modes']
        fi, fj = r['freqs']
        mech_str = ','.join([m[0] for m in r['mechanisms']])
        print(f"{rank:>4} {f'{i}×{j}':>12} {f'{fi:.0f},{fj:.0f}':>20} {mech_str:>15} {abs(r['C_0']):>12.4e} {r['peak_height']:>+12.4e}")
    
    # === 可視化 ===
    print("\n" + "=" * 70)
    print("可視化を生成中...")
    print("=" * 70)
    
    # Figure 1: ピーク高さ比較（棒グラフ）
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    ax = axes[0]
    heights_3rd = [r['peak_height'] for r in results_3rd]
    labels_3rd = [f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}" for r in results_3rd]
    colors_3rd = ['red' if h > 0 else 'blue' for h in heights_3rd]
    ax.bar(range(len(heights_3rd)), heights_3rd, color=colors_3rd, alpha=0.7)
    ax.set_xticks(range(len(heights_3rd)))
    ax.set_xticklabels(labels_3rd, rotation=90, ha='center', fontsize=7)
    ax.set_ylabel('Peak Height at 524 cm⁻¹')
    ax.set_title(f'3rd Order Correlation <Δq̃_i(t)Δq̃_j(t)Δq̃_k(0)>: Peak Heights at {TARGET_FREQ:.0f} cm⁻¹ (Normalized)')
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    ax = axes[1]
    heights_4th = [r['peak_height'] for r in results_4th]
    labels_4th = [f"{r['modes'][0]}-{r['modes'][1]}" for r in results_4th]
    colors_4th = ['red' if h > 0 else 'blue' for h in heights_4th]
    ax.bar(range(len(heights_4th)), heights_4th, color=colors_4th, alpha=0.7)
    ax.set_xticks(range(len(heights_4th)))
    ax.set_xticklabels(labels_4th, rotation=90, ha='center', fontsize=7)
    ax.set_ylabel('Peak Height at 524 cm⁻¹')
    ax.set_title(f'4th Order Correlation <Δq̃_i(t)Δq̃_j(t)Δq̃_i(0)Δq̃_j(0)>: Peak Heights at {TARGET_FREQ:.0f} cm⁻¹ (Normalized)')
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_normalized_peaks.png', dpi=150)
    print("Saved: correlation_524cm_normalized_peaks.png")
    plt.close()
    
    # Figure 2: スペクトル重ね合わせ
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax = axes[0]
    for idx, r in enumerate(results_3rd[:15]):
        alpha = 1.0 - idx * 0.05
        label = f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}" if idx < 6 else None
        ax.plot(freq_grid, r['spectrum'], alpha=alpha, linewidth=1.2, label=label)
    ax.axvline(x=TARGET_FREQ, color='k', linestyle='--', linewidth=2)
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Spectral Density')
    ax.set_title('3rd Order Correlations (Top 15, Normalized)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    ax = axes[1]
    for idx, r in enumerate(results_4th[:15]):
        alpha = 1.0 - idx * 0.05
        label = f"{r['modes'][0]}-{r['modes'][1]}" if idx < 6 else None
        ax.plot(freq_grid, r['spectrum'], alpha=alpha, linewidth=1.2, label=label)
    ax.axvline(x=TARGET_FREQ, color='k', linestyle='--', linewidth=2)
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Spectral Density')
    ax.set_title('4th Order Correlations (Top 15, Normalized)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_normalized_spectra.png', dpi=150)
    print("Saved: correlation_524cm_normalized_spectra.png")
    plt.close()
    
    # Figure 3: 3次 vs 4次 比較
    fig, ax = plt.subplots(figsize=(10, 8))
    
    peak_3rd = [abs(r['peak_height']) for r in results_3rd]
    peak_4th = [abs(r['peak_height']) for r in results_4th]
    
    ax.scatter(range(len(peak_3rd)), peak_3rd, s=80, c='blue', alpha=0.7, label='3rd order')
    ax.scatter(range(len(peak_4th)), peak_4th, s=80, c='red', alpha=0.7, label='4th order')
    ax.set_xlabel('Rank')
    ax.set_ylabel('|Peak Height| at 524 cm⁻¹')
    ax.set_title('3rd vs 4th Order Correlation: Peak Magnitude Comparison (Normalized)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_normalized_comparison.png', dpi=150)
    print("Saved: correlation_524cm_normalized_comparison.png")
    plt.close()
    
    # Figure 4: メカニズム別
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    mech_3rd = {}
    for r in results_3rd:
        for mech, _ in r['mechanisms']:
            if mech not in mech_3rd:
                mech_3rd[mech] = []
            mech_3rd[mech].append(abs(r['peak_height']))
    
    ax = axes[0]
    mech_names = list(mech_3rd.keys())
    mech_max = [max(mech_3rd[m]) for m in mech_names]
    mech_sum = [sum(mech_3rd[m]) for m in mech_names]
    x = np.arange(len(mech_names))
    ax.bar(x - 0.2, mech_max, 0.4, label='Max', alpha=0.8)
    ax.bar(x + 0.2, mech_sum, 0.4, label='Sum', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mech_names, rotation=45, ha='right')
    ax.set_ylabel('|Peak Height|')
    ax.set_title('3rd Order: By Mechanism (Normalized)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    mech_4th = {}
    for r in results_4th:
        for mech, _ in r['mechanisms']:
            if mech not in mech_4th:
                mech_4th[mech] = []
            mech_4th[mech].append(abs(r['peak_height']))
    
    ax = axes[1]
    mech_names = list(mech_4th.keys())
    mech_max = [max(mech_4th[m]) for m in mech_names]
    mech_sum = [sum(mech_4th[m]) for m in mech_names]
    x = np.arange(len(mech_names))
    ax.bar(x - 0.2, mech_max, 0.4, label='Max', alpha=0.8)
    ax.bar(x + 0.2, mech_sum, 0.4, label='Sum', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mech_names, rotation=45, ha='right')
    ax.set_ylabel('|Peak Height|')
    ax.set_title('4th Order: By Mechanism (Normalized)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_normalized_mechanisms.png', dpi=150)
    print("Saved: correlation_524cm_normalized_mechanisms.png")
    plt.close()
    
    # サマリー
    print("\n" + "=" * 70)
    print("Summary (Normalized Displacement)")
    print("=" * 70)
    
    max_3rd = max(abs(r['peak_height']) for r in results_3rd)
    max_4th = max(abs(r['peak_height']) for r in results_4th)
    
    print(f"3次相関の最大ピーク: {max_3rd:.4e}")
    top = results_3rd[0]
    print(f"  Mode {top['modes']} ({top['freqs'][0]:.0f}, {top['freqs'][1]:.0f}, {top['freqs'][2]:.0f} cm⁻¹)")
    print(f"  Mechanism: {top['mechanisms']}")
    
    print(f"\n4次相関の最大ピーク: {max_4th:.4e}")
    top = results_4th[0]
    print(f"  Mode {top['modes']} ({top['freqs'][0]:.0f}, {top['freqs'][1]:.0f} cm⁻¹)")
    print(f"  Mechanism: {top['mechanisms']}")
    
    ratio = max_4th / max_3rd if max_3rd > 0 else float('inf')
    print(f"\n4次/3次 比: {ratio:.4f}")
    if ratio > 1:
        print("→ 4次相関が支配的")
    else:
        print("→ 3次相関が支配的")
    
    print("\n" + "=" * 70)
    print("解析完了！")
    print("=" * 70)

if __name__ == '__main__':
    main()
