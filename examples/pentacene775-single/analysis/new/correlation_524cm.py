#!/usr/bin/env python3
"""
524 cm⁻¹ピークに寄与する可能性のあるすべての3次・4次相関の解析

モード周波数の和・差のすべての組み合わせを考慮：
- 3次相関: ν_i, ν_j, ν_k, ν_i±ν_j, ν_i±ν_k, ν_j±ν_k
- 4次相関: さらに ν_l を含む組み合わせ
"""

import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement, product
import sys

# ターゲット周波数と許容幅
TARGET_FREQ = 524.0  # cm⁻¹
TOLERANCE = 15.0  # ±15 cm⁻¹

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
    
    def close(self):
        self.log.close()

# ログ設定
sys.stdout = Logger('correlation_524cm.log')

def load_mode_data():
    """モードデータの読み込み"""
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
    
    coords_file = 'mode_coords_new.txt'
    data = np.loadtxt(coords_file)
    mode_coords = data[:, 1:]  # 最初の列は時間
    
    return mode_coords, mode_info

def find_target_combinations(mode_info, freq_min=200, freq_max=1600):
    """
    ターゲット周波数(524 cm⁻¹)にピークを生じる可能性のあるモード組み合わせを探索
    
    3次相関 <q_i(t)q_j(t)q_k(0)> のFTで現れるピーク:
    - ν_k (直接)
    - |ν_i - ν_j| (差周波数)
    - ν_i + ν_j (和周波数、通常は弱い)
    
    4次相関 <q_i(t)q_j(t)q_k(0)q_l(0)> のFTで現れるピーク:
    - ν_k, ν_l (直接)
    - |ν_k - ν_l|, ν_k + ν_l
    - |ν_i - ν_j| (t側の差)
    """
    # 対象モードを抽出
    valid_modes = [(idx, freq) for idx, freq in mode_info.items() 
                   if freq_min <= freq <= freq_max]
    valid_modes.sort(key=lambda x: x[1])
    
    print(f"対象モード数: {len(valid_modes)}")
    print(f"周波数範囲: {valid_modes[0][1]:.1f} - {valid_modes[-1][1]:.1f} cm⁻¹\n")
    
    combinations_3rd = []  # (i, j, k, mechanism, expected_freq)
    combinations_4th = []  # (i, j, k, l, mechanism, expected_freq)
    
    # 3次相関の組み合わせ
    print("=" * 70)
    print(f"3次相関: {TARGET_FREQ:.0f} cm⁻¹ に寄与する組み合わせを探索中...")
    print("=" * 70)
    
    for i, freq_i in valid_modes:
        for j, freq_j in valid_modes:
            if j < i:
                continue
            for k, freq_k in valid_modes:
                mechanisms = []
                
                # ν_k が直接ターゲット
                if abs(freq_k - TARGET_FREQ) <= TOLERANCE:
                    mechanisms.append(('ν_k', freq_k))
                
                # |ν_i - ν_j| がターゲット (差周波数)
                diff_ij = abs(freq_i - freq_j)
                if abs(diff_ij - TARGET_FREQ) <= TOLERANCE:
                    mechanisms.append((f'|ν_i-ν_j|', diff_ij))
                
                # ν_i + ν_j がターゲット (和周波数)
                sum_ij = freq_i + freq_j
                if abs(sum_ij - TARGET_FREQ) <= TOLERANCE:
                    mechanisms.append((f'ν_i+ν_j', sum_ij))
                
                # |ν_i - ν_k| がターゲット
                diff_ik = abs(freq_i - freq_k)
                if abs(diff_ik - TARGET_FREQ) <= TOLERANCE and i != k:
                    mechanisms.append((f'|ν_i-ν_k|', diff_ik))
                
                # |ν_j - ν_k| がターゲット
                diff_jk = abs(freq_j - freq_k)
                if abs(diff_jk - TARGET_FREQ) <= TOLERANCE and j != k:
                    mechanisms.append((f'|ν_j-ν_k|', diff_jk))
                
                for mech, exp_freq in mechanisms:
                    combinations_3rd.append((i, j, k, mech, exp_freq, freq_i, freq_j, freq_k))
    
    print(f"3次相関の候補数: {len(combinations_3rd)}")
    
    # 4次相関の組み合わせ (計算コスト削減のため、主要なメカニズムに限定)
    print("\n" + "=" * 70)
    print(f"4次相関: {TARGET_FREQ:.0f} cm⁻¹ に寄与する組み合わせを探索中...")
    print("=" * 70)
    
    # 4次は組み合わせが膨大なので、主要パターンに限定
    # C4 = <q_i(t)q_j(t)q_k(0)q_l(0)>
    for i, freq_i in valid_modes:
        for j, freq_j in valid_modes:
            if j < i:
                continue
            for k, freq_k in valid_modes:
                for l, freq_l in valid_modes:
                    if l < k:
                        continue
                    
                    mechanisms = []
                    
                    # ν_k + ν_l がターゲット
                    sum_kl = freq_k + freq_l
                    if abs(sum_kl - TARGET_FREQ) <= TOLERANCE:
                        mechanisms.append(('ν_k+ν_l', sum_kl))
                    
                    # |ν_k - ν_l| がターゲット
                    diff_kl = abs(freq_k - freq_l)
                    if abs(diff_kl - TARGET_FREQ) <= TOLERANCE:
                        mechanisms.append(('|ν_k-ν_l|', diff_kl))
                    
                    # |ν_i - ν_j| がターゲット
                    diff_ij = abs(freq_i - freq_j)
                    if abs(diff_ij - TARGET_FREQ) <= TOLERANCE:
                        mechanisms.append(('|ν_i-ν_j|', diff_ij))
                    
                    # ν_i + ν_j がターゲット
                    sum_ij = freq_i + freq_j
                    if abs(sum_ij - TARGET_FREQ) <= TOLERANCE:
                        mechanisms.append(('ν_i+ν_j', sum_ij))
                    
                    for mech, exp_freq in mechanisms:
                        combinations_4th.append((i, j, k, l, mech, exp_freq, 
                                                  freq_i, freq_j, freq_k, freq_l))
    
    print(f"4次相関の候補数: {len(combinations_4th)}")
    
    return combinations_3rd, combinations_4th, valid_modes

def compute_static_3rd(mode_coords, i, j, k):
    """3次相関の静的値 <q_i q_j q_k> を計算（t=0）"""
    qi = mode_coords[:, i-7]
    qj = mode_coords[:, j-7]
    qk = mode_coords[:, k-7]
    return np.mean(qi * qj * qk)

def compute_static_4th(mode_coords, i, j, k, l):
    """4次相関の静的値 <q_i q_j q_k q_l> を計算（t=0）"""
    qi = mode_coords[:, i-7]
    qj = mode_coords[:, j-7]
    qk = mode_coords[:, k-7]
    ql = mode_coords[:, l-7]
    return np.mean(qi * qj * qk * ql)

def compute_correlation_3rd(mode_coords, i, j, k, max_lag=3000):
    """3次相関関数 C3(t) = <q_i(t)q_j(t)q_k(0)> を計算"""
    n_frames = mode_coords.shape[0]
    qi = mode_coords[:, i-7]  # インデックス調整
    qj = mode_coords[:, j-7]
    qk = mode_coords[:, k-7]
    
    C3 = np.zeros(max_lag)
    for t in range(max_lag):
        if t < n_frames:
            product = qi[t:] * qj[t:] * qk[:n_frames-t]
            C3[t] = np.mean(product)
    
    return C3

def compute_correlation_4th(mode_coords, i, j, k, l, max_lag=3000):
    """4次相関関数 C4(t) = <q_i(t)q_j(t)q_k(0)q_l(0)> を計算"""
    n_frames = mode_coords.shape[0]
    qi = mode_coords[:, i-7]
    qj = mode_coords[:, j-7]
    qk = mode_coords[:, k-7]
    ql = mode_coords[:, l-7]
    
    C4 = np.zeros(max_lag)
    for t in range(max_lag):
        if t < n_frames:
            product = qi[t:] * qj[t:] * qk[:n_frames-t] * ql[:n_frames-t]
            C4[t] = np.mean(product)
    
    return C4

def compute_spectrum(C, dt=4.0, freq_grid=None):
    """相関関数のフーリエ変換（コサイン変換）"""
    if freq_grid is None:
        freq_grid = np.linspace(0, 1500, 1501)
    
    n_points = len(C)
    t = np.arange(n_points) * dt  # fs
    
    # cm⁻¹への変換係数
    c_cm = 2.998e10  # cm/s
    fs_to_s = 1e-15
    
    spectrum = np.zeros(len(freq_grid))
    for idx, freq in enumerate(freq_grid):
        omega = 2 * np.pi * freq * c_cm * fs_to_s
        integrand = C * np.cos(omega * t)
        spectrum[idx] = simpson(integrand, x=t)
    
    return freq_grid, spectrum

def get_peak_at_target(freq_grid, spectrum, target=TARGET_FREQ, window=TOLERANCE):
    """ターゲット周波数付近のピーク高さを取得"""
    mask = (freq_grid >= target - window) & (freq_grid <= target + window)
    if not np.any(mask):
        return 0.0, target
    
    sub_spec = spectrum[mask]
    sub_freq = freq_grid[mask]
    
    # 最大絶対値を持つ点
    max_idx = np.argmax(np.abs(sub_spec))
    return sub_spec[max_idx], sub_freq[max_idx]

def main():
    print("=" * 70)
    print(f"524 cm⁻¹ ピークに寄与する相関関数の包括的解析")
    print(f"ターゲット: {TARGET_FREQ:.0f} ± {TOLERANCE:.0f} cm⁻¹")
    print("=" * 70)
    
    # データ読み込み
    mode_coords, mode_info = load_mode_data()
    print(f"\nモード座標: {mode_coords.shape[0]} frames × {mode_coords.shape[1]} modes")
    
    # 候補組み合わせの探索
    combs_3rd, combs_4th, valid_modes = find_target_combinations(mode_info)
    
    # パラメータ
    dt = 4.0  # fs
    max_lag = 3000
    freq_grid = np.linspace(400, 650, 251)  # 524付近に集中
    
    # ==== 3次相関: まずt=0の値でスクリーニング ====
    print("\n" + "=" * 70)
    print("3次相関関数: t=0でのスクリーニング")
    print("=" * 70)
    
    unique_3rd = {}  # 重複排除
    
    for idx, (i, j, k, mech, exp_freq, fi, fj, fk) in enumerate(combs_3rd):
        key = (i, j, k)
        if key in unique_3rd:
            unique_3rd[key]['mechanisms'].append((mech, exp_freq))
            continue
        
        unique_3rd[key] = {
            'modes': (i, j, k),
            'freqs': (fi, fj, fk),
            'mechanisms': [(mech, exp_freq)]
        }
    
    print(f"ユニークな3次相関の組み合わせ: {len(unique_3rd)}")
    
    # t=0の静的値を計算
    print("静的相関値 <q_i q_j q_k> を計算中...")
    static_3rd = []
    for key, info in unique_3rd.items():
        i, j, k = info['modes']
        c0 = compute_static_3rd(mode_coords, i, j, k)
        static_3rd.append({
            'key': key,
            'info': info,
            'c0': c0
        })
    
    # |<q_i q_j q_k>| でソートして上位200を選択
    static_3rd.sort(key=lambda x: abs(x['c0']), reverse=True)
    top_3rd = static_3rd[:200]
    
    print(f"上位200の |<q_i q_j q_k>| 範囲: {abs(top_3rd[0]['c0']):.4e} ~ {abs(top_3rd[-1]['c0']):.4e}")
    
    # 上位200についてスペクトル計算
    print("\n" + "=" * 70)
    print("3次相関関数: スペクトル計算 (上位200)")
    print("=" * 70)
    
    results_3rd = []
    for idx, item in enumerate(top_3rd):
        i, j, k = item['key']
        info = item['info']
        fi, fj, fk = info['freqs']
        
        if (idx + 1) % 20 == 0:
            print(f"  Computing {idx+1}/200...")
        
        C3 = compute_correlation_3rd(mode_coords, i, j, k, max_lag)
        freq, spec = compute_spectrum(C3, dt, freq_grid)
        peak_height, peak_freq = get_peak_at_target(freq, spec)
        
        results_3rd.append({
            'modes': (i, j, k),
            'freqs': (fi, fj, fk),
            'mechanisms': info['mechanisms'],
            'C3_0': item['c0'],
            'peak_height': peak_height,
            'peak_freq': peak_freq,
            'spectrum': spec.copy(),
            'correlation': C3.copy()
        })
    
    # ==== 4次相関: まずt=0の値でスクリーニング ====
    print("\n" + "=" * 70)
    print("4次相関関数: t=0でのスクリーニング")
    print("=" * 70)
    
    unique_4th = {}
    
    for idx, (i, j, k, l, mech, exp_freq, fi, fj, fk, fl) in enumerate(combs_4th):
        key = (i, j, k, l)
        if key in unique_4th:
            unique_4th[key]['mechanisms'].append((mech, exp_freq))
            continue
        
        unique_4th[key] = {
            'modes': (i, j, k, l),
            'freqs': (fi, fj, fk, fl),
            'mechanisms': [(mech, exp_freq)]
        }
    
    print(f"ユニークな4次相関の組み合わせ: {len(unique_4th)}")
    
    # t=0の静的値を計算
    print("静的相関値 <q_i q_j q_k q_l> を計算中...")
    static_4th = []
    for key, info in unique_4th.items():
        i, j, k, l = key
        c0 = compute_static_4th(mode_coords, i, j, k, l)
        static_4th.append({
            'key': key,
            'info': info,
            'c0': c0
        })
    
    # |<q_i q_j q_k q_l>| でソートして上位200を選択
    static_4th.sort(key=lambda x: abs(x['c0']), reverse=True)
    top_4th = static_4th[:200]
    
    print(f"上位200の |<q_i q_j q_k q_l>| 範囲: {abs(top_4th[0]['c0']):.4e} ~ {abs(top_4th[-1]['c0']):.4e}")
    
    # 上位200についてスペクトル計算
    print("\n" + "=" * 70)
    print("4次相関関数: スペクトル計算 (上位200)")
    print("=" * 70)
    
    results_4th = []
    for idx, item in enumerate(top_4th):
        i, j, k, l = item['key']
        info = item['info']
        fi, fj, fk, fl = info['freqs']
        
        if (idx + 1) % 20 == 0:
            print(f"  Computing {idx+1}/200...")
        
        C4 = compute_correlation_4th(mode_coords, i, j, k, l, max_lag)
        freq, spec = compute_spectrum(C4, dt, freq_grid)
        peak_height, peak_freq = get_peak_at_target(freq, spec)
        
        results_4th.append({
            'modes': (i, j, k, l),
            'freqs': (fi, fj, fk, fl),
            'mechanisms': info['mechanisms'],
            'C4_0': item['c0'],
            'peak_height': peak_height,
            'peak_freq': peak_freq,
            'spectrum': spec.copy(),
            'correlation': C4.copy()
        })
    
    # 結果のソート（ピーク高さの絶対値で）
    results_3rd.sort(key=lambda x: abs(x['peak_height']), reverse=True)
    results_4th.sort(key=lambda x: abs(x['peak_height']), reverse=True)
    
    # 結果の表示
    print("\n" + "=" * 70)
    print(f"3次相関: {TARGET_FREQ:.0f} cm⁻¹ 付近のピーク (上位30)")
    print("=" * 70)
    print(f"{'Rank':>4} {'Modes':>20} {'Freqs(cm⁻¹)':>25} {'Mechanism':>15} {'|C(0)|':>12} {'Peak Height':>12} {'Peak Freq':>10}")
    print("-" * 110)
    
    for rank, r in enumerate(results_3rd[:30], 1):
        i, j, k = r['modes']
        fi, fj, fk = r['freqs']
        mech_str = ','.join([m[0] for m in r['mechanisms']])
        print(f"{rank:>4} {f'{i}×{j}×{k}':>20} {f'{fi:.0f},{fj:.0f},{fk:.0f}':>25} {mech_str:>15} {abs(r['C3_0']):>12.4e} {r['peak_height']:>+12.4e} {r['peak_freq']:>10.1f}")
    
    print("\n" + "=" * 70)
    print(f"4次相関: {TARGET_FREQ:.0f} cm⁻¹ 付近のピーク (上位30)")
    print("=" * 70)
    print(f"{'Rank':>4} {'Modes':>25} {'Freqs(cm⁻¹)':>30} {'Mechanism':>15} {'|C(0)|':>12} {'Peak Height':>12} {'Peak Freq':>10}")
    print("-" * 120)
    
    for rank, r in enumerate(results_4th[:30], 1):
        i, j, k, l = r['modes']
        fi, fj, fk, fl = r['freqs']
        mech_str = ','.join([m[0] for m in r['mechanisms']])
        print(f"{rank:>4} {f'{i}×{j}×{k}×{l}':>25} {f'{fi:.0f},{fj:.0f},{fk:.0f},{fl:.0f}':>30} {mech_str:>15} {abs(r['C4_0']):>12.4e} {r['peak_height']:>+12.4e} {r['peak_freq']:>10.1f}")
    
    # === 可視化 ===
    print("\n" + "=" * 70)
    print("可視化を生成中...")
    print("=" * 70)
    
    # Figure 1: 3次と4次のピーク高さ比較（棒グラフ）
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # 3次相関
    ax = axes[0]
    n_show = min(25, len(results_3rd))
    heights_3rd = [r['peak_height'] for r in results_3rd[:n_show]]
    labels_3rd = [f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}" for r in results_3rd[:n_show]]
    colors_3rd = ['red' if h > 0 else 'blue' for h in heights_3rd]
    
    ax.bar(range(n_show), heights_3rd, color=colors_3rd, alpha=0.7)
    ax.set_xticks(range(n_show))
    ax.set_xticklabels(labels_3rd, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Peak Height at 524 cm⁻¹')
    ax.set_title(f'3rd Order Correlation: Peak Heights at {TARGET_FREQ:.0f} cm⁻¹')
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # 4次相関
    ax = axes[1]
    n_show = min(25, len(results_4th))
    heights_4th = [r['peak_height'] for r in results_4th[:n_show]]
    labels_4th = [f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}-{r['modes'][3]}" for r in results_4th[:n_show]]
    colors_4th = ['red' if h > 0 else 'blue' for h in heights_4th]
    
    ax.bar(range(n_show), heights_4th, color=colors_4th, alpha=0.7)
    ax.set_xticks(range(n_show))
    ax.set_xticklabels(labels_4th, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Peak Height at 524 cm⁻¹')
    ax.set_title(f'4th Order Correlation: Peak Heights at {TARGET_FREQ:.0f} cm⁻¹')
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_peaks.png', dpi=150)
    print("Saved: correlation_524cm_peaks.png")
    plt.close()
    
    # Figure 2: スペクトル比較（上位6つずつ）
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    
    for idx in range(6):
        ax = axes[0, idx] if idx < 3 else axes[1, idx-3]
        
        if idx < len(results_3rd):
            r = results_3rd[idx]
            ax.plot(freq_grid, r['spectrum'], 'b-', linewidth=1.5, label='3rd order')
        if idx < len(results_4th):
            r4 = results_4th[idx]
            ax2 = ax.twinx()
            ax2.plot(freq_grid, r4['spectrum'], 'r-', linewidth=1.5, label='4th order')
            ax2.set_ylabel('4th order', color='r')
        
        ax.axvline(x=TARGET_FREQ, color='g', linestyle='--', alpha=0.7, label=f'{TARGET_FREQ:.0f} cm⁻¹')
        ax.set_xlabel('Frequency (cm⁻¹)')
        ax.set_ylabel('3rd order', color='b')
        
        if idx < len(results_3rd):
            r = results_3rd[idx]
            modes_3 = f"3rd: {r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}"
        else:
            modes_3 = ""
        if idx < len(results_4th):
            r4 = results_4th[idx]
            modes_4 = f"4th: {r4['modes'][0]}-{r4['modes'][1]}-{r4['modes'][2]}-{r4['modes'][3]}"
        else:
            modes_4 = ""
        ax.set_title(f"Rank {idx+1}\n{modes_3}\n{modes_4}", fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_spectra.png', dpi=150)
    print("Saved: correlation_524cm_spectra.png")
    plt.close()
    
    # Figure 3: メカニズム別の寄与
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 3次相関のメカニズム別集計
    mech_3rd = {}
    for r in results_3rd:
        for mech, _ in r['mechanisms']:
            if mech not in mech_3rd:
                mech_3rd[mech] = []
            mech_3rd[mech].append(abs(r['peak_height']))
    
    ax = axes[0]
    mech_names = list(mech_3rd.keys())
    mech_max = [max(mech_3rd[m]) if mech_3rd[m] else 0 for m in mech_names]
    mech_mean = [np.mean(mech_3rd[m]) if mech_3rd[m] else 0 for m in mech_names]
    
    x = np.arange(len(mech_names))
    width = 0.35
    ax.bar(x - width/2, mech_max, width, label='Max', alpha=0.8)
    ax.bar(x + width/2, mech_mean, width, label='Mean', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mech_names, rotation=45, ha='right')
    ax.set_ylabel('|Peak Height|')
    ax.set_title('3rd Order: Contribution by Mechanism')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4次相関のメカニズム別集計
    mech_4th = {}
    for r in results_4th:
        for mech, _ in r['mechanisms']:
            if mech not in mech_4th:
                mech_4th[mech] = []
            mech_4th[mech].append(abs(r['peak_height']))
    
    ax = axes[1]
    mech_names = list(mech_4th.keys())
    mech_max = [max(mech_4th[m]) if mech_4th[m] else 0 for m in mech_names]
    mech_mean = [np.mean(mech_4th[m]) if mech_4th[m] else 0 for m in mech_names]
    
    x = np.arange(len(mech_names))
    ax.bar(x - width/2, mech_max, width, label='Max', alpha=0.8)
    ax.bar(x + width/2, mech_mean, width, label='Mean', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(mech_names, rotation=45, ha='right')
    ax.set_ylabel('|Peak Height|')
    ax.set_title('4th Order: Contribution by Mechanism')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_mechanisms.png', dpi=150)
    print("Saved: correlation_524cm_mechanisms.png")
    plt.close()
    
    # Figure 4: 3次 vs 4次 散布図
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 上位の結果をプロット
    peak_3rd = [abs(r['peak_height']) for r in results_3rd[:50]]
    peak_4th = [abs(r['peak_height']) for r in results_4th[:50]]
    
    max_val = max(max(peak_3rd) if peak_3rd else 0, max(peak_4th) if peak_4th else 0)
    
    ax.scatter(range(len(peak_3rd)), peak_3rd, s=80, c='blue', alpha=0.7, label='3rd order')
    ax.scatter(range(len(peak_4th)), peak_4th, s=80, c='red', alpha=0.7, label='4th order')
    ax.set_xlabel('Rank')
    ax.set_ylabel('|Peak Height| at 524 cm⁻¹')
    ax.set_title('3rd vs 4th Order Correlation: Peak Magnitude Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_3rd_vs_4th.png', dpi=150)
    print("Saved: correlation_524cm_3rd_vs_4th.png")
    plt.close()
    
    # Figure 5: 全スペクトルの重ね合わせ
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 3次相関（上位20）
    ax = axes[0]
    for idx, r in enumerate(results_3rd[:20]):
        alpha = 1.0 - idx * 0.04
        ax.plot(freq_grid, r['spectrum'], alpha=alpha, linewidth=1.0,
                label=f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}" if idx < 5 else None)
    ax.axvline(x=TARGET_FREQ, color='k', linestyle='--', linewidth=2, label=f'{TARGET_FREQ:.0f} cm⁻¹')
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Spectral Density')
    ax.set_title('3rd Order Correlations (Top 20)')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # 4次相関（上位20）
    ax = axes[1]
    for idx, r in enumerate(results_4th[:20]):
        alpha = 1.0 - idx * 0.04
        ax.plot(freq_grid, r['spectrum'], alpha=alpha, linewidth=1.0,
                label=f"{r['modes'][0]}-{r['modes'][1]}-{r['modes'][2]}-{r['modes'][3]}" if idx < 5 else None)
    ax.axvline(x=TARGET_FREQ, color='k', linestyle='--', linewidth=2, label=f'{TARGET_FREQ:.0f} cm⁻¹')
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Spectral Density')
    ax.set_title('4th Order Correlations (Top 20)')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_524cm_spectra_overlay.png', dpi=150)
    print("Saved: correlation_524cm_spectra_overlay.png")
    plt.close()
    
    # サマリー
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    
    if results_3rd:
        max_3rd = max(abs(r['peak_height']) for r in results_3rd)
        print(f"3次相関の最大ピーク高さ: {max_3rd:.4e}")
        top_3rd = results_3rd[0]
        print(f"  最大寄与: Mode {top_3rd['modes']} ({top_3rd['freqs'][0]:.0f}, {top_3rd['freqs'][1]:.0f}, {top_3rd['freqs'][2]:.0f} cm⁻¹)")
        print(f"  メカニズム: {top_3rd['mechanisms']}")
    
    if results_4th:
        max_4th = max(abs(r['peak_height']) for r in results_4th)
        print(f"\n4次相関の最大ピーク高さ: {max_4th:.4e}")
        top_4th = results_4th[0]
        print(f"  最大寄与: Mode {top_4th['modes']} ({top_4th['freqs'][0]:.0f}, {top_4th['freqs'][1]:.0f}, {top_4th['freqs'][2]:.0f}, {top_4th['freqs'][3]:.0f} cm⁻¹)")
        print(f"  メカニズム: {top_4th['mechanisms']}")
    
    if results_3rd and results_4th:
        ratio = max_4th / max_3rd if max_3rd > 0 else float('inf')
        print(f"\n4次/3次 比: {ratio:.4f}")
        if ratio > 1:
            print("→ 4次相関が支配的")
        else:
            print("→ 3次相関が支配的")
    
    print("\n" + "=" * 70)
    print("解析完了！")
    print(f"ログ保存先: correlation_524cm.log")
    print("=" * 70)

if __name__ == '__main__':
    main()
