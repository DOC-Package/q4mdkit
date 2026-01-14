#!/usr/bin/env python
"""
cdftbciファイルからJ, E1, E2を抽出し、
E1, E2についてはqm_energy.datのエネルギーを引いて揺らぎを計算する。
"""

import numpy as np
import glob
import os

def read_qm_energy(filepath):
    """qm_energy.datを読み込む"""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                time_fs = float(parts[1])
                energy_au = float(parts[2])
                data.append((time_fs, energy_au))
    return data

def read_cdftbci(filepath):
    """cdftbciファイルを読み込む"""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                frame = int(parts[0])
                time_fs = float(parts[1])
                j_lowdin = float(parts[2])  # meV
                e1 = float(parts[3])        # Ha
                e2 = float(parts[4])        # Ha
                data.append((frame, time_fs, j_lowdin, e1, e2))
    return data

def main():
    result_dir = os.path.dirname(os.path.abspath(__file__))
    
    # qm_energy.datを読み込む
    qm_energy_path = os.path.join(result_dir, 'qm_energy.dat')
    qm_data = read_qm_energy(qm_energy_path)
    print(f"Read {len(qm_data)} entries from qm_energy.dat")
    
    # qm_energyをtime_fsでインデックス化 (4fsステップ)
    qm_energy_dict = {time_fs: energy for time_fs, energy in qm_data}
    
    # cdftbciファイルを順番に読み込む
    # ファイル名のパターン: cdftbci{i}-{j}.dat (i=1-4, j=1-4)
    cdftbci_files = []
    for i in range(1, 5):
        for j in range(1, 5):
            filepath = os.path.join(result_dir, f'cdftbci{i}-{j}.dat')
            if os.path.exists(filepath):
                cdftbci_files.append(filepath)
    
    print(f"Found {len(cdftbci_files)} cdftbci files")
    
    # 全データを結合
    all_cdftbci_data = []
    global_frame = 0
    
    for filepath in cdftbci_files:
        file_data = read_cdftbci(filepath)
        print(f"Read {len(file_data)} entries from {os.path.basename(filepath)}")
        
        for local_frame, time_fs, j_lowdin, e1, e2 in file_data:
            # グローバルフレーム番号から時間を計算
            global_time_fs = global_frame * 4.0
            all_cdftbci_data.append((global_frame, global_time_fs, j_lowdin, e1, e2))
            global_frame += 1
    
    print(f"Total cdftbci entries: {len(all_cdftbci_data)}")
    
    # J, dE1, dE2 を計算して出力
    output_path = os.path.join(result_dir, 'cdftbci_extracted.dat')
    
    with open(output_path, 'w') as f:
        f.write("# CDFTB-CI Extracted Data with Energy Fluctuations\n")
        f.write("# Frame  Time(fs)   J_lowdin(meV)      dE1(Ha)           dE2(Ha)\n")
        f.write("# dE1 = E1(CDFTB-CI) - E(QM/MM)\n")
        f.write("# dE2 = E2(CDFTB-CI) - E(QM/MM)\n")
        
        matched_count = 0
        for frame, time_fs, j_lowdin, e1, e2 in all_cdftbci_data:
            if time_fs in qm_energy_dict:
                qm_energy = qm_energy_dict[time_fs]
                de1 = e1 - qm_energy  # 揺らぎ (Ha)
                de2 = e2 - qm_energy  # 揺らぎ (Ha)
                f.write(f"{frame:6d} {time_fs:10.2f}   {j_lowdin:12.4f}   {de1:16.10f}   {de2:16.10f}\n")
                matched_count += 1
            else:
                print(f"Warning: No QM energy found for time {time_fs} fs")
    
    print(f"Wrote {matched_count} entries to {output_path}")

if __name__ == "__main__":
    main()
