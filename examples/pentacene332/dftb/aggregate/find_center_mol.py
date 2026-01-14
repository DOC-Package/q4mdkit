#!/usr/bin/env python3
"""
PDBファイルから中心分子を特定し、原子インデックスを出力

36分子の集合体の場合、セルの中心に最も近い分子を特定する
"""

import numpy as np
from pathlib import Path


def find_center_molecule(pdb_file: str):
    """中心分子を特定"""
    molecules = {}  # mol_id -> [(element, x, y, z), ...]
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                mol_id = int(line[22:26])
                element = line[76:78].strip()
                if not element:
                    atom_name = line[12:16].strip()
                    element = ''.join(c for c in atom_name if c.isalpha())[:1]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                if mol_id not in molecules:
                    molecules[mol_id] = []
                molecules[mol_id].append((element, x, y, z))
    
    # 各分子の重心を計算
    centroids = {}
    for mol_id, atoms in molecules.items():
        coords = np.array([(x, y, z) for _, x, y, z in atoms])
        centroids[mol_id] = np.mean(coords, axis=0)
    
    # 全体の重心を計算
    all_centroids = np.array(list(centroids.values()))
    system_center = np.mean(all_centroids, axis=0)
    
    print(f"分子数: {len(molecules)}")
    print(f"各分子の原子数: {[len(atoms) for atoms in molecules.values()]}")
    print(f"システム中心: {system_center}")
    
    # 中心に最も近い分子を特定
    min_dist = float('inf')
    center_mol_id = None
    for mol_id, centroid in centroids.items():
        dist = np.linalg.norm(centroid - system_center)
        if dist < min_dist:
            min_dist = dist
            center_mol_id = mol_id
    
    print(f"\n中心分子: 分子 {center_mol_id}")
    print(f"中心分子の重心: {centroids[center_mol_id]}")
    print(f"システム中心からの距離: {min_dist:.3f} Å")
    
    # 原子インデックスを計算（1-based）
    atoms_per_mol = len(molecules[1])
    start_idx = (center_mol_id - 1) * atoms_per_mol + 1
    end_idx = center_mol_id * atoms_per_mol
    
    print(f"\n中心分子の原子インデックス: {start_idx} - {end_idx}")
    print(f"(DFTB+形式: {start_idx}:{end_idx})")
    
    # 各分子の情報を表示
    print("\n=== 全分子の重心座標 ===")
    for mol_id in sorted(molecules.keys()):
        centroid = centroids[mol_id]
        dist = np.linalg.norm(centroid - system_center)
        marker = " <-- CENTER" if mol_id == center_mol_id else ""
        print(f"分子 {mol_id:2d}: ({centroid[0]:7.2f}, {centroid[1]:7.2f}, {centroid[2]:7.2f}) 距離: {dist:6.2f}{marker}")
    
    return center_mol_id, start_idx, end_idx


if __name__ == "__main__":
    pdb_file = Path(__file__).parent.parent / "pentacene.pdb"
    find_center_molecule(str(pdb_file))
