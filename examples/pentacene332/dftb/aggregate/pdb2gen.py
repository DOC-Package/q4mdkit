#!/usr/bin/env python3
"""PDBファイルからDFTB+用のGEN形式ファイルを生成"""

import sys
from pathlib import Path


def pdb_to_gen(pdb_file: str, gen_file: str, periodic: bool = False):
    """
    PDBファイルをDFTB+のGEN形式に変換

    Parameters
    ----------
    pdb_file : str
        入力PDBファイルのパス
    gen_file : str
        出力GENファイルのパス
    periodic : bool
        周期境界条件を使用するかどうか
    """
    atoms = []
    elements = []
    lattice = None

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                # 格子定数を読み取り（オプション）
                parts = line.split()
                a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
                alpha, beta, gamma = float(parts[4]), float(parts[5]), float(parts[6])
                # 簡易的な格子ベクトル計算（直方体近似）
                import math
                alpha_rad = math.radians(alpha)
                beta_rad = math.radians(beta)
                gamma_rad = math.radians(gamma)
                
                # 格子ベクトルの計算
                ax = a
                ay = 0.0
                az = 0.0
                
                bx = b * math.cos(gamma_rad)
                by = b * math.sin(gamma_rad)
                bz = 0.0
                
                cx = c * math.cos(beta_rad)
                cy = c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad)
                cz = math.sqrt(c**2 - cx**2 - cy**2)
                
                lattice = [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
                
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                # 原子情報を読み取り
                element = line[76:78].strip()
                if not element:
                    # 元素記号が無い場合、原子名から推定
                    atom_name = line[12:16].strip()
                    element = ''.join(c for c in atom_name if c.isalpha())[:2]
                    if len(element) > 1 and element[1].islower():
                        element = element[:2]
                    else:
                        element = element[0]
                
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                atoms.append((element, x, y, z))
                if element not in elements:
                    elements.append(element)

    n_atoms = len(atoms)
    
    with open(gen_file, 'w') as f:
        # ヘッダー行: 原子数と構造タイプ
        if periodic and lattice:
            f.write(f"{n_atoms} S\n")  # S = supercell (周期的)
        else:
            f.write(f"{n_atoms} C\n")  # C = cluster (非周期的)
        
        # 元素タイプのリスト
        f.write(" " + " ".join(elements) + "\n")
        
        # 原子座標
        for i, (elem, x, y, z) in enumerate(atoms, 1):
            elem_idx = elements.index(elem) + 1
            f.write(f"{i:5d} {elem_idx:3d}  {x:15.8f} {y:15.8f} {z:15.8f}\n")
        
        # 周期境界条件の場合、格子ベクトルを追加
        if periodic and lattice:
            f.write(f"  {0.0:15.8f} {0.0:15.8f} {0.0:15.8f}\n")  # origin
            for vec in lattice:
                f.write(f"  {vec[0]:15.8f} {vec[1]:15.8f} {vec[2]:15.8f}\n")

    print(f"変換完了: {pdb_file} -> {gen_file}")
    print(f"原子数: {n_atoms}")
    print(f"元素: {elements}")


if __name__ == "__main__":
    # デフォルトパス
    script_dir = Path(__file__).parent
    pdb_file = script_dir.parent / "pentacene.pdb"
    gen_file = script_dir / "pentacene.gen"
    
    # コマンドライン引数があれば使用
    if len(sys.argv) >= 2:
        pdb_file = sys.argv[1]
    if len(sys.argv) >= 3:
        gen_file = sys.argv[2]
    
    # 周期境界条件を使わない（クラスター計算）
    pdb_to_gen(str(pdb_file), str(gen_file), periodic=False)
