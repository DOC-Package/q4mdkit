"""
結晶構造の生成と管理
"""
import numpy as np
from typing import Tuple, List


class CrystalStructure:
    """結晶構造を表すクラス"""
    
    def __init__(self, lattice_constant: float, num_cells: Tuple[int, int, int]):
        """
        Args:
            lattice_constant: 格子定数 (nm)
            num_cells: 各方向のセル数 (nx, ny, nz)
        """
        self.lattice_constant = lattice_constant
        self.num_cells = num_cells
        self.positions = None
        
    def generate_fcc_lattice(self) -> np.ndarray:
        """面心立方格子(FCC)を生成"""
        a = self.lattice_constant
        nx, ny, nz = self.num_cells
        
        # FCC基本単位の相対座標
        base_positions = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5]
        ])
        
        positions = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    cell_origin = np.array([i, j, k]) * a
                    for base_pos in base_positions:
                        pos = cell_origin + base_pos * a
                        positions.append(pos)
        
        self.positions = np.array(positions)
        return self.positions
    
    def get_box_vectors(self) -> np.ndarray:
        """ボックスベクトルを取得"""
        a = self.lattice_constant
        nx, ny, nz = self.num_cells
        return np.array([
            [nx * a, 0, 0],
            [0, ny * a, 0],
            [0, 0, nz * a]
        ])
    
    def get_num_atoms(self) -> int:
        """原子数を取得"""
        return 4 * self.num_cells[0] * self.num_cells[1] * self.num_cells[2]
