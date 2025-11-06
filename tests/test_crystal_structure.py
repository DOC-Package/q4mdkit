"""
結晶構造モジュールのテスト
"""
import sys
from pathlib import Path
import numpy as np

# srcディレクトリをパスに追加
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from modules.crystal_structure import CrystalStructure


def test_fcc_lattice_generation():
    """FCC格子生成のテスト"""
    crystal = CrystalStructure(lattice_constant=0.4, num_cells=(2, 2, 2))
    positions = crystal.generate_fcc_lattice()
    
    # 原子数の確認
    expected_atoms = 4 * 2 * 2 * 2  # 32原子
    assert len(positions) == expected_atoms, f"期待: {expected_atoms}, 実際: {len(positions)}"
    
    # 座標の形状確認
    assert positions.shape == (expected_atoms, 3), "座標の形状が不正"
    
    print("✓ FCC格子生成テスト合格")


def test_box_vectors():
    """ボックスベクトルのテスト"""
    crystal = CrystalStructure(lattice_constant=0.5, num_cells=(3, 3, 3))
    box_vectors = crystal.get_box_vectors()
    
    expected = np.array([
        [1.5, 0, 0],
        [0, 1.5, 0],
        [0, 0, 1.5]
    ])
    
    assert np.allclose(box_vectors, expected), "ボックスベクトルが不正"
    print("✓ ボックスベクトルテスト合格")


def test_num_atoms():
    """原子数計算のテスト"""
    crystal = CrystalStructure(lattice_constant=0.4, num_cells=(2, 3, 4))
    num_atoms = crystal.get_num_atoms()
    
    expected = 4 * 2 * 3 * 4  # 96原子
    assert num_atoms == expected, f"期待: {expected}, 実際: {num_atoms}"
    print("✓ 原子数計算テスト合格")


if __name__ == "__main__":
    test_fcc_lattice_generation()
    test_box_vectors()
    test_num_atoms()
    print("\n全テスト合格!")
