"""
Tests for crystal structure module
"""
import sys
from pathlib import Path
import numpy as np

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from common.crystal_structure import CrystalStructure


def test_fcc_lattice_generation():
    """Test FCC lattice generation."""
    crystal = CrystalStructure(lattice_constant=0.4, num_cells=(2, 2, 2))
    positions = crystal.generate_fcc_lattice()
    
    # Check number of atoms
    expected_atoms = 4 * 2 * 2 * 2  # 32 atoms
    assert len(positions) == expected_atoms, f"Expected: {expected_atoms}, got: {len(positions)}"
    
    # Check coordinate shape
    assert positions.shape == (expected_atoms, 3), "Coordinate shape is invalid"
    
    print("✓ FCC lattice generation test passed")


def test_box_vectors():
    """Test box vectors."""
    crystal = CrystalStructure(lattice_constant=0.5, num_cells=(3, 3, 3))
    box_vectors = crystal.get_box_vectors()
    
    expected = np.array([
        [1.5, 0, 0],
        [0, 1.5, 0],
        [0, 0, 1.5]
    ])
    
    assert np.allclose(box_vectors, expected), "Box vectors are invalid"
    print("✓ Box vectors test passed")


def test_num_atoms():
    """Test atom count calculation."""
    crystal = CrystalStructure(lattice_constant=0.4, num_cells=(2, 3, 4))
    num_atoms = crystal.get_num_atoms()
    
    expected = 4 * 2 * 3 * 4  # 96 atoms
    assert num_atoms == expected, f"Expected: {expected}, got: {num_atoms}"
    print("✓ Atom count test passed")


if __name__ == "__main__":
    test_fcc_lattice_generation()
    test_box_vectors()
    test_num_atoms()
    print("\nAll tests passed!")
