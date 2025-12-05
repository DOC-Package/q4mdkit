"""
Crystal structure generation and management
"""
import numpy as np
from typing import Tuple, List


class CrystalStructure:
    """Class representing a crystal structure."""
    
    def __init__(self, lattice_constant: float, num_cells: Tuple[int, int, int]):
        """
        Args:
            lattice_constant: Lattice constant (nm)
            num_cells: Number of cells in each direction (nx, ny, nz)
        """
        self.lattice_constant = lattice_constant
        self.num_cells = num_cells
        self.positions = None
        
    def generate_fcc_lattice(self) -> np.ndarray:
        """Generate a face-centered cubic (FCC) lattice."""
        a = self.lattice_constant
        nx, ny, nz = self.num_cells
        
    # Relative coordinates of the FCC basis
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
        """Return box vectors."""
        a = self.lattice_constant
        nx, ny, nz = self.num_cells
        return np.array([
            [nx * a, 0, 0],
            [0, ny * a, 0],
            [0, 0, nz * a]
        ])
    
    def get_num_atoms(self) -> int:
        """Return the number of atoms."""
        return 4 * self.num_cells[0] * self.num_cells[1] * self.num_cells[2]

class UnitCell:
    """Class representing a unit cell."""
    
    def __init__(self, lattice_constant: float):
        """
        Args:
            lattice_constant: Lattice constant (nm)
        """
        self.lattice_constant = lattice_constant
        
    def get_fcc_positions(self) -> np.ndarray:
        """Return atomic positions within an FCC unit cell."""
        a = self.lattice_constant
        return np.array([
            [0.0, 0.0, 0.0],
            [0.5 * a, 0.5 * a, 0.0],
            [0.5 * a, 0.0, 0.5 * a],
            [0.0, 0.5 * a, 0.5 * a]
        ])