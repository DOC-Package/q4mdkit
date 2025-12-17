"""
DFTB+ Result Reader Module (Unrestricted/Spin-Polarized Version)

This module provides functions to read matrices and eigenvectors from 
spin-polarized DFTB+ output files:

- Hamiltonian: hamsqr1.dat (α/up spin), hamsqr2.dat (β/down spin)
- Overlap: oversqr.dat (same for both spins)
- Eigenvectors: eigenvec.out (contains both up and down sections)
- Eigenvalues: band.out (contains SPIN 1 and SPIN 2 sections)
- Constraint potential: final_Vc.dat

File format for spin-polarized eigenvec.out:
    Coefficients and Mulliken populations of the atomic orbitals
    
    Eigenvector:   1    (up)
        1 C   s             0.001000   0.000002
              p_y           0.000111   0.000000
              ...
    Eigenvector:   2    (up)
        ...
    (after all up eigenvectors)
    Eigenvector:   1    (down)
        ...
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Dict, List, Optional
import re

def read_dftb_matrix(filepath: str | Path) -> Tuple[np.ndarray, int]:
    """
    Read a square matrix from DFTB+ output file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the matrix file (hamsqr1.dat, hamsqr2.dat, or oversqr.dat).
        
    Returns
    -------
    matrix : np.ndarray
        The square matrix (n_orbitals x n_orbitals).
    n_orbitals : int
        Number of orbitals (matrix dimension).
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Matrix file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Parse header: line 0 is comment, line 1 has values
    header_line = lines[1].strip()
    parts = header_line.split()
    
    if len(parts) < 3:
        raise ValueError(f"Invalid header format in {filepath}")
    
    is_real = parts[0].upper() == 'T'
    n_orbitals = int(parts[1])
    n_kpoints = int(parts[2])
    
    if not is_real:
        raise ValueError("Complex matrices are not supported yet")
    
    if n_kpoints != 1:
        raise ValueError(f"Only single k-point supported, got {n_kpoints}")
    
    # Skip header lines (5 lines total)
    data_start = 5
    
    # Read matrix data
    matrix = np.zeros((n_orbitals, n_orbitals), dtype=np.float64)
    
    for i in range(n_orbitals):
        line_idx = data_start + i
        if line_idx >= len(lines):
            raise ValueError(f"Unexpected end of file at row {i}")
        
        row_data = lines[line_idx].split()
        if len(row_data) != n_orbitals:
            raise ValueError(
                f"Row {i} has {len(row_data)} elements, expected {n_orbitals}"
            )
        
        matrix[i, :] = [float(x) for x in row_data]
    
    return matrix, n_orbitals


def read_hamiltonian_spin_polarized(
    dirpath: str | Path
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Read Hamiltonian matrices for both spins from a spin-polarized calculation.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing hamsqr1.dat and hamsqr2.dat.
        
    Returns
    -------
    H_alpha : np.ndarray
        α (up) spin Hamiltonian matrix in Hartree.
    H_beta : np.ndarray
        β (down) spin Hamiltonian matrix in Hartree.
    n_orbitals : int
        Number of orbitals.
    """
    dirpath = Path(dirpath)
    
    # hamsqr1.dat = α (up) spin
    filepath_alpha = dirpath / "hamsqr1.dat"
    H_alpha, n_orbitals = read_dftb_matrix(filepath_alpha)
    
    # hamsqr2.dat = β (down) spin
    filepath_beta = dirpath / "hamsqr2.dat"
    if filepath_beta.exists():
        H_beta, _ = read_dftb_matrix(filepath_beta)
    else:
        # Restricted calculation - use same Hamiltonian
        H_beta = H_alpha.copy()
    
    return H_alpha, H_beta, n_orbitals


def read_overlap(dirpath: str | Path) -> Tuple[np.ndarray, int]:
    """
    Read overlap matrix from a DFTB+ calculation directory.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing oversqr.dat.
        
    Returns
    -------
    overlap : np.ndarray
        The overlap matrix.
    n_orbitals : int
        Number of orbitals.
    """
    filepath = Path(dirpath) / "oversqr.dat"
    return read_dftb_matrix(filepath)


# =============================================================================
# Eigenvalue I/O (band.out) - Spin-Polarized
# =============================================================================

def read_eigenvalues_spin_polarized(
    filepath: str | Path
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read eigenvalues and occupations from spin-polarized DFTB+ band.out file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to band.out file.
        
    Returns
    -------
    eigenvalues_alpha : np.ndarray
        α (up) eigenvalues in eV.
    occupations_alpha : np.ndarray
        α occupation numbers.
    eigenvalues_beta : np.ndarray
        β (down) eigenvalues in eV.
    occupations_beta : np.ndarray
        β occupation numbers.
        
    Notes
    -----
    File format:
        KPT            1  SPIN            1  KWEIGHT    1.0000000000000000
             1    -21.3112  1.00000
             ...
        KPT            1  SPIN            2  KWEIGHT    1.0000000000000000
             1    -21.3110  1.00000
             ...
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Band file not found: {filepath}")
    
    eigenvalues_alpha = []
    occupations_alpha = []
    eigenvalues_beta = []
    occupations_beta = []
    
    current_spin = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('KPT'):
                # Parse spin from header: "KPT  1  SPIN  1  KWEIGHT ..."
                parts = line.split()
                try:
                    spin_idx = parts.index('SPIN') + 1
                    current_spin = int(parts[spin_idx])
                except (ValueError, IndexError):
                    current_spin = 1  # Default to spin 1
                continue
            
            parts = line.split()
            if len(parts) >= 3:
                try:
                    eig = float(parts[1])
                    occ = float(parts[2])
                    
                    if current_spin == 1:
                        eigenvalues_alpha.append(eig)
                        occupations_alpha.append(occ)
                    elif current_spin == 2:
                        eigenvalues_beta.append(eig)
                        occupations_beta.append(occ)
                except ValueError:
                    continue
    
    # If no SPIN 2 found, it's a restricted calculation
    if len(eigenvalues_beta) == 0:
        eigenvalues_beta = eigenvalues_alpha.copy()
        occupations_beta = occupations_alpha.copy()
    
    return (np.array(eigenvalues_alpha), np.array(occupations_alpha),
            np.array(eigenvalues_beta), np.array(occupations_beta))


def read_eigenvalues_from_dir_spin_polarized(
    dirpath: str | Path
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read eigenvalues from a spin-polarized DFTB+ calculation directory.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing band.out.
        
    Returns
    -------
    eigenvalues_alpha : np.ndarray
        α eigenvalues in eV.
    occupations_alpha : np.ndarray
        α occupation numbers.
    eigenvalues_beta : np.ndarray
        β eigenvalues in eV.
    occupations_beta : np.ndarray
        β occupation numbers.
    """
    filepath = Path(dirpath) / "band.out"
    return read_eigenvalues_spin_polarized(filepath)


# =============================================================================
# Eigenvector I/O from eigenvec.out - Spin-Polarized
# =============================================================================

def get_orbital_info_from_eigenvec(filepath: str | Path) -> List[Tuple[int, str, str]]:
    """
    Extract orbital information (atom index, element, orbital type) from eigenvec.out.
    
    Parameters
    ----------
    filepath : str or Path
        Path to eigenvec.out file.
        
    Returns
    -------
    orbital_info : list of tuples
        List of (atom_index, element, orbital_type) for each orbital.
        atom_index is 1-based.
    """
    filepath = Path(filepath)
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find first eigenvector section (up)
    start_line = None
    for i, line in enumerate(lines):
        if 'Eigenvector:' in line and '(up)' in line:
            start_line = i
            break
    
    if start_line is None:
        # Try without spin label (restricted calculation)
        for i, line in enumerate(lines):
            if 'Eigenvector:' in line:
                start_line = i
                break
    
    if start_line is None:
        raise ValueError("No eigenvector found in file")
    
    # Find next eigenvector section to determine range
    end_line = len(lines)
    for i in range(start_line + 1, len(lines)):
        if 'Eigenvector:' in lines[i]:
            end_line = i
            break
    
    # Parse orbital info
    orbital_info = []
    current_atom = None
    current_element = None
    
    for i in range(start_line + 1, end_line):
        line = lines[i]
        if not line.strip():
            continue
        
        parts = line.split()
        if len(parts) < 3:
            continue
        
        # Check if this line starts with atom index
        try:
            atom_idx = int(parts[0])
            current_atom = atom_idx
            current_element = parts[1]
            orbital_type = parts[2]
        except ValueError:
            # Line starts with orbital type (continuation of previous atom)
            orbital_type = parts[0]
        
        if current_atom is not None:
            orbital_info.append((current_atom, current_element, orbital_type))
    
    return orbital_info


def get_atom_orbital_map(
    orbital_info: List[Tuple[int, str, str]]
) -> Dict[int, List[int]]:
    """
    Create a mapping from atom index to orbital indices.
    
    Parameters
    ----------
    orbital_info : list of tuples
        List of (atom_index, element, orbital_type) for each orbital.
        
    Returns
    -------
    atom_to_orbitals : dict
        Dictionary mapping atom index (1-based) to list of orbital indices (0-based).
    """
    atom_to_orbitals = {}
    
    for orbital_idx, (atom_idx, element, orbital_type) in enumerate(orbital_info):
        if atom_idx not in atom_to_orbitals:
            atom_to_orbitals[atom_idx] = []
        atom_to_orbitals[atom_idx].append(orbital_idx)
    
    return atom_to_orbitals


def read_eigenvectors_from_eigenvec_out_spin_polarized(
    filepath: str | Path
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Read eigenvector coefficients from spin-polarized eigenvec.out file.
    
    Reads both (up) and (down) spin eigenvectors from the file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to eigenvec.out file.
        
    Returns
    -------
    C_alpha : np.ndarray
        α (up) eigenvector matrix (n_orbitals x n_eigvec).
        Each column is an eigenvector.
    C_beta : np.ndarray
        β (down) eigenvector matrix (n_orbitals x n_eigvec).
    n_orbitals : int
        Number of orbitals.
    n_eigvec : int
        Number of eigenvectors per spin.
        
    Notes
    -----
    File format:
        Eigenvector:   1    (up)
        
            1 C   s             0.001000   0.000002
                  p_y           0.000111   0.000000
                  p_z          -0.000025   0.000000
                  p_x          -0.000160   0.000000
        
            2 C   s             0.000596   0.000001
            ...
    
    The coefficient is the second to last column, Mulliken population is last.
    Lines starting with atom index have format: index element orbital coeff pop
    Continuation lines have format: orbital coeff pop (with leading spaces)
    """
    filepath = Path(filepath)
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find all eigenvector sections
    eigenvec_lines_up = []
    eigenvec_lines_down = []
    
    for i, line in enumerate(lines):
        if 'Eigenvector:' in line:
            if '(up)' in line:
                eigenvec_lines_up.append(i)
            elif '(down)' in line:
                eigenvec_lines_down.append(i)
            else:
                # Restricted calculation (no spin label)
                eigenvec_lines_up.append(i)
    
    n_eigvec = len(eigenvec_lines_up)
    
    if n_eigvec == 0:
        raise ValueError("No eigenvectors found in file")
    
    def parse_eigenvector(start_line: int, end_line: int) -> List[float]:
        """Parse one eigenvector section and return coefficients."""
        coeffs = []
        for i in range(start_line + 1, min(end_line, len(lines))):
            line = lines[i]
            
            # Skip empty lines
            if not line.strip():
                continue
            
            # Stop at next eigenvector
            if 'Eigenvector:' in line:
                break
            
            parts = line.split()
            if len(parts) < 2:
                continue
            
            # Check if this is a data line (coefficient and population at end)
            try:
                coeff = float(parts[-2])  # Coefficient is second to last
                float(parts[-1])  # Mulliken population is last
                coeffs.append(coeff)
            except (ValueError, IndexError):
                continue
        
        return coeffs
    
    # Parse first eigenvector to determine number of orbitals
    end_line = eigenvec_lines_up[1] if len(eigenvec_lines_up) > 1 else (
        eigenvec_lines_down[0] if len(eigenvec_lines_down) > 0 else len(lines)
    )
    first_coeffs = parse_eigenvector(eigenvec_lines_up[0], end_line)
    n_orbitals = len(first_coeffs)
    
    if n_orbitals == 0:
        raise ValueError("No orbitals found in eigenvector file")
    
    # Read all α (up) eigenvectors
    C_alpha = np.zeros((n_orbitals, n_eigvec), dtype=np.float64)
    
    for eig_idx, start_line in enumerate(eigenvec_lines_up):
        # Determine end line
        if eig_idx + 1 < len(eigenvec_lines_up):
            end_line = eigenvec_lines_up[eig_idx + 1]
        elif len(eigenvec_lines_down) > 0:
            end_line = eigenvec_lines_down[0]
        else:
            end_line = len(lines)
        
        coeffs = parse_eigenvector(start_line, end_line)
        if len(coeffs) == n_orbitals:
            C_alpha[:, eig_idx] = coeffs
        else:
            print(f"Warning: Eigenvector {eig_idx+1} (up) has {len(coeffs)} orbitals, expected {n_orbitals}")
    
    # Read all β (down) eigenvectors
    if len(eigenvec_lines_down) > 0:
        n_eigvec_down = len(eigenvec_lines_down)
        C_beta = np.zeros((n_orbitals, n_eigvec_down), dtype=np.float64)
        
        for eig_idx, start_line in enumerate(eigenvec_lines_down):
            # Determine end line
            if eig_idx + 1 < len(eigenvec_lines_down):
                end_line = eigenvec_lines_down[eig_idx + 1]
            else:
                end_line = len(lines)
            
            coeffs = parse_eigenvector(start_line, end_line)
            if len(coeffs) == n_orbitals:
                C_beta[:, eig_idx] = coeffs
            else:
                print(f"Warning: Eigenvector {eig_idx+1} (down) has {len(coeffs)} orbitals, expected {n_orbitals}")
    else:
        # Restricted calculation - same eigenvectors for both spins
        C_beta = C_alpha.copy()
    
    return C_alpha, C_beta, n_orbitals, n_eigvec


def read_eigenvectors_from_dir_spin_polarized(
    dirpath: str | Path,
    filename: str = "eigenvec.out"
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Read eigenvector matrices from a spin-polarized DFTB+ calculation directory.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing eigenvec.out.
    filename : str
        Name of the eigenvector file (default: eigenvec.out).
        
    Returns
    -------
    C_alpha : np.ndarray
        α eigenvector matrix (n_orbitals x n_eigvec).
    C_beta : np.ndarray
        β eigenvector matrix (n_orbitals x n_eigvec).
    n_orbitals : int
        Number of orbitals.
    n_eigvec : int
        Number of eigenvectors per spin.
    """
    filepath = Path(dirpath) / filename
    return read_eigenvectors_from_eigenvec_out_spin_polarized(filepath)


# =============================================================================
# Fragment Weight Matrix Construction
# =============================================================================

def build_fragment_weight_matrix(
    S_AO: np.ndarray,
    fragment_atoms: List[int],
    n_atoms: int,
    atom_to_orbitals: Optional[Dict[int, List[int]]] = None
) -> np.ndarray:
    """
    Build fragment weight matrix for Mulliken population analysis.
    
    w_A = P_A @ S_AO where P_A is the projection onto fragment A orbitals.
    
    For Mulliken partitioning: w_μν = δ_{μ∈A} S_μν
    
    Parameters
    ----------
    S_AO : np.ndarray
        AO overlap matrix.
    fragment_atoms : list of int
        List of atom indices (1-based) belonging to the fragment.
    n_atoms : int
        Total number of atoms.
    atom_to_orbitals : dict, optional
        Dictionary mapping atom index to orbital indices.
        
    Returns
    -------
    w : np.ndarray
        Fragment weight matrix.
    """
    n_orbitals = S_AO.shape[0]
    w = np.zeros((n_orbitals, n_orbitals), dtype=np.float64)
    
    if atom_to_orbitals is None:
        # Assume simple mapping (equal orbitals per atom)
        orbitals_per_atom = n_orbitals // n_atoms
        atom_to_orbitals = {}
        for atom_idx in range(1, n_atoms + 1):
            start = (atom_idx - 1) * orbitals_per_atom
            end = atom_idx * orbitals_per_atom
            atom_to_orbitals[atom_idx] = list(range(start, end))
    
    # Build weight matrix: w_μν = δ_{μ∈A} S_μν
    for atom_idx in fragment_atoms:
        if atom_idx in atom_to_orbitals:
            for orbital_idx in atom_to_orbitals[atom_idx]:
                w[orbital_idx, :] = S_AO[orbital_idx, :]
    
    return w


# =============================================================================
# Constraint Potential I/O
# =============================================================================

def read_constraint_potential(dirpath: str | Path) -> float:
    """
    Read constraint potential from final_Vc.dat.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing final_Vc.dat.
        
    Returns
    -------
    Vc : float
        Constraint potential in Hartree.
        
    Notes
    -----
    File format:
        # Final constraint potential Vc and deviation (N_calc - N_target)
        # Generated by DFTB+ with electronic constraints
        # Status: CONVERGED
        #
         # Index                 Vc [Ha]               Deviation
             1   1.061753957522534E-01   9.341008028229680E-05
    """
    filepath = Path(dirpath) / "final_Vc.dat"
    
    if not filepath.exists():
        raise FileNotFoundError(f"Constraint potential file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comment lines
            if line.startswith('#') or not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # Index, Vc, Deviation format
                    Vc = float(parts[1])
                    return Vc
                except (ValueError, IndexError):
                    continue
    
    raise ValueError("Could not parse constraint potential from final_Vc.dat")


# =============================================================================
# Energy I/O
# =============================================================================

def read_total_energy(dirpath: str | Path) -> float:
    """
    Read total energy from detailed.out.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to the directory containing detailed.out.
        
    Returns
    -------
    energy : float
        Total energy in Hartree.
    """
    filepath = Path(dirpath) / "detailed.out"
    
    if not filepath.exists():
        raise FileNotFoundError(f"Detailed output file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        for line in f:
            if 'Total energy' in line and 'H' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == 'H':
                        return float(parts[i - 1])
    
    raise ValueError("Total energy not found in detailed.out")


# =============================================================================
# Convenience function for loading all spin-polarized data
# =============================================================================

def load_spin_polarized_calculation(
    dirpath: str | Path
) -> dict:
    """
    Load all data from a spin-polarized DFTB+ calculation.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to DFTB+ output directory.
        
    Returns
    -------
    data : dict
        Dictionary containing:
        - 'H_alpha': α Hamiltonian matrix
        - 'H_beta': β Hamiltonian matrix
        - 'S': Overlap matrix
        - 'C_alpha': α eigenvector matrix
        - 'C_beta': β eigenvector matrix
        - 'eig_alpha': α eigenvalues
        - 'occ_alpha': α occupation numbers
        - 'eig_beta': β eigenvalues
        - 'occ_beta': β occupation numbers
        - 'n_alpha': Number of occupied α orbitals
        - 'n_beta': Number of occupied β orbitals
        - 'n_orbitals': Number of AOs
        - 'Vc': Constraint potential (if available)
        - 'energy': Total energy (if available)
    """
    dirpath = Path(dirpath)
    
    data = {}
    
    # Read Hamiltonians
    H_alpha, H_beta, n_orbitals = read_hamiltonian_spin_polarized(dirpath)
    data['H_alpha'] = H_alpha
    data['H_beta'] = H_beta
    data['n_orbitals'] = n_orbitals
    
    # Read overlap
    S, _ = read_overlap(dirpath)
    data['S'] = S
    
    # Read eigenvectors
    C_alpha, C_beta, _, _ = read_eigenvectors_from_dir_spin_polarized(dirpath)
    data['C_alpha'] = C_alpha
    data['C_beta'] = C_beta
    
    # Read eigenvalues and occupations
    eig_alpha, occ_alpha, eig_beta, occ_beta = read_eigenvalues_from_dir_spin_polarized(dirpath)
    data['eig_alpha'] = eig_alpha
    data['occ_alpha'] = occ_alpha
    data['eig_beta'] = eig_beta
    data['occ_beta'] = occ_beta
    
    # Count occupied orbitals
    data['n_alpha'] = np.sum(occ_alpha > 0.5).astype(int)
    data['n_beta'] = np.sum(occ_beta > 0.5).astype(int)
    
    # Read constraint potential (optional)
    try:
        data['Vc'] = read_constraint_potential(dirpath)
    except FileNotFoundError:
        data['Vc'] = None
    
    # Read total energy (optional)
    try:
        data['energy'] = read_total_energy(dirpath)
    except (FileNotFoundError, ValueError):
        data['energy'] = None
    
    return data


# =============================================================================
# Main: Test the module
# =============================================================================
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        dirpath = Path(sys.argv[1])
    else:
        dirpath = Path(__file__).parent / "output" / "frame_00000" / "fragment1"
    
    print("=" * 70)
    print("DFTB+ Result Reader Test (Spin-Polarized)")
    print("=" * 70)
    print(f"\nReading from: {dirpath}")
    
    if not dirpath.exists():
        print(f"Error: Directory not found: {dirpath}")
        sys.exit(1)
    
    # Load all data
    print("\nLoading spin-polarized calculation data...")
    data = load_spin_polarized_calculation(dirpath)
    
    print(f"\nResults:")
    print(f"  Number of AOs: {data['n_orbitals']}")
    print(f"  Occupied α orbitals: {data['n_alpha']}")
    print(f"  Occupied β orbitals: {data['n_beta']}")
    print(f"  Total electrons: {data['n_alpha'] + data['n_beta']}")
    
    print(f"\n  Hamiltonian matrices:")
    print(f"    H_α shape: {data['H_alpha'].shape}")
    print(f"    H_β shape: {data['H_beta'].shape}")
    print(f"    H_α[0,0] = {data['H_alpha'][0,0]:.6f} Ha")
    print(f"    H_β[0,0] = {data['H_beta'][0,0]:.6f} Ha")
    
    print(f"\n  Eigenvector matrices:")
    print(f"    C_α shape: {data['C_alpha'].shape}")
    print(f"    C_β shape: {data['C_beta'].shape}")
    
    print(f"\n  Eigenvalues (HOMO):")
    homo_idx_alpha = data['n_alpha'] - 1
    homo_idx_beta = data['n_beta'] - 1
    print(f"    ε_HOMO^α = {data['eig_alpha'][homo_idx_alpha]:.4f} eV")
    print(f"    ε_HOMO^β = {data['eig_beta'][homo_idx_beta]:.4f} eV")
    
    print(f"\n  Eigenvalues (LUMO):")
    lumo_idx_alpha = data['n_alpha']
    lumo_idx_beta = data['n_beta']
    print(f"    ε_LUMO^α = {data['eig_alpha'][lumo_idx_alpha]:.4f} eV")
    print(f"    ε_LUMO^β = {data['eig_beta'][lumo_idx_beta]:.4f} eV")
    
    if data['Vc'] is not None:
        print(f"\n  Constraint potential: V_c = {data['Vc']:.6f} Ha")
    
    if data['energy'] is not None:
        print(f"  Total energy: E = {data['energy']:.6f} Ha")
    
    # Check orthonormality of eigenvectors
    print("\n  Checking eigenvector orthonormality...")
    S = data['S']
    C_alpha = data['C_alpha']
    C_beta = data['C_beta']
    
    n_check = min(5, data['n_alpha'])
    
    # C^T S C should be identity for occupied orbitals
    identity_check_alpha = C_alpha[:, :n_check].T @ S @ C_alpha[:, :n_check]
    identity_check_beta = C_beta[:, :n_check].T @ S @ C_beta[:, :n_check]
    
    print(f"    α: ||C^T S C - I||_max = {np.max(np.abs(identity_check_alpha - np.eye(n_check))):.2e}")
    print(f"    β: ||C^T S C - I||_max = {np.max(np.abs(identity_check_beta - np.eye(n_check))):.2e}")
    
    # HOMO orbital info
    print(f"\n  HOMO orbital localization:")
    C_homo_alpha = C_alpha[:, homo_idx_alpha]
    C_homo_beta = C_beta[:, homo_idx_beta]
    
    # Get orbital info for atom mapping
    eigenvec_file = dirpath / "eigenvec.out"
    orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
    atom_to_orbitals = get_atom_orbital_map(orbital_info)
    n_atoms = len(atom_to_orbitals)
    n_atoms_half = n_atoms // 2
    
    # Calculate Mulliken population on each fragment
    pop_frag1_alpha = 0.0
    pop_frag2_alpha = 0.0
    
    for atom_idx in range(1, n_atoms_half + 1):
        for orb_idx in atom_to_orbitals[atom_idx]:
            for mu in range(data['n_orbitals']):
                pop_frag1_alpha += C_homo_alpha[orb_idx] * S[orb_idx, mu] * C_homo_alpha[mu]
    
    for atom_idx in range(n_atoms_half + 1, n_atoms + 1):
        for orb_idx in atom_to_orbitals[atom_idx]:
            for mu in range(data['n_orbitals']):
                pop_frag2_alpha += C_homo_alpha[orb_idx] * S[orb_idx, mu] * C_homo_alpha[mu]
    
    print(f"    HOMO^α on Fragment 1: {pop_frag1_alpha * 100:.1f}%")
    print(f"    HOMO^α on Fragment 2: {pop_frag2_alpha * 100:.1f}%")
