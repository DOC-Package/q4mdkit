"""
CDFTB Hamiltonian Verification Module

This module provides functions to verify that the CDFTB Hamiltonian
H' = H + Vc * wc satisfies the eigenvalue equation H'C = SCE.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict

from cdftbci import (
    CDFTBResults,
    read_fragment_matrices,
)
from cdftb_result_reader import (
    build_fragment_weight_matrix,
    get_orbital_info_from_eigenvec,
    get_atom_orbital_map,
    read_eigenvalues_from_dir,
)


def verify_cdftb_hamiltonian(
    results: CDFTBResults,
    constrained_atoms: list,
    eigenvec_file: str | Path,
    tolerance: float = 2e-4,
    verbose: bool = True
) -> bool:
    """
    Verify that H'C = SCE holds for CDFTB Hamiltonian H' = H + Vc * w.
    
    Parameters
    ----------
    results : CDFTBResults
        Container with H, S, eigenvectors, eigenvalues, and Vc.
    constrained_atoms : list
        List of atom indices (1-based) that are constrained.
    eigenvec_file : str or Path
        Path to eigenvec.out for orbital info.
    tolerance : float
        Maximum allowed residual for verification to pass.
    verbose : bool
        If True, print verification details.
        
    Returns
    -------
    is_valid : bool
        True if H'C = SCE is satisfied within tolerance.
    """
    H = results.hamiltonian
    S = results.overlap
    C = results.eigenvectors
    eigenvalues = results.eigenvalues
    Vc = results.Vc
    
    # Check required fields
    if C is None:
        if verbose:
            print("Error: eigenvectors not loaded")
        return False
    
    if eigenvalues is None:
        if verbose:
            print("Error: eigenvalues not loaded")
        return False
    
    if Vc is None:
        if verbose:
            print("Error: Vc not loaded")
        return False
    
    # Get orbital info for weight matrix
    eigenvec_file = Path(eigenvec_file)
    orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
    atom_to_orbitals = get_atom_orbital_map(orbital_info)
    n_atoms = len(atom_to_orbitals)
    
    # Build weight matrix
    w = build_fragment_weight_matrix(S, constrained_atoms, n_atoms,
                                     atom_to_orbitals=atom_to_orbitals)
    
    # Build effective Hamiltonian: H' = H + Vc * w
    H_prime = H + Vc * w
    
    n = H.shape[0]
    
    # Verify H'C = SCE using DFTB+ eigenvalues
    # Compute residuals: ||H'C_i - e_i * SC_i||
    residuals = np.zeros(n)
    for i in range(n):
        c_i = C[:, i]
        e_i = eigenvalues[i]
        residual = np.linalg.norm(H_prime @ c_i - e_i * S @ c_i)
        residuals[i] = residual
    
    max_residual = np.max(residuals)
    mean_residual = np.mean(residuals)
    
    # Check orthonormality of eigenvectors: C^T S C should be identity
    CtSC = C.T @ S @ C
    orthogonality_error = np.max(np.abs(CtSC - np.eye(n)))
    
    is_valid = max_residual < tolerance
    
    if verbose:
        print(f"CDFTB Hamiltonian Verification (H' = H + Vc*w)")
        print(f"  Vc = {Vc:.6f} Ha")
        print(f"  Constrained atoms: {constrained_atoms[0]}-{constrained_atoms[-1]}")
        print(f"  Eigenvalues from band.out (converted to Ha)")
        print(f"  Residual ||H'C - SCe||:")
        print(f"    Max: {max_residual:.6e}")
        print(f"    Mean: {mean_residual:.6e}")
        print(f"  Orthonormality error: {orthogonality_error:.6e}")
        print(f"  Tolerance: {tolerance:.6e}")
        print(f"  Result: {'PASS' if is_valid else 'FAIL'}")
    
    return is_valid


def verify_cdftb_from_dir(
    output_dir: str | Path,
    frame_id: int,
    fragment_id: int,
    tolerance: float = 2e-4,
    verbose: bool = True
) -> bool:
    """
    Verify CDFTB Hamiltonian for a specific frame and fragment.
    
    Parameters
    ----------
    output_dir : str or Path
        Base output directory.
    frame_id : int
        Frame number.
    fragment_id : int
        Fragment number (1 or 2).
    tolerance : float
        Maximum allowed residual for verification to pass.
    verbose : bool
        If True, print verification details.
        
    Returns
    -------
    is_valid : bool
        True if H'C = SCE is satisfied within tolerance.
    """
    output_dir = Path(output_dir)
    fragment_dir = output_dir / f"frame_{frame_id:05d}" / f"fragment{fragment_id}"
    
    # Read results
    results = read_fragment_matrices(output_dir, frame_id, fragment_id)
    
    # Get orbital info
    eigenvec_file = fragment_dir / "eigenvec.out"
    orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
    atom_to_orbitals = get_atom_orbital_map(orbital_info)
    n_atoms = len(atom_to_orbitals)
    n_atoms_half = n_atoms // 2  # 36 atoms per pentacene
    
    # Determine constrained atoms based on fragment_id
    # fragment 1: atoms 1-36, fragment 2: atoms 37-72
    if fragment_id == 1:
        constrained_atoms = list(range(1, n_atoms_half + 1))
    else:
        constrained_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
    
    return verify_cdftb_hamiltonian(
        results, constrained_atoms, eigenvec_file, tolerance, verbose
    )


# =============================================================================
# Main: Test verification
# =============================================================================
if __name__ == "__main__":
    from pathlib import Path
    
    output_dir = Path(__file__).parent / "output"
    
    print("=" * 60)
    print("CDFTB Hamiltonian Verification")
    print("=" * 60)
    
    # Verify both fragments
    for frag_id in [1, 2]:
        print(f"\n--- Frame 0, Fragment {frag_id} ---")
        is_valid = verify_cdftb_from_dir(output_dir, 0, frag_id)
        print()
