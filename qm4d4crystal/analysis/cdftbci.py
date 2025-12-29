"""
CDFTB-CI Hamiltonian Construction Module (Unrestricted/Open-Shell Version)

This module provides functions to construct the 2x2 Hamiltonian matrix
for CDFTB-CI calculations with spin-polarized (unrestricted) wavefunctions.

For open-shell systems, α and β orbitals are treated separately:
- MO overlap: O^BA = O^BA_α ⊗ O^BA_β (determinant is product of α and β)
- State overlap: S_AB = det(O^BA_α) × det(O^BA_β)
- Weight overlap: W_BA = S_AB × [Tr(O^{-1,α}_BA Ω^α_BA) + Tr(O^{-1,β}_BA Ω^β_BA)]

References:
    - Wu and Van Voorhis, J. Chem. Phys. 125, 164105 (2006)
    - Wu et al., J. Chem. Phys. 127, 164119 (2007)
    - Oberhofer and Blumberger, J. Chem. Phys. 133, 244105 (2010)
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class UnrestrictedOrbitalData:
    """Container for unrestricted orbital data (α and β separately)."""
    C_alpha: np.ndarray  # α MO coefficients (n_AO x n_MO)
    C_beta: np.ndarray   # β MO coefficients (n_AO x n_MO)
    occ_alpha: np.ndarray  # α occupation numbers
    occ_beta: np.ndarray   # β occupation numbers
    n_alpha: int  # Number of α electrons
    n_beta: int   # Number of β electrons


@dataclass
class UnrestrictedCDFTBCIHamiltonian:
    """Container for CDFTB-CI Hamiltonian (unrestricted case)."""
    H: np.ndarray  # 2x2 Hamiltonian matrix
    S: np.ndarray  # 2x2 Overlap matrix (state overlap)
    E_A: float  # Energy of state A
    E_B: float  # Energy of state B
    H_AB: float  # Coupling element
    S_AB: float  # State overlap ⟨Φ^B|Φ^A⟩
    S_AB_alpha: float  # α contribution to state overlap
    S_AB_beta: float   # β contribution to state overlap
    W_BA: float  # Weight matrix overlap ⟨Φ^B|w^A|Φ^A⟩
    W_AB: float  # Weight matrix overlap ⟨Φ^A|w^B|Φ^B⟩
    W_BA_alpha: float  # α contribution to W_BA
    W_BA_beta: float   # β contribution to W_BA
    W_AB_alpha: float  # α contribution to W_AB
    W_AB_beta: float   # β contribution to W_AB
    V_A: float = 0.0  # Constraint potential for state A
    V_B: float = 0.0  # Constraint potential for state B
    N_A: float = 0.0  # Target population for constraint A
    N_B: float = 0.0  # Target population for constraint B
    
    def save(self, filepath: str | Path, frame_id: int = 0, time_fs: float = 0.0) -> None:
        """
        Save Hamiltonian parameters to file.
        
        Parameters
        ----------
        filepath : str or Path
            Output file path.
        frame_id : int
            Frame number.
        time_fs : float
            Time in femtoseconds.
        """
        filepath = Path(filepath)
        
        # Check if file exists to determine if we need header
        write_header = not filepath.exists()
        
        with open(filepath, 'a') as f:
            if write_header:
                f.write("# CDFTB-CI Hamiltonian Parameters\n")
                f.write("# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_A(Ha)  V_B(Ha)  N_A  N_B  "
                        "S_AB  S_AB_alpha  S_AB_beta  W_BA  W_AB  H_AB(Ha)  J(meV)\n")
            
            # Compute J
            J = self.H_AB - self.S_AB * (self.E_A + self.E_B) / 2
            J_meV = J * 27211.386
            
            f.write(f"{frame_id:6d}  {time_fs:8.2f}  {self.E_A:14.10f}  {self.E_B:14.10f}  "
                    f"{self.V_A:12.8f}  {self.V_B:12.8f}  {self.N_A:6.1f}  {self.N_B:6.1f}  "
                    f"{self.S_AB:12.8f}  {self.S_AB_alpha:12.8f}  {self.S_AB_beta:12.8f}  "
                    f"{self.W_BA:12.8f}  {self.W_AB:12.8f}  {self.H_AB:14.10f}  {J_meV:10.4f}\n")


def compute_mo_overlap_matrix_spin(
    C_A: np.ndarray,
    C_B: np.ndarray,
    S_AO: np.ndarray,
    n_occ: int
) -> np.ndarray:
    """
    Compute MO overlap matrix O^BA for one spin channel.
    
    O^BA_ij = ⟨φ_i^B|φ_j^A⟩ = Σ_μ Σ_ν C^B_μi S_μν C^A_νj
    
    Parameters
    ----------
    C_A : np.ndarray
        MO coefficients for state A (n_AO x n_MO).
    C_B : np.ndarray
        MO coefficients for state B (n_AO x n_MO).
    S_AO : np.ndarray
        AO overlap matrix (n_AO x n_AO).
    n_occ : int
        Number of occupied orbitals for this spin.
        
    Returns
    -------
    O_BA : np.ndarray
        MO overlap matrix (n_occ x n_occ).
    """
    if n_occ == 0:
        return np.array([[1.0]])  # Identity for empty subspace
    
    # Extract occupied orbitals
    C_A_occ = C_A[:, :n_occ]
    C_B_occ = C_B[:, :n_occ]
    
    # O^BA = C_B^T @ S_AO @ C_A
    O_BA = C_B_occ.T @ S_AO @ C_A_occ
    
    return O_BA


def compute_omega_matrix_spin(
    C_A: np.ndarray,
    C_B: np.ndarray,
    w: np.ndarray,
    n_occ: int
) -> np.ndarray:
    """
    Compute Ω^BA matrix for one spin channel.
    
    Ω^BA_ij = ⟨φ_i^B|w|φ_j^A⟩ = Σ_μ Σ_ν C^B_μi w_μν C^A_νj
    
    Parameters
    ----------
    C_A : np.ndarray
        MO coefficients for state A (n_AO x n_MO).
    C_B : np.ndarray
        MO coefficients for state B (n_AO x n_MO).
    w : np.ndarray
        Fragment weight matrix (n_AO x n_AO).
    n_occ : int
        Number of occupied orbitals for this spin.
        
    Returns
    -------
    Omega_BA : np.ndarray
        Ω^BA matrix (n_occ x n_occ).
    """
    if n_occ == 0:
        return np.array([[0.0]])
    
    # Extract occupied orbitals
    C_A_occ = C_A[:, :n_occ]
    C_B_occ = C_B[:, :n_occ]
    
    # Ω^BA = C_B^T @ w @ C_A
    Omega_BA = C_B_occ.T @ w @ C_A_occ
    
    return Omega_BA


def compute_state_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray
) -> Tuple[float, float, float]:
    """
    Compute state overlap for unrestricted wavefunctions.
    
    S_AB = det(O^BA_α) × det(O^BA_β)
    
    Parameters
    ----------
    O_BA_alpha : np.ndarray
        α MO overlap matrix.
    O_BA_beta : np.ndarray
        β MO overlap matrix.
        
    Returns
    -------
    S_AB : float
        Total state overlap.
    S_AB_alpha : float
        α contribution (determinant).
    S_AB_beta : float
        β contribution (determinant).
    """
    S_AB_alpha = np.linalg.det(O_BA_alpha)
    S_AB_beta = np.linalg.det(O_BA_beta)
    S_AB = S_AB_alpha * S_AB_beta
    
    return S_AB, S_AB_alpha, S_AB_beta


def compute_weight_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray,
    Omega_BA_alpha: np.ndarray,
    Omega_BA_beta: np.ndarray
) -> Tuple[float, float, float]:
    """
    Compute weight matrix overlap for unrestricted wavefunctions.
    
    For unrestricted:
    W_BA = ⟨Φ^B|w^A|Φ^A⟩ = S_AB × [Tr(O^{-1,α}_BA Ω^α_BA) + Tr(O^{-1,β}_BA Ω^β_BA)]
    
    Parameters
    ----------
    O_BA_alpha : np.ndarray
        α MO overlap matrix.
    O_BA_beta : np.ndarray
        β MO overlap matrix.
    Omega_BA_alpha : np.ndarray
        α Ω matrix.
    Omega_BA_beta : np.ndarray
        β Ω matrix.
        
    Returns
    -------
    W_BA : float
        Total weight overlap.
    W_BA_alpha : float
        α contribution to the trace.
    W_BA_beta : float
        β contribution to the trace.
    """
    # Compute determinants
    det_alpha = np.linalg.det(O_BA_alpha)
    det_beta = np.linalg.det(O_BA_beta)
    S_AB = det_alpha * det_beta
    
    # Compute trace terms for each spin
    # α contribution
    if O_BA_alpha.shape[0] > 0:
        O_inv_alpha = np.linalg.inv(O_BA_alpha)
        trace_alpha = np.trace(O_inv_alpha @ Omega_BA_alpha)
    else:
        trace_alpha = 0.0
    
    # β contribution
    if O_BA_beta.shape[0] > 0:
        O_inv_beta = np.linalg.inv(O_BA_beta)
        trace_beta = np.trace(O_inv_beta @ Omega_BA_beta)
    else:
        trace_beta = 0.0
    
    # Total weight overlap
    # W_BA = S_AB × (trace_α + trace_β)
    W_BA_alpha = S_AB * trace_alpha
    W_BA_beta = S_AB * trace_beta
    W_BA = W_BA_alpha + W_BA_beta
    
    return W_BA, W_BA_alpha, W_BA_beta


def compute_coupling_element_unrestricted(
    E_A: float,
    E_B: float,
    V_A: float,
    V_B: float,
    N_A: float,
    N_B: float,
    S_AB: float,
    W_BA: float,
    W_AB: float
) -> float:
    """
    Compute the coupling element H_AB (same formula as restricted).
    
    H_AB = 1/2 (E_A + E_B + N_A V_A + N_B V_B) S_AB - 1/2 (V_A W_BA + V_B W_AB)
    
    Parameters
    ----------
    E_A : float
        Total energy of state A (CDFTB energy).
    E_B : float
        Total energy of state B (CDFTB energy).
    V_A : float
        Constraint potential for state A.
    V_B : float
        Constraint potential for state B.
    N_A : float
        Target population for constraint A.
    N_B : float
        Target population for constraint B.
    S_AB : float
        State overlap ⟨Φ^B|Φ^A⟩.
    W_BA : float
        Weight overlap ⟨Φ^B|w^A|Φ^A⟩.
    W_AB : float
        Weight overlap ⟨Φ^A|w^B|Φ^B⟩.
        
    Returns
    -------
    H_AB : float
        Coupling element.
    """
    term1 = 0.5 * (E_A + E_B + N_A * V_A + N_B * V_B) * S_AB
    term2 = 0.5 * (V_A * W_BA + V_B * W_AB)
    
    H_AB = term1 - term2
    
    return H_AB


def build_cdftbci_hamiltonian_unrestricted(
    orb_A: UnrestrictedOrbitalData,
    orb_B: UnrestrictedOrbitalData,
    S_AO: np.ndarray,
    w_A: np.ndarray,
    w_B: np.ndarray,
    E_A: float,
    E_B: float,
    V_A: float,
    V_B: float,
    N_A: float,
    N_B: float
) -> UnrestrictedCDFTBCIHamiltonian:
    """
    Build the 2x2 CDFTB-CI Hamiltonian matrix for unrestricted wavefunctions.
    
    Parameters
    ----------
    orb_A : UnrestrictedOrbitalData
        Orbital data for state A (α and β MO coefficients and occupations).
    orb_B : UnrestrictedOrbitalData
        Orbital data for state B.
    S_AO : np.ndarray
        AO overlap matrix (n_AO x n_AO).
    w_A : np.ndarray
        Fragment weight matrix for state A.
    w_B : np.ndarray
        Fragment weight matrix for state B.
    E_A : float
        Total energy of state A.
    E_B : float
        Total energy of state B.
    V_A : float
        Constraint potential for state A.
    V_B : float
        Constraint potential for state B.
    N_A : float
        Target population for constraint A.
    N_B : float
        Target population for constraint B.
        
    Returns
    -------
    result : UnrestrictedCDFTBCIHamiltonian
        Container with Hamiltonian and related quantities.
    """
    # Compute MO overlap matrices for α and β
    O_BA_alpha = compute_mo_overlap_matrix_spin(
        orb_A.C_alpha, orb_B.C_alpha, S_AO, orb_A.n_alpha
    )
    O_BA_beta = compute_mo_overlap_matrix_spin(
        orb_A.C_beta, orb_B.C_beta, S_AO, orb_A.n_beta
    )
    
    # Transpose for O_AB
    O_AB_alpha = O_BA_alpha.T.conj()
    O_AB_beta = O_BA_beta.T.conj()
    
    # Compute state overlap
    S_AB, S_AB_alpha, S_AB_beta = compute_state_overlap_unrestricted(
        O_BA_alpha, O_BA_beta
    )
    
    # Compute Ω matrices for W_BA (using w_A)
    Omega_BA_alpha = compute_omega_matrix_spin(
        orb_A.C_alpha, orb_B.C_alpha, w_A, orb_A.n_alpha
    )
    Omega_BA_beta = compute_omega_matrix_spin(
        orb_A.C_beta, orb_B.C_beta, w_A, orb_A.n_beta
    )
    
    # Compute Ω matrices for W_AB (using w_B)
    Omega_AB_alpha = compute_omega_matrix_spin(
        orb_B.C_alpha, orb_A.C_alpha, w_B, orb_B.n_alpha
    )
    Omega_AB_beta = compute_omega_matrix_spin(
        orb_B.C_beta, orb_A.C_beta, w_B, orb_B.n_beta
    )
    
    # Compute weight overlaps
    W_BA, W_BA_alpha, W_BA_beta = compute_weight_overlap_unrestricted(
        O_BA_alpha, O_BA_beta, Omega_BA_alpha, Omega_BA_beta
    )
    W_AB, W_AB_alpha, W_AB_beta = compute_weight_overlap_unrestricted(
        O_AB_alpha, O_AB_beta, Omega_AB_alpha, Omega_AB_beta
    )
    
    # Compute coupling element
    H_AB = compute_coupling_element_unrestricted(
        E_A, E_B, V_A, V_B, N_A, N_B, S_AB, W_BA, W_AB
    )
    
    # Build 2x2 matrices
    H = np.array([
        [E_A, H_AB],
        [H_AB, E_B]
    ])
    
    S = np.array([
        [1.0, S_AB],
        [S_AB, 1.0]
    ])
    
    return UnrestrictedCDFTBCIHamiltonian(
        H=H,
        S=S,
        E_A=E_A,
        E_B=E_B,
        H_AB=H_AB,
        S_AB=S_AB,
        S_AB_alpha=S_AB_alpha,
        S_AB_beta=S_AB_beta,
        W_BA=W_BA,
        W_AB=W_AB,
        W_BA_alpha=W_BA_alpha,
        W_BA_beta=W_BA_beta,
        W_AB_alpha=W_AB_alpha,
        W_AB_beta=W_AB_beta,
        V_A=V_A,
        V_B=V_B,
        N_A=N_A,
        N_B=N_B
    )


def solve_cdftbci_unrestricted(H: np.ndarray, S: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve the generalized eigenvalue problem for CDFTB-CI.
    
    HC = SCE
    
    Parameters
    ----------
    H : np.ndarray
        2x2 Hamiltonian matrix.
    S : np.ndarray
        2x2 Overlap matrix.
        
    Returns
    -------
    eigenvalues : np.ndarray
        Eigenvalues (adiabatic state energies).
    eigenvectors : np.ndarray
        Eigenvectors (mixing coefficients).
    """
    from scipy.linalg import eigh
    
    eigenvalues, eigenvectors = eigh(H, S)
    
    return eigenvalues, eigenvectors


def compute_transfer_integral_unrestricted(
    H: np.ndarray,
    S: np.ndarray,
    method: str = "direct"
) -> float:
    """
    Compute the effective transfer integral.
    
    Parameters
    ----------
    H : np.ndarray
        2x2 Hamiltonian matrix.
    S : np.ndarray
        2x2 Overlap matrix.
    method : str
        Method: "direct" or "lowdin".
        
    Returns
    -------
    J : float
        Transfer integral in the same units as H.
    """
    E_A = H[0, 0]
    E_B = H[1, 1]
    H_AB = H[0, 1]
    S_AB = S[0, 1]
    
    J_direct = H_AB - S_AB * (E_A + E_B) / 2
    
    if method == "direct":
        return J_direct
    elif method == "lowdin":
        J_lowdin = J_direct / (1 - S_AB**2)
        return J_lowdin
    else:
        raise ValueError(f"Unknown method: {method}")


# =============================================================================
# I/O Functions for Unrestricted Calculations
# =============================================================================

def read_eigenvalues_spin_polarized(filepath: str | Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read eigenvalues and occupations from DFTB+ band.out file (spin-polarized).
    
    Parameters
    ----------
    filepath : str or Path
        Path to band.out file.
        
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
        
    Notes
    -----
    File format for spin-polarized:
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
                # Parse spin from header
                parts = line.split()
                spin_idx = parts.index('SPIN') + 1
                current_spin = int(parts[spin_idx])
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
        # Return alpha values for both (restricted case)
        return (np.array(eigenvalues_alpha), np.array(occupations_alpha),
                np.array(eigenvalues_alpha), np.array(occupations_alpha))
    
    return (np.array(eigenvalues_alpha), np.array(occupations_alpha),
            np.array(eigenvalues_beta), np.array(occupations_beta))


def read_eigenvector_matrix_spin_polarized(filepath: str | Path) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Read eigenvector matrices from spin-polarized DFTB+ calculation.
    
    For spin-polarized calculations, DFTB+ may output separate eigenvector
    matrices for α and β spins, or a combined file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to eigenvec_matrix.dat file.
        
    Returns
    -------
    C_alpha : np.ndarray
        α eigenvector matrix (n_AO x n_MO).
    C_beta : np.ndarray
        β eigenvector matrix (n_AO x n_MO).
    n_orbitals : int
        Number of AOs.
        
    Notes
    -----
    The format depends on DFTB+ version and settings.
    If separate files exist (eigenvec_matrix_spin1.dat, eigenvec_matrix_spin2.dat),
    those will be read. Otherwise, assumes the standard file contains both.
    """
    filepath = Path(filepath)
    dirpath = filepath.parent
    
    # Check for separate spin files
    spin1_file = dirpath / "eigenvec_matrix_spin1.dat"
    spin2_file = dirpath / "eigenvec_matrix_spin2.dat"
    
    if spin1_file.exists() and spin2_file.exists():
        # Read separate files
        C_alpha, n_orbitals = _read_eigenvector_matrix_single(spin1_file)
        C_beta, _ = _read_eigenvector_matrix_single(spin2_file)
        return C_alpha, C_beta, n_orbitals
    
    # Try to read from combined file
    # DFTB+ may output both spins in one file with doubled columns
    if not filepath.exists():
        raise FileNotFoundError(f"Eigenvector matrix file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header_line = lines[1].strip()
    parts = header_line.split()
    
    is_real = parts[0].upper() == 'T'
    n_orbitals = int(parts[1])
    n_eigvec = int(parts[2])
    
    if not is_real:
        raise ValueError("Complex eigenvectors not supported")
    
    # Find eigenvector section
    eigenvec_start = None
    for i, line in enumerate(lines):
        if '# EIGENVECTORS' in line or '# SPIN' in line:
            eigenvec_start = i + 1
            break
    
    if eigenvec_start is None:
        eigenvec_start = 2
    
    # Check if this is a combined spin file (n_eigvec = 2 * n_orbitals)
    if n_eigvec == 2 * n_orbitals:
        # Combined file: first n_orbitals columns are α, next n_orbitals are β
        eigenvectors = np.zeros((n_orbitals, n_eigvec), dtype=np.float64)
        
        for i in range(n_orbitals):
            line_idx = eigenvec_start + i
            row_data = lines[line_idx].split()
            eigenvectors[i, :] = [float(x) for x in row_data]
        
        C_alpha = eigenvectors[:, :n_orbitals]
        C_beta = eigenvectors[:, n_orbitals:]
        
        return C_alpha, C_beta, n_orbitals
    else:
        # Single spin file or restricted calculation
        C, n_orbitals = _read_eigenvector_matrix_single(filepath)
        return C, C.copy(), n_orbitals


def _read_eigenvector_matrix_single(filepath: str | Path) -> Tuple[np.ndarray, int]:
    """Read a single eigenvector matrix file."""
    filepath = Path(filepath)
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    header_line = lines[1].strip()
    parts = header_line.split()
    
    n_orbitals = int(parts[1])
    n_eigvec = int(parts[2])
    
    eigenvec_start = None
    for i, line in enumerate(lines):
        if '# EIGENVECTORS' in line:
            eigenvec_start = i + 1
            break
    
    if eigenvec_start is None:
        eigenvec_start = 2
    
    eigenvectors = np.zeros((n_orbitals, n_eigvec), dtype=np.float64)
    
    for i in range(n_orbitals):
        line_idx = eigenvec_start + i
        row_data = lines[line_idx].split()
        eigenvectors[i, :] = [float(x) for x in row_data]
    
    return eigenvectors, n_orbitals


def load_unrestricted_orbital_data(
    dirpath: str | Path,
    eigenvec_filename: str = "eigenvec_matrix.dat",
    band_filename: str = "band.out"
) -> UnrestrictedOrbitalData:
    """
    Load unrestricted orbital data from DFTB+ output directory.
    
    Parameters
    ----------
    dirpath : str or Path
        Path to DFTB+ output directory.
    eigenvec_filename : str
        Name of eigenvector matrix file.
    band_filename : str
        Name of band (eigenvalue) file.
        
    Returns
    -------
    data : UnrestrictedOrbitalData
        Container with α and β orbital data.
    """
    dirpath = Path(dirpath)
    
    # Read eigenvectors
    eigenvec_file = dirpath / eigenvec_filename
    C_alpha, C_beta, n_orbitals = read_eigenvector_matrix_spin_polarized(eigenvec_file)
    
    # Read eigenvalues and occupations
    band_file = dirpath / band_filename
    eig_alpha, occ_alpha, eig_beta, occ_beta = read_eigenvalues_spin_polarized(band_file)
    
    # Count occupied orbitals
    n_alpha = np.sum(occ_alpha > 0.5).astype(int)
    n_beta = np.sum(occ_beta > 0.5).astype(int)
    
    return UnrestrictedOrbitalData(
        C_alpha=C_alpha,
        C_beta=C_beta,
        occ_alpha=occ_alpha,
        occ_beta=occ_beta,
        n_alpha=n_alpha,
        n_beta=n_beta
    )


# =============================================================================
# CDFTB-CI Analysis with Online Calculation
# =============================================================================

import dftbplus
import hsd
import mdtraj as md
import shutil
import yaml
from dataclasses import field
from typing import List, Dict, Any

from .cdftb import (
    CDFTBConfig,
    FragmentConfig,
    load_config,
    load_qm_indices,
    load_mm_indices,
    load_pccharges,
    iter_qm_coordinates,
    extract_atom_types_from_hsd,
    save_qm_coords_xyz,
    run_dftb_in_subprocess,
    parse_atom_range,
    ANG_PER_NM,
    BOHR_PER_ANG,
)

from .cdftb_result_reader import (
    load_spin_polarized_calculation_for_ci,
    get_orbital_info_from_eigenvec,
    get_atom_orbital_map,
    build_fragment_weight_matrix,
)

from .spin import compute_fragment_spin_population


@dataclass
class CDFTBCIConfig(CDFTBConfig):
    """Configuration for CDFTB-CI analysis."""
    # CI settings
    ci_enabled: bool = True
    ci_N_A: float = 101.0  # Target population for constraint A
    ci_N_B: float = 101.0  # Target population for constraint B
    
    # Output mode
    output_mode: str = "online_ci"  # "online_ci" or "store_frames"
    work_directory: str = "work"
    
    # CI output files
    ci_output_file: Optional[Path] = None
    ci_sub_file: Optional[Path] = None
    
    # Spin output file
    spin_output_file: Optional[Path] = None
    
    # Retry log file
    retry_log_file: Optional[Path] = None
    
    # Use previous frame's charges as initial guess
    use_previous_charges: bool = False


def load_cdftbci_config(config_path: Path) -> CDFTBCIConfig:
    """Load configuration from YAML file with CI settings."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse input paths (relative to config file) - all required
    input_cfg = data.get('input', {})
    if 'trajectory' not in input_cfg:
        raise ValueError("Missing required input: trajectory")
    if 'topology' not in input_cfg:
        raise ValueError("Missing required input: topology")
    if 'qm_atoms' not in input_cfg:
        raise ValueError("Missing required input: qm_atoms")
    if 'pccharges_template' not in input_cfg:
        raise ValueError("Missing required input: pccharges_template")
    if 'hsd_template' not in input_cfg:
        raise ValueError("Missing required input: hsd_template")
    
    traj_path = base_dir / input_cfg['trajectory']
    topology_path = base_dir / input_cfg['topology']
    qm_atoms_file = base_dir / input_cfg['qm_atoms']
    pccharges_template = base_dir / input_cfg['pccharges_template']
    hsd_template = base_dir / input_cfg['hsd_template']
    
    # Parse output paths
    output_cfg = data.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output')
    energy_filename = output_cfg.get('energy_file', 'energies.dat')
    energy_file = output_dir / energy_filename
    charge_filename = output_cfg.get('charge_file', 'charges.dat')
    charge_file = output_dir / charge_filename
    
    # Output mode settings
    output_mode = output_cfg.get('mode', 'online_ci')
    work_directory = output_cfg.get('work_directory', 'work')
    
    # Parse frame settings
    frames_cfg = data.get('frames', {})
    start_frame = frames_cfg.get('start', 0)
    n_frames = frames_cfg.get('n_frames', None)
    t0_fs = frames_cfg.get('t0_fs', 0.0)
    dt_fs = frames_cfg.get('dt_fs', 4.0)
    
    # Parse DFTB+ settings
    dftb_cfg = data.get('dftb', {})
    dftb_library_path = dftb_cfg.get('library_path', '/home/takahashi/opt/dftb+/lib/libdftbplus.so')
    num_threads = dftb_cfg.get('num_threads', None)
    timeout = dftb_cfg.get('timeout', 300)
    
    # Parse fragment definitions
    fragments = []
    for frag_data in data.get('fragments', []):
        atom_range = frag_data['atom_range']
        if 'charge_sum_range' in frag_data:
            charge_sum_range = frag_data['charge_sum_range']
        else:
            charge_sum_range = parse_atom_range(atom_range)
        # Parse initial_charges path if provided
        initial_charges = None
        if 'initial_charges' in frag_data:
            initial_charges = base_dir / frag_data['initial_charges']
        frag = FragmentConfig(
            name=frag_data['name'],
            atom_range=atom_range,
            charge_sum_range=charge_sum_range,
            initial_charges=initial_charges
        )
        fragments.append(frag)
    
    # Parse CI settings
    ci_cfg = data.get('cdftb_ci', {})
    ci_enabled = ci_cfg.get('enabled', True)
    ci_N_A = ci_cfg.get('N_A', 101.0)
    ci_N_B = ci_cfg.get('N_B', 101.0)
    
    # Parse SCC settings
    scc_cfg = data.get('scc', {})
    use_previous_charges = scc_cfg.get('use_previous_charges', False)
    
    # CI output files
    ci_output_file = output_dir / output_cfg.get('ci_file', 'cdftbci.dat')
    ci_sub_file = output_dir / output_cfg.get('ci_sub_file', 'cdftbci_sub.dat')
    
    # Spin output file
    spin_output_file = output_dir / output_cfg.get('spin_file', 'spin.dat')
    
    # Retry log file
    retry_log_file = output_dir / output_cfg.get('retry_log_file', 'retry_log.dat')
    
    return CDFTBCIConfig(
        traj_path=traj_path,
        topology_path=topology_path,
        qm_atoms_file=qm_atoms_file,
        pccharges_template=pccharges_template,
        hsd_template=hsd_template,
        output_dir=output_dir,
        energy_file=energy_file,
        charge_file=charge_file,
        start_frame=start_frame,
        n_frames=n_frames,
        t0_fs=t0_fs,
        dt_fs=dt_fs,
        dftb_library_path=dftb_library_path,
        num_threads=num_threads,
        timeout=timeout,
        fragments=fragments,
        ci_enabled=ci_enabled,
        ci_N_A=ci_N_A,
        ci_N_B=ci_N_B,
        output_mode=output_mode,
        work_directory=work_directory,
        ci_output_file=ci_output_file,
        ci_sub_file=ci_sub_file,
        spin_output_file=spin_output_file,
        retry_log_file=retry_log_file,
        use_previous_charges=use_previous_charges,
    )


def setup_work_directory(
    work_dir: Path,
    qm_coords_ang: np.ndarray,
    mm_coords_ang: np.ndarray,
    mm_charges: np.ndarray,
    frame_id: int,
    time_fs: float,
    atom_types: List[str],
) -> None:
    """
    Set up or update the working directory with current frame data.
    
    This function reuses the same directory for each frame, only updating
    the coordinate files.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Save QM coordinates in XYZ format
    save_qm_coords_xyz(work_dir, frame_id, time_fs, qm_coords_ang, atom_types)
    
    # Save PCcharges.dat with updated MM coordinates
    pc_file = work_dir / "PCcharges.dat"
    pc_data = np.column_stack([mm_coords_ang, mm_charges])
    np.savetxt(pc_file, pc_data, fmt="%20.10f %20.10f %20.10f %10.4f")


def setup_fragment_work_directory(
    work_dir: Path,
    fragment_name: str,
    constrained_atoms: str,
    hsd_template: Path,
    read_initial_charges: bool = False,
    disable_constraint: bool = False,
) -> Path:
    """
    Set up fragment subdirectory within work directory.
    
    Parameters
    ----------
    work_dir : Path
        Parent working directory.
    fragment_name : str
        Name of the fragment (used for subdirectory name).
    constrained_atoms : str
        DFTB+ style atom range for constraint, e.g., "1:36".
    hsd_template : Path
        Path to HSD template file.
    read_initial_charges : bool
        Whether to read initial charges from charges.dat.
    disable_constraint : bool
        If True, disable the electronic constraint (for fallback calculations).
    """
    frag_dir = work_dir / fragment_name
    frag_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy and modify HSD template
    with open(hsd_template, 'r') as f:
        data = hsd.load(f)
    
    # Remove Geometry from data (will be added as text later)
    if 'Geometry' in data:
        del data['Geometry']
    
    # Update PCcharges path (use parent directory's file)
    data['Hamiltonian']['DFTB']['ElectricField']['PointCharges']['CoordsAndCharges']['DirectRead']['File'] = '../PCcharges.dat'
    
    # Disable constraint if requested (for fallback calculations)
    if disable_constraint:
        if 'ElectronicConstraints' in data['Hamiltonian']['DFTB']:
            del data['Hamiltonian']['DFTB']['ElectronicConstraints']
    else:
        # Modify ElectronicConstraints to constrain specific fragment
        if 'ElectronicConstraints' in data['Hamiltonian']['DFTB']:
            constraints = data['Hamiltonian']['DFTB']['ElectronicConstraints']['Constraints']
            if 'MullikenPopulation' in constraints:
                mulliken = constraints['MullikenPopulation']
                if isinstance(mulliken, list):
                    for constraint in mulliken:
                        constraint['Atoms'] = constrained_atoms
                else:
                    mulliken['Atoms'] = constrained_atoms
    
    # Modify InitialSpins.Atoms to guide spin localization
    if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
        spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
        if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
            atom_spin = spin_pol['Colinear']['InitialSpins']['AtomSpin']
            atom_spin['Atoms'] = constrained_atoms
    
    # Set ReadInitialCharges if using previous frame's charges
    if read_initial_charges:
        data['Hamiltonian']['DFTB']['ReadInitialCharges'] = 'Yes'
        # Need ReadChargesAsText to read text format charges.dat
        if 'Options' not in data:
            data['Options'] = {}
        data['Options']['ReadChargesAsText'] = 'Yes'
        # Remove InitialSpins when using ReadInitialCharges (they conflict)
        if 'SpinPolarisation' in data['Hamiltonian']['DFTB']:
            spin_pol = data['Hamiltonian']['DFTB']['SpinPolarisation']
            if 'Colinear' in spin_pol and 'InitialSpins' in spin_pol['Colinear']:
                del spin_pol['Colinear']['InitialSpins']
    
    # Save modified HSD
    hsd_file = frag_dir / "dftb_in.hsd"
    with open(hsd_file, 'w') as f:
        # Write Geometry section with file include directive first
        f.write('Geometry {\n')
        f.write('  xyzFormat {\n')
        f.write('    <<< "../qm_coords.xyz"\n')
        f.write('  }\n')
        f.write('}\n\n')
        # Then write the rest of the HSD
        hsd.dump(data, f)
    
    return frag_dir


def compute_cdftbci_for_frame(
    work_dir: Path,
    frag1_name: str,
    frag2_name: str,
    E_A: float,
    E_B: float,
    N_A: float,
    N_B: float,
) -> Tuple[Optional[UnrestrictedCDFTBCIHamiltonian], Optional[float], Optional[float], Optional[str]]:
    """
    Compute CDFTB-CI quantities for current frame using data in work directory.
    """
    frag1_dir = work_dir / frag1_name
    frag2_dir = work_dir / frag2_name
    
    try:
        # Load spin-polarized calculation results (without Hamiltonian matrices)
        data_A = load_spin_polarized_calculation_for_ci(frag1_dir)
        data_B = load_spin_polarized_calculation_for_ci(frag2_dir)
        
        # Create UnrestrictedOrbitalData objects
        orb_A = UnrestrictedOrbitalData(
            C_alpha=data_A['C_alpha'],
            C_beta=data_A['C_beta'],
            occ_alpha=data_A['occ_alpha'],
            occ_beta=data_A['occ_beta'],
            n_alpha=data_A['n_alpha'],
            n_beta=data_A['n_beta']
        )
        orb_B = UnrestrictedOrbitalData(
            C_alpha=data_B['C_alpha'],
            C_beta=data_B['C_beta'],
            occ_alpha=data_B['occ_alpha'],
            occ_beta=data_B['occ_beta'],
            n_alpha=data_B['n_alpha'],
            n_beta=data_B['n_beta']
        )
        
        # Read AO overlap matrix
        S_AO = data_A['S']
        n_orbitals = data_A['n_orbitals']
        
        # Build fragment weight matrices
        eigenvec_file = frag1_dir / "eigenvec.out"
        orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
        atom_to_orbitals = get_atom_orbital_map(orbital_info)
        n_atoms = len(atom_to_orbitals)
        n_atoms_half = n_atoms // 2
        
        frag_A_atoms = list(range(1, n_atoms_half + 1))
        frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
        
        w_A = build_fragment_weight_matrix(
            S_AO, frag_A_atoms, n_atoms, atom_to_orbitals=atom_to_orbitals
        )
        w_B = build_fragment_weight_matrix(
            S_AO, frag_B_atoms, n_atoms, atom_to_orbitals=atom_to_orbitals
        )
        
        # Get constraint potentials
        V_A = data_A['Vc']
        V_B = data_B['Vc']
        
        # Build CDFTB-CI Hamiltonian
        ham = build_cdftbci_hamiltonian_unrestricted(
            orb_A, orb_B, S_AO, w_A, w_B,
            E_A, E_B, V_A, V_B, N_A, N_B
        )
        
        # Compute transfer integrals
        J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
        J_lowdin = compute_transfer_integral_unrestricted(ham.H, ham.S, method="lowdin")
        
        return ham, J_direct, J_lowdin, None
        
    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        return None, None, None, error_msg


def compute_spin_for_fragment(
    frag_dir: Path,
    n_atoms_half: int
) -> Tuple[Optional[Dict[str, Dict[str, float]]], Optional[str]]:
    """
    Compute spin expectation values for both fragments from a CDFTB calculation.
    """
    try:
        # Load spin-polarized data
        data = load_spin_polarized_calculation_for_ci(frag_dir)
        
        # Get orbital-to-atom mapping
        eigenvec_file = frag_dir / "eigenvec.out"
        orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
        atom_to_orbitals = get_atom_orbital_map(orbital_info)
        n_atoms = len(atom_to_orbitals)
        
        # Define fragment atoms (1-indexed)
        frag_A_atoms = list(range(1, n_atoms_half + 1))
        frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
        
        # Compute spin populations
        spin_frag1 = compute_fragment_spin_population(
            data['C_alpha'], data['C_beta'], data['S'],
            data['n_alpha'], data['n_beta'],
            frag_A_atoms, atom_to_orbitals
        )
        spin_frag2 = compute_fragment_spin_population(
            data['C_alpha'], data['C_beta'], data['S'],
            data['n_alpha'], data['n_beta'],
            frag_B_atoms, atom_to_orbitals
        )
        
        return {'frag1': spin_frag1, 'frag2': spin_frag2}, None
        
    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        return None, error_msg


def run_cdftb_with_retry(
    work_dir: Path,
    frag: 'FragmentConfig',
    qm_coords_bohr: np.ndarray,
    hsd_template: Path,
    dftb_library_path: str,
    num_threads: Optional[int],
    timeout: int,
    use_initial_charges: bool,
    ci_enabled: bool,
) -> Tuple[float, List[float], bool, str]:
    """
    Run CDFTB calculation with automatic retry on convergence failure.
    
    Retry strategy:
    1. First attempt: Run with specified initial charges setting
    2. If failed: Run DFTB to generate initial charges, then retry CDFTB with those charges
    3. If still failed: Run without constraint and without reading initial charges
    
    Parameters
    ----------
    work_dir : Path
        Working directory.
    frag : FragmentConfig
        Fragment configuration.
    qm_coords_bohr : np.ndarray
        QM coordinates in Bohr.
    hsd_template : Path
        Path to HSD template file.
    dftb_library_path : str
        Path to DFTB+ library.
    num_threads : int or None
        Number of OpenMP threads.
    timeout : int
        Timeout in seconds.
    use_initial_charges : bool
        Whether to read initial charges for first attempt.
    ci_enabled : bool
        Whether CDFTB-CI is enabled (determines if WriteHS is needed).
        
    Returns
    -------
    energy : float
        Energy (nan if all attempts failed).
    all_frag_charges : list of float
        Mulliken charges for all fragments (nan if failed).
    success : bool
        Whether calculation succeeded.
    retry_info : str
        Information about retry attempts ("", "no_init_charges", "all_failed:...").
    """
    frag_dir = work_dir / frag.name
    retry_info = ""
    
    # First attempt: normal CDFTB
    frag_dir = setup_fragment_work_directory(
        work_dir, frag.name, frag.atom_range, hsd_template,
        read_initial_charges=use_initial_charges
    )
    
    energy, mcharge, error = run_dftb_in_subprocess(
        qm_coords_bohr, frag_dir, write_hs=False,
        dftb_library_path=dftb_library_path,
        num_threads=num_threads,
        timeout=timeout
    )
    
    if error is None:
        # Success on first attempt
        return energy, mcharge, True, retry_info
    
    # Second attempt: Run with constraint but without reading initial charges
    print(f"    [Retry] First attempt failed ({error[:50]}...), trying without initial charges")
    
    frag_dir = setup_fragment_work_directory(
        work_dir, frag.name, frag.atom_range, hsd_template,
        read_initial_charges=False,
        disable_constraint=False
    )
    
    energy, mcharge, error = run_dftb_in_subprocess(
        qm_coords_bohr, frag_dir, write_hs=False,
        dftb_library_path=dftb_library_path,
        num_threads=num_threads,
        timeout=timeout
    )
    
    if error is None:
        retry_info = "no_init_charges"
        return energy, mcharge, True, retry_info
    
    # All attempts failed
    return float("nan"), None, False, f"all_failed:{error}"


def run_cdftbci_analysis(config_path: Path) -> None:
    """
    Run CDFTB-CI analysis with online calculation.
    
    This is the main entry point for CDFTB-CI workflow.
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file.
    """
    # Load configuration
    config = load_cdftbci_config(config_path)
    
    print("=" * 70)
    print("CDFTB-CI Analysis")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Trajectory: {config.traj_path}")
    print(f"Topology: {config.topology_path}")
    print(f"Output directory: {config.output_dir}")
    print(f"CI enabled: {config.ci_enabled}")
    print("=" * 70)
    
    # Load QM/MM indices
    qm_indices = load_qm_indices(config.qm_atoms_file)
    print(f"Loaded {len(qm_indices)} QM atom indices from {config.qm_atoms_file}")
    
    traj_info = md.load(str(config.traj_path), top=str(config.topology_path), frame=0)
    total_atoms = traj_info.n_atoms
    mm_indices = load_mm_indices(qm_indices, total_atoms)
    print(f"MM atoms: {len(mm_indices)}")
    
    # Load point charges
    mm_charges = load_pccharges(config.pccharges_template)
    print(f"Loaded {len(mm_charges)} point charges from {config.pccharges_template}")
    
    # Extract atom types from HSD template
    atom_types = extract_atom_types_from_hsd(config.hsd_template)
    print(f"Extracted {len(atom_types)} atom types from HSD template")
    
    # Print fragment information
    print(f"\nFragments ({len(config.fragments)}):")
    for frag in config.fragments:
        print(f"  - {frag.name}: atoms {frag.atom_range}, charge sum range {frag.charge_sum_range}")
    
    if config.ci_enabled:
        print(f"\nCDFTB-CI settings:")
        print(f"  N_A = {config.ci_N_A}")
        print(f"  N_B = {config.ci_N_B}")
    
    if config.use_previous_charges:
        print(f"\nSCC settings:")
        print(f"  Use previous charges: Yes")
    print()
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up work directory path
    work_dir = config.output_dir / config.work_directory
    
    # Build headers for output files
    energy_header_parts = ["# Frame", "Time[fs]"]
    for frag in config.fragments:
        energy_header_parts.append(f"Energy_{frag.name}[a.u.]")
    energy_header = "  ".join(energy_header_parts) + "\n"
    
    charge_header_lines = []
    charge_header_lines.append("# Mulliken charges from constrained DFT calculations")
    charge_header_lines.append("# Qn|calcm = charge on fragment n when constraint applied to fragment m")
    charge_header_parts = [f"{'# Frame':>7s}", f"{'Time[fs]':>12s}"]
    for i, frag in enumerate(config.fragments, 1):
        for j, frag2 in enumerate(config.fragments, 1):
            charge_header_parts.append(f"{'Q' + str(j) + '|calc' + str(i):>12s}")
    charge_header_lines.append("  ".join(charge_header_parts))
    charge_header = "\n".join(charge_header_lines) + "\n"
    
    # Open output files
    energy_file = open(config.energy_file, "w")
    energy_file.write(energy_header)
    
    charge_file = open(config.charge_file, "w")
    charge_file.write(charge_header)
    
    ci_file = None
    ci_sub_file = None
    spin_file = None
    retry_log_file = None
    
    # Open retry log file
    retry_log_file = open(config.retry_log_file, "w")
    retry_log_file.write("# Retry Log for CDFTB Calculations\n")
    retry_log_file.write("# Retry info: '-' = success on first attempt,\n")
    retry_log_file.write("#             'no_init_charges' = succeeded without initial charges, 'all_failed:...' = all attempts failed\n")
    frag_names = [f.name for f in config.fragments]
    retry_header = "# Frame  Time(fs)  " + "  ".join([f"retry_{name}" for name in frag_names]) + "\n"
    retry_log_file.write(retry_header)
    
    if config.ci_enabled:
        ci_file = open(config.ci_output_file, "w")
        ci_file.write("# CDFTB-CI Results\n")
        ci_file.write("# Frame  Time(fs)   J_lowdin(meV)        E1(Ha)        E2(Ha)        dE(eV)\n")
        
        ci_sub_file = open(config.ci_sub_file, "w")
        ci_sub_file.write("# CDFTB-CI Sub Values\n")
        ci_sub_file.write("# Frame  Time(fs)      E_A(Ha)      E_B(Ha)        H_AB(Ha)        J_direct(meV)        "
                         "V_A(Ha)        V_B(Ha)      N_A      N_B        S_AB        "
                         "S_AB_alpha      S_AB_beta        W_BA        W_BA_alpha      W_BA_beta        "
                         "W_AB        W_AB_alpha      W_AB_beta\n")
    
    # Always open spin output file (useful even without CI)
    spin_file = open(config.spin_output_file, "w")
    spin_file.write("# Spin Populations from CDFTB Calculations\n")
    spin_file.write("# Spin density = N_α - N_β on each fragment\n")
    spin_file.write("# For each constraint state (A or B):\n")
    spin_file.write("#   frag1 = first fragment atoms, frag2 = second fragment atoms\n")
    frag1_name = config.fragments[0].name if len(config.fragments) > 0 else "frag1"
    frag2_name = config.fragments[1].name if len(config.fragments) > 1 else "frag2"
    spin_file.write(f"# Frame  Time(fs)  "
                   f"{frag1_name}_A_nalpha  {frag1_name}_A_nbeta  {frag1_name}_A_spin  "
                   f"{frag2_name}_A_nalpha  {frag2_name}_A_nbeta  {frag2_name}_A_spin  "
                   f"{frag1_name}_B_nalpha  {frag1_name}_B_nbeta  {frag1_name}_B_spin  "
                   f"{frag2_name}_B_nalpha  {frag2_name}_B_nbeta  {frag2_name}_B_spin\n")
    
    try:
        processed_count = 0
        
        for frame_id, time_fs, qm_coords_bohr, mm_coords_ang in iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices, mm_indices
        ):
            # Skip frames before start_frame
            if frame_id < config.start_frame:
                continue
            
            # Stop if we've processed enough frames
            if config.n_frames is not None and processed_count >= config.n_frames:
                break
            
            # Convert QM coords to Angstrom
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Calculate time using user-defined t0 and dt
            time_val = config.t0_fs + frame_id * config.dt_fs
            time_str = f"t = {time_val:.3f} fs"
            
            print(f"Frame {frame_id:05d} ({time_str})")
            
            # Set up work directory with current coordinates
            setup_work_directory(
                work_dir, qm_coords_ang, mm_coords_ang, mm_charges,
                frame_id, time_fs, atom_types
            )
            
            # Process each fragment
            energies = []
            charges = []
            cdftb_success = True
            retry_infos = []
            
            for frag in config.fragments:
                frag_dir = work_dir / frag.name
                charges_file = frag_dir / "charges.dat"
                
                # Determine if we should use initial charges
                use_initial_charges = False
                
                if config.use_previous_charges:
                    if processed_count == 0:
                        # First frame: use initial_charges from config if provided
                        if frag.initial_charges is not None and frag.initial_charges.exists():
                            use_initial_charges = True
                            # Copy initial charges to work directory
                            frag_dir.mkdir(parents=True, exist_ok=True)
                            shutil.copy(frag.initial_charges, charges_file)
                    else:
                        # Subsequent frames: use previous frame's charges if they exist
                        if charges_file.exists():
                            use_initial_charges = True
                
                # Run CDFTB calculation with automatic retry
                energy, mcharge, success, retry_info = run_cdftb_with_retry(
                    work_dir, frag, qm_coords_bohr, config.hsd_template,
                    config.dftb_library_path, config.num_threads, config.timeout,
                    use_initial_charges, config.ci_enabled
                )
                retry_infos.append(retry_info)
                
                if not success:
                    print(f"  {frag.name}: ERROR - {retry_info}")
                    frag_charge = float("nan")
                    energy = float("nan")
                    all_frag_charges = [float("nan")] * len(config.fragments)
                    cdftb_success = False
                else:
                    start_idx, end_idx = frag.charge_sum_range
                    frag_charge = float(np.sum(mcharge[start_idx:end_idx]))
                    all_frag_charges = []
                    for f in config.fragments:
                        s, e = f.charge_sum_range
                        all_frag_charges.append(float(np.sum(mcharge[s:e])))
                    charges_str = ", ".join([f"Q_{f.name}={q:+.6f}" for f, q in zip(config.fragments, all_frag_charges)])
                    retry_str = f" [{retry_info}]" if retry_info else ""
                    print(f"  {frag.name}: E = {energy:.10f} a.u., {charges_str}{retry_str}")
                    
                    # Run WriteHS calculation to output oversqr.dat for CDFTB-CI
                    # Note: DFTB+ crashes after WriteHS output, which is expected behavior
                    if config.ci_enabled:
                        frag_dir = work_dir / frag.name
                        _, _, error_hs = run_dftb_in_subprocess(
                            qm_coords_bohr, frag_dir, write_hs=True,
                            dftb_library_path=config.dftb_library_path,
                            num_threads=config.num_threads,
                            timeout=config.timeout
                        )
                        # Ignore WriteHS errors as DFTB+ crashes after output (expected)
                
                energies.append(energy)
                charges.append(all_frag_charges)
            
            # Save CDFTB results
            energy_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for e in energies:
                energy_parts.append(f"{e:18.10f}")
            energy_file.write("  ".join(energy_parts) + "\n")
            energy_file.flush()
            
            charge_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for qs in charges:
                for q in qs:
                    charge_parts.append(f"{q:+12.6f}")
            charge_file.write("  ".join(charge_parts) + "\n")
            charge_file.flush()
            
            # Save retry info
            retry_parts = [f"{frame_id:5d}", f"{time_val:12.3f}"]
            for info in retry_infos:
                # Use '-' for empty string (first attempt success)
                retry_parts.append(f"{info if info else '-':>20s}")
            retry_log_file.write("  ".join(retry_parts) + "\n")
            retry_log_file.flush()
            
            # Compute CDFTB-CI if enabled and CDFTB was successful
            if config.ci_enabled and cdftb_success and len(config.fragments) == 2:
                E_A = energies[0]
                E_B = energies[1]
                
                ham, J_direct, J_lowdin, ci_error = compute_cdftbci_for_frame(
                    work_dir,
                    config.fragments[0].name,
                    config.fragments[1].name,
                    E_A, E_B,
                    config.ci_N_A, config.ci_N_B
                )
                
                if ci_error:
                    print(f"  CDFTB-CI: ERROR - {ci_error}")
                    # Write NaN values
                    ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {'nan':>12s}  "
                                  f"{'nan':>16s}  {'nan':>16s}  {'nan':>10s}\n")
                    ci_sub_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 17) + "\n")
                else:
                    # Solve eigenvalue problem
                    eigenvalues, _ = solve_cdftbci_unrestricted(ham.H, ham.S)
                    
                    J_direct_meV = J_direct * 27211.386
                    J_lowdin_meV = J_lowdin * 27211.386
                    dE_eV = (eigenvalues[1] - eigenvalues[0]) * 27.211386
                    
                    print(f"  CDFTB-CI: S_AB={ham.S_AB:.6f}, J={J_lowdin_meV:.2f} meV, ΔE={dE_eV:.4f} eV")
                    
                    # Write CI results
                    ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {J_lowdin_meV:12.4f}  "
                                  f"{eigenvalues[0]:16.10f}  {eigenvalues[1]:16.10f}  {dE_eV:10.6f}\n")
                    ci_file.flush()
                    
                    # Write CI sub values
                    ci_sub_file.write(
                        f"{frame_id:5d}  {time_val:8.2f}  {E_A:14.10f}  {E_B:14.10f}  "
                        f"{ham.H_AB:14.10f}  {J_direct_meV:12.4f}  "
                        f"{ham.V_A:14.10f}  {ham.V_B:14.10f}  "
                        f"{ham.N_A:6.1f}  {ham.N_B:6.1f}  {ham.S_AB:14.10f}  "
                        f"{ham.S_AB_alpha:14.10f}  {ham.S_AB_beta:14.10f}  "
                        f"{ham.W_BA:14.10f}  {ham.W_BA_alpha:14.10f}  {ham.W_BA_beta:14.10f}  "
                        f"{ham.W_AB:14.10f}  {ham.W_AB_alpha:14.10f}  {ham.W_AB_beta:14.10f}\n"
                    )
                    ci_sub_file.flush()
            
            # Compute spin populations for both constraint states
            if cdftb_success and len(config.fragments) == 2:
                # Get number of atoms per fragment from atom_types
                n_atoms = len(atom_types)
                n_atoms_half = n_atoms // 2
                
                # State A: constraint on fragment 1
                spin_A, spin_A_error = compute_spin_for_fragment(
                    work_dir / config.fragments[0].name,
                    n_atoms_half
                )
                
                # State B: constraint on fragment 2
                spin_B, spin_B_error = compute_spin_for_fragment(
                    work_dir / config.fragments[1].name,
                    n_atoms_half
                )
                
                if spin_A_error or spin_B_error:
                    if spin_A_error:
                        print(f"  Spin (state A): ERROR - {spin_A_error[:100]}")
                    if spin_B_error:
                        print(f"  Spin (state B): ERROR - {spin_B_error[:100]}")
                    spin_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 12) + "\n")
                else:
                    # Print spin populations
                    print(f"  Spin (state A): frag1={spin_A['frag1']['spin_density']:+.4f}, "
                          f"frag2={spin_A['frag2']['spin_density']:+.4f}")
                    print(f"  Spin (state B): frag1={spin_B['frag1']['spin_density']:+.4f}, "
                          f"frag2={spin_B['frag2']['spin_density']:+.4f}")
                    
                    # Write spin values
                    spin_file.write(
                        f"{frame_id:5d}  {time_val:8.2f}  "
                        f"{spin_A['frag1']['n_alpha_frag']:14.8f}  "
                        f"{spin_A['frag1']['n_beta_frag']:14.8f}  "
                        f"{spin_A['frag1']['spin_density']:14.8f}  "
                        f"{spin_A['frag2']['n_alpha_frag']:14.8f}  "
                        f"{spin_A['frag2']['n_beta_frag']:14.8f}  "
                        f"{spin_A['frag2']['spin_density']:14.8f}  "
                        f"{spin_B['frag1']['n_alpha_frag']:14.8f}  "
                        f"{spin_B['frag1']['n_beta_frag']:14.8f}  "
                        f"{spin_B['frag1']['spin_density']:14.8f}  "
                        f"{spin_B['frag2']['n_alpha_frag']:14.8f}  "
                        f"{spin_B['frag2']['n_beta_frag']:14.8f}  "
                        f"{spin_B['frag2']['spin_density']:14.8f}\n"
                    )
                    spin_file.flush()
            
            processed_count += 1
            print()
    
    finally:
        energy_file.close()
        charge_file.close()
        if ci_file:
            ci_file.close()
        if ci_sub_file:
            ci_sub_file.close()
        if spin_file:
            spin_file.close()
        if retry_log_file:
            retry_log_file.close()
        
        print("=" * 70)
        print(f"Energies saved to {config.energy_file}")
        print(f"Charges saved to {config.charge_file}")
        print(f"Spin populations saved to {config.spin_output_file}")
        print(f"Retry log saved to {config.retry_log_file}")
        if config.ci_enabled:
            print(f"CDFTB-CI results saved to {config.ci_output_file}")
            print(f"CDFTB-CI sub values saved to {config.ci_sub_file}")
        print("=" * 70)
