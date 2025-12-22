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
# Main: Test the module
# =============================================================================
if __name__ == "__main__":
    import numpy as np
    from pathlib import Path
    from cdftb_result_reader import (
        load_spin_polarized_calculation,
        get_orbital_info_from_eigenvec,
        get_atom_orbital_map,
        build_fragment_weight_matrix,
    )
    
    output_dir = Path(__file__).parent / "output"
    
    print("=" * 70)
    print("CDFTB-CI Hamiltonian Construction Test (Unrestricted)")
    print("=" * 70)
    
    # Read energies from energies.dat
    energies_file = output_dir / "energies.dat"
    if not energies_file.exists():
        print(f"Energies file not found: {energies_file}")
        exit(1)
    
    energies_data = np.loadtxt(energies_file, comments='#')
    if energies_data.ndim == 1:
        energies_data = energies_data.reshape(1, -1)
    
    frame_id = int(energies_data[0, 0])
    E_A = energies_data[0, 2]
    E_B = energies_data[0, 3]
    
    # Fragment directories
    frag1_dir = output_dir / f"frame_{frame_id:05d}" / "fragment1"
    frag2_dir = output_dir / f"frame_{frame_id:05d}" / "fragment2"
    
    print(f"\nFrame {frame_id}:")
    print(f"  Fragment 1 (State A): E_A = {E_A:.6f} Ha = {E_A * 27.211386:.4f} eV")
    print(f"  Fragment 2 (State B): E_B = {E_B:.6f} Ha = {E_B * 27.211386:.4f} eV")
    
    # Load spin-polarized data
    print("\nLoading spin-polarized orbital data...")
    
    try:
        data_A = load_spin_polarized_calculation(frag1_dir)
        data_B = load_spin_polarized_calculation(frag2_dir)
        
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
        
        print(f"  State A: n_α = {orb_A.n_alpha}, n_β = {orb_A.n_beta}")
        print(f"  State B: n_α = {orb_B.n_alpha}, n_β = {orb_B.n_beta}")
        
        # Read AO overlap
        S_AO = data_A['S']
        n_orbitals = data_A['n_orbitals']
        print(f"  Number of AOs: {n_orbitals}")
        
        # Build weight matrices
        eigenvec_file = frag1_dir / "eigenvec.out"
        orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
        atom_to_orbitals = get_atom_orbital_map(orbital_info)
        n_atoms = len(atom_to_orbitals)
        n_atoms_half = n_atoms // 2
        
        frag_A_atoms = list(range(1, n_atoms_half + 1))
        frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
        
        w_A = build_fragment_weight_matrix(S_AO, frag_A_atoms, n_atoms,
                                           atom_to_orbitals=atom_to_orbitals)
        w_B = build_fragment_weight_matrix(S_AO, frag_B_atoms, n_atoms,
                                           atom_to_orbitals=atom_to_orbitals)
        
        # Get constraint potentials
        V_A = data_A['Vc']
        V_B = data_B['Vc']
        
        # Target populations
        N_A = 101  # Target for pentacene cation
        N_B = 101
        
        print(f"\n  V_A = {V_A:.6f} Ha")
        print(f"  V_B = {V_B:.6f} Ha")
        
        # Build Hamiltonian
        print("\nBuilding CDFTB-CI Hamiltonian (unrestricted)...")
        
        ham = build_cdftbci_hamiltonian_unrestricted(
            orb_A, orb_B, S_AO, w_A, w_B,
            E_A, E_B, V_A, V_B, N_A, N_B
        )
        
        print(f"\nResults:")
        print(f"  State overlap S_AB = {ham.S_AB:.6f}")
        print(f"    S_AB_α = {ham.S_AB_alpha:.6f}")
        print(f"    S_AB_β = {ham.S_AB_beta:.6f}")
        print(f"  Weight overlap W_BA = {ham.W_BA:.6f}")
        print(f"    W_BA_α = {ham.W_BA_alpha:.6f}")
        print(f"    W_BA_β = {ham.W_BA_beta:.6f}")
        print(f"  Weight overlap W_AB = {ham.W_AB:.6f}")
        print(f"  Coupling H_AB = {ham.H_AB:.6f} Ha")
        print(f"           H_AB = {ham.H_AB * 27.211386:.4f} eV")
        
        print(f"\n  Hamiltonian matrix H:")
        print(f"    [{ham.H[0,0]:.6f}  {ham.H[0,1]:.6f}]")
        print(f"    [{ham.H[1,0]:.6f}  {ham.H[1,1]:.6f}]")
        
        print(f"\n  Overlap matrix S:")
        print(f"    [{ham.S[0,0]:.6f}  {ham.S[0,1]:.6f}]")
        print(f"    [{ham.S[1,0]:.6f}  {ham.S[1,1]:.6f}]")
        
        # Solve eigenvalue problem
        eigenvalues, eigenvectors = solve_cdftbci_unrestricted(ham.H, ham.S)
        
        print(f"\n  Adiabatic state energies:")
        print(f"    E_1 = {eigenvalues[0]:.6f} Ha = {eigenvalues[0] * 27.211386:.4f} eV")
        print(f"    E_2 = {eigenvalues[1]:.6f} Ha = {eigenvalues[1] * 27.211386:.4f} eV")
        print(f"    ΔE = {(eigenvalues[1] - eigenvalues[0]):.6f} Ha = {(eigenvalues[1] - eigenvalues[0]) * 27.211386:.4f} eV")
        
        # Compute transfer integral
        J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
        
        print(f"\n  Transfer integral (direct method):")
        print(f"    J = {J_direct:.6f} Ha = {J_direct * 27211.386:.2f} meV")
        
        # Compare α and β HOMO-HOMO overlaps
        print("\n  HOMO-HOMO overlaps (single orbital):")
        homo_idx_alpha = orb_A.n_alpha - 1
        homo_idx_beta = orb_A.n_beta - 1
        
        # Single orbital overlaps
        C_A_homo_alpha = orb_A.C_alpha[:, homo_idx_alpha]
        C_B_homo_alpha = orb_B.C_alpha[:, homo_idx_alpha]
        S_homo_alpha = C_B_homo_alpha @ S_AO @ C_A_homo_alpha
        
        C_A_homo_beta = orb_A.C_beta[:, homo_idx_beta]
        C_B_homo_beta = orb_B.C_beta[:, homo_idx_beta]
        S_homo_beta = C_B_homo_beta @ S_AO @ C_A_homo_beta
        
        print(f"    S_HOMO^α = {S_homo_alpha:.6f}")
        print(f"    S_HOMO^β = {S_homo_beta:.6f}")
        
        # Energy eigenvalues
        print(f"\n  HOMO/LUMO energies:")
        print(f"    State A: ε_HOMO^α = {data_A['eig_alpha'][homo_idx_alpha]:.4f} eV, "
              f"ε_HOMO^β = {data_A['eig_beta'][homo_idx_beta]:.4f} eV")
        print(f"    State B: ε_HOMO^α = {data_B['eig_alpha'][homo_idx_alpha]:.4f} eV, "
              f"ε_HOMO^β = {data_B['eig_beta'][homo_idx_beta]:.4f} eV")
        
    except Exception as e:
        print(f"\nError: {e}")
        print("\nNote: This module requires spin-polarized DFTB+ calculations.")
        print("Run DFTB+ with SpinPolarisation enabled and check output files.")
        import traceback
        traceback.print_exc()


def run_cdftbci_analysis(config_path: Path) -> None:
    """
    Run CDFTB-CI analysis on existing CDFTB output.
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file.
    """
    import yaml
    from .cdftb_result_reader import (
        load_spin_polarized_calculation,
        get_orbital_info_from_eigenvec,
        get_atom_orbital_map,
        build_fragment_weight_matrix,
    )
    
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Get output directory
    output_cfg = config.get('output', {})
    output_dir = base_dir / output_cfg.get('directory', 'output')
    
    # Get fragment info
    fragments = config.get('fragments', [])
    if len(fragments) != 2:
        raise ValueError("CDFTB-CI requires exactly 2 fragments")
    
    frag1_name = fragments[0]['name']
    frag2_name = fragments[1]['name']
    
    print("=" * 70)
    print("CDFTB-CI Analysis")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Output directory: {output_dir}")
    print(f"Fragment A: {frag1_name}")
    print(f"Fragment B: {frag2_name}")
    print("=" * 70)
    
    # Read energies from energies.dat
    energies_file = output_dir / "energies.dat"
    if not energies_file.exists():
        raise FileNotFoundError(
            f"Energies file not found: {energies_file}\n"
            "Please run run_cdftb.py first to generate CDFTB output."
        )
    
    energies_data = np.loadtxt(energies_file, comments='#')
    if energies_data.ndim == 1:
        energies_data = energies_data.reshape(1, -1)
    
    n_frames = energies_data.shape[0]
    print(f"\nFound {n_frames} frames in energies.dat")
    
    # Prepare output file
    ci_output_file = output_dir / "cdftbci.dat"
    
    with open(ci_output_file, 'w') as f:
        f.write("# CDFTB-CI Results\n")
        f.write("# Frame  Time(fs)      E_A(Ha)      E_B(Ha)        V_A(Ha)        V_B(Ha)        S_AB        "
                "H_AB(Ha)        J_direct(meV)   J_lowdin(meV)        E1(Ha)        E2(Ha)        dE(eV)\n")
    
    # Process each frame
    for i in range(n_frames):
        frame_id = int(energies_data[i, 0])
        time_fs = energies_data[i, 1]
        E_A = energies_data[i, 2]
        E_B = energies_data[i, 3]
        
        print(f"\nFrame {frame_id:05d} (t = {time_fs:.1f} fs)")
        
        # Fragment directories
        frag1_dir = output_dir / f"frame_{frame_id:05d}" / frag1_name
        frag2_dir = output_dir / f"frame_{frame_id:05d}" / frag2_name
        
        if not frag1_dir.exists() or not frag2_dir.exists():
            print(f"  WARNING: Fragment directories not found, skipping frame")
            continue
        
        try:
            # Load spin-polarized data
            data_A = load_spin_polarized_calculation(frag1_dir)
            data_B = load_spin_polarized_calculation(frag2_dir)
            
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
            
            # Read AO overlap
            S_AO = data_A['S']
            n_orbitals = data_A['n_orbitals']
            
            # Build weight matrices
            eigenvec_file = frag1_dir / "eigenvec.out"
            orbital_info = get_orbital_info_from_eigenvec(eigenvec_file)
            atom_to_orbitals = get_atom_orbital_map(orbital_info)
            n_atoms = len(atom_to_orbitals)
            n_atoms_half = n_atoms // 2
            
            frag_A_atoms = list(range(1, n_atoms_half + 1))
            frag_B_atoms = list(range(n_atoms_half + 1, n_atoms + 1))
            
            w_A = build_fragment_weight_matrix(S_AO, frag_A_atoms, n_atoms,
                                               atom_to_orbitals=atom_to_orbitals)
            w_B = build_fragment_weight_matrix(S_AO, frag_B_atoms, n_atoms,
                                               atom_to_orbitals=atom_to_orbitals)
            
            # Get constraint potentials
            V_A = data_A['Vc']
            V_B = data_B['Vc']
            
            # Target populations (Mulliken population for constrained fragment)
            # Read from config or use default
            ci_cfg = config.get('cdftb_ci', {})
            N_A = ci_cfg.get('N_A', 101)  # Default for pentacene cation
            N_B = ci_cfg.get('N_B', 101)
            
            # Debug output
            print(f"  n_alpha_A={orb_A.n_alpha}, n_beta_A={orb_A.n_beta}")
            print(f"  V_A={V_A:.6f} Ha, V_B={V_B:.6f} Ha")
            print(f"  N_A={N_A}, N_B={N_B}")
            
            # Build Hamiltonian
            ham = build_cdftbci_hamiltonian_unrestricted(
                orb_A, orb_B, S_AO, w_A, w_B,
                E_A, E_B, V_A, V_B, N_A, N_B
            )
            
            # Solve eigenvalue problem
            eigenvalues, eigenvectors = solve_cdftbci_unrestricted(ham.H, ham.S)
            
            # Compute transfer integrals
            J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
            J_lowdin = compute_transfer_integral_unrestricted(ham.H, ham.S, method="lowdin")
            
            # Convert to meV
            J_direct_meV = J_direct * 27211.386
            J_lowdin_meV = J_lowdin * 27211.386
            dE_eV = (eigenvalues[1] - eigenvalues[0]) * 27.211386
            
            print(f"  S_AB = {ham.S_AB:.6f} (α: {ham.S_AB_alpha:.6f}, β: {ham.S_AB_beta:.6f})")
            print(f"  W_BA = {ham.W_BA:.6f}, W_AB = {ham.W_AB:.6f}")
            print(f"  H_AB = {ham.H_AB:.6f} Ha")
            print(f"  J (direct) = {J_direct_meV:.2f} meV")
            print(f"  J (Löwdin) = {J_lowdin_meV:.2f} meV")
            print(f"  E1 = {eigenvalues[0]:.6f} Ha, E2 = {eigenvalues[1]:.6f} Ha")
            print(f"  ΔE = {dE_eV:.4f} eV")
            
            # Write results
            with open(ci_output_file, 'a') as f:
                f.write(f"{frame_id:5d}  {time_fs:8.2f}  {E_A:14.10f}  {E_B:14.10f}  "
                        f"{V_A:12.10f}  {V_B:12.10f}  {ham.S_AB:12.10f}  {ham.H_AB:14.10f}  "
                        f"{J_direct_meV:12.10f}  {J_lowdin_meV:12.10f}  "
                        f"{eigenvalues[0]:14.10f}  {eigenvalues[1]:14.10f}  {dE_eV:12.6f}\n")
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "=" * 70)
    print(f"Results saved to {ci_output_file}")
    print("=" * 70)
