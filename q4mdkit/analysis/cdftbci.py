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
    V_A: float = 0.0  # Charge constraint potential for state A (V_N,A)
    V_B: float = 0.0  # Charge constraint potential for state B (V_N,B)
    N_A: float = 0.0  # Target population for constraint A
    N_B: float = 0.0  # Target population for constraint B
    # Spin constraint fields
    V_M_A: float = 0.0  # Spin constraint potential for state A
    V_M_B: float = 0.0  # Spin constraint potential for state B
    M_A: float = 0.0  # Target spin for constraint A
    M_B: float = 0.0  # Target spin for constraint B
    W_M_BA: float = 0.0  # Spin weight overlap ⟨Φ^B|w^A_spin|Φ^A⟩
    W_M_AB: float = 0.0  # Spin weight overlap ⟨Φ^A|w^B_spin|Φ^B⟩
    
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
        
        # Check if spin constraint is used
        has_spin_constraint = (self.V_M_A != 0.0 or self.V_M_B != 0.0 or 
                               self.M_A != 0.0 or self.M_B != 0.0)
        
        with open(filepath, 'a') as f:
            if write_header:
                f.write("# CDFTB-CI Hamiltonian Parameters\n")
                if has_spin_constraint:
                    f.write("# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_N_A(Ha)  V_N_B(Ha)  N_A  N_B  "
                            "V_M_A(Ha)  V_M_B(Ha)  M_A  M_B  "
                            "S_AB  S_AB_alpha  S_AB_beta  W_BA  W_AB  W_M_BA  W_M_AB  H_AB(Ha)  J(meV)\n")
                else:
                    f.write("# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_A(Ha)  V_B(Ha)  N_A  N_B  "
                            "S_AB  S_AB_alpha  S_AB_beta  W_BA  W_AB  H_AB(Ha)  J(meV)\n")
            
            # Compute J
            J = self.H_AB - self.S_AB * (self.E_A + self.E_B) / 2
            J_meV = J * 27211.386
            
            if has_spin_constraint:
                f.write(f"{frame_id:6d}  {time_fs:8.2f}  {self.E_A:14.10f}  {self.E_B:14.10f}  "
                        f"{self.V_A:12.8f}  {self.V_B:12.8f}  {self.N_A:6.1f}  {self.N_B:6.1f}  "
                        f"{self.V_M_A:12.8f}  {self.V_M_B:12.8f}  {self.M_A:6.1f}  {self.M_B:6.1f}  "
                        f"{self.S_AB:12.8f}  {self.S_AB_alpha:12.8f}  {self.S_AB_beta:12.8f}  "
                        f"{self.W_BA:12.8f}  {self.W_AB:12.8f}  {self.W_M_BA:12.8f}  {self.W_M_AB:12.8f}  "
                        f"{self.H_AB:14.10f}  {J_meV:10.4f}\n")
            else:
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
    W_AB: float,
    V_M_A: float = 0.0,
    V_M_B: float = 0.0,
    M_A: float = 0.0,
    M_B: float = 0.0,
    W_M_BA: float = 0.0,
    W_M_AB: float = 0.0
) -> float:
    """
    Compute the coupling element H_AB with optional spin constraint.
    
    Without spin constraint:
        H_AB = 1/2 (E_A + E_B + N_A V_N,A + N_B V_N,B) S_AB 
               - 1/2 (V_N,A W_N,BA + V_N,B W_N,AB)
    
    With spin constraint:
        H_AB = 1/2 (E_A + E_B + N_A V_N,A + M_A V_M,A + N_B V_N,B + M_B V_M,B) S_AB 
               - 1/2 (V_N,A W_N,BA + V_M,A W_M,BA + V_N,B W_N,AB + V_M,B W_M,AB)
    
    Parameters
    ----------
    E_A : float
        Total energy of state A (CDFTB energy).
    E_B : float
        Total energy of state B (CDFTB energy).
    V_A : float
        Charge constraint potential for state A (V_N,A).
    V_B : float
        Charge constraint potential for state B (V_N,B).
    N_A : float
        Target population for constraint A.
    N_B : float
        Target population for constraint B.
    S_AB : float
        State overlap ⟨Φ^B|Φ^A⟩.
    W_BA : float
        Charge weight overlap ⟨Φ^B|w^A_N|Φ^A⟩.
    W_AB : float
        Charge weight overlap ⟨Φ^A|w^B_N|Φ^B⟩.
    V_M_A : float, optional
        Spin constraint potential for state A (V_M,A). Default is 0.0.
    V_M_B : float, optional
        Spin constraint potential for state B (V_M,B). Default is 0.0.
    M_A : float, optional
        Target spin population for constraint A. Default is 0.0.
    M_B : float, optional
        Target spin population for constraint B. Default is 0.0.
    W_M_BA : float, optional
        Spin weight overlap ⟨Φ^B|w^A_M|Φ^A⟩. Default is 0.0.
    W_M_AB : float, optional
        Spin weight overlap ⟨Φ^A|w^B_M|Φ^B⟩. Default is 0.0.
        
    Returns
    -------
    H_AB : float
        Coupling element.
    """
    # Charge constraint terms
    term1 = 0.5 * (E_A + E_B + N_A * V_A + N_B * V_B) * S_AB
    term2 = 0.5 * (V_A * W_BA + V_B * W_AB)
    
    # Spin constraint terms (if present)
    term1 += 0.5 * (M_A * V_M_A + M_B * V_M_B) * S_AB
    term2 += 0.5 * (V_M_A * W_M_BA + V_M_B * W_M_AB)
    
    H_AB = term1 - term2
    
    return H_AB


def compute_spin_weight_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray,
    Omega_BA_alpha: np.ndarray,
    Omega_BA_beta: np.ndarray
) -> float:
    """
    Compute spin weight overlap W_M,BA = ⟨Φ^B|w^A_spin|Φ^A⟩.
    
    For spin constraint, the weight overlap is:
        W_M,BA = S_AB × [Tr(O^{-1,α}_BA Ω^α_BA) - Tr(O^{-1,β}_BA Ω^β_BA)]
    
    Note the MINUS sign for beta, unlike charge weight overlap which uses plus.
    This is because spin population = N_alpha - N_beta.
    
    Parameters
    ----------
    O_BA_alpha : np.ndarray
        α MO overlap matrix (n_occ × n_occ).
    O_BA_beta : np.ndarray
        β MO overlap matrix (n_occ × n_occ).
    Omega_BA_alpha : np.ndarray
        α Ω matrix for fragment weight.
    Omega_BA_beta : np.ndarray
        β Ω matrix for fragment weight.
        
    Returns
    -------
    W_M_BA : float
        Spin weight overlap.
    """
    # Compute state overlap
    S_AB, S_AB_alpha, S_AB_beta = compute_state_overlap_unrestricted(
        O_BA_alpha, O_BA_beta
    )
    
    # α contribution
    if O_BA_alpha.shape[0] > 0:
        O_inv_alpha = np.linalg.inv(O_BA_alpha)
        trace_alpha = np.trace(O_inv_alpha @ Omega_BA_alpha)
    else:
        trace_alpha = 0.0
    
    # β contribution (negative for spin)
    if O_BA_beta.shape[0] > 0:
        O_inv_beta = np.linalg.inv(O_BA_beta)
        trace_beta = np.trace(O_inv_beta @ Omega_BA_beta)
    else:
        trace_beta = 0.0
    
    # Spin weight overlap: W_M_BA = S_AB × (trace_α - trace_β)
    W_M_BA = S_AB * (trace_alpha - trace_beta)
    
    return W_M_BA


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
    N_B: float,
    V_M_A: float = 0.0,
    V_M_B: float = 0.0,
    M_A: float = 0.0,
    M_B: float = 0.0
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
        Charge constraint potential for state A (V_N,A).
    V_B : float
        Charge constraint potential for state B (V_N,B).
    N_A : float
        Target population for constraint A.
    N_B : float
        Target population for constraint B.
    V_M_A : float, optional
        Spin constraint potential for state A. Default is 0.0.
    V_M_B : float, optional
        Spin constraint potential for state B. Default is 0.0.
    M_A : float, optional
        Target spin population for constraint A. Default is 0.0.
    M_B : float, optional
        Target spin population for constraint B. Default is 0.0.
        
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
    
    # Compute charge weight overlaps
    W_BA, W_BA_alpha, W_BA_beta = compute_weight_overlap_unrestricted(
        O_BA_alpha, O_BA_beta, Omega_BA_alpha, Omega_BA_beta
    )
    W_AB, W_AB_alpha, W_AB_beta = compute_weight_overlap_unrestricted(
        O_AB_alpha, O_AB_beta, Omega_AB_alpha, Omega_AB_beta
    )
    
    # Compute spin weight overlaps (only if spin constraint is used)
    W_M_BA = 0.0
    W_M_AB = 0.0
    if V_M_A != 0.0 or V_M_B != 0.0 or M_A != 0.0 or M_B != 0.0:
        W_M_BA = compute_spin_weight_overlap_unrestricted(
            O_BA_alpha, O_BA_beta, Omega_BA_alpha, Omega_BA_beta
        )
        W_M_AB = compute_spin_weight_overlap_unrestricted(
            O_AB_alpha, O_AB_beta, Omega_AB_alpha, Omega_AB_beta
        )
    
    # Compute coupling element
    H_AB = compute_coupling_element_unrestricted(
        E_A, E_B, V_A, V_B, N_A, N_B, S_AB, W_BA, W_AB,
        V_M_A, V_M_B, M_A, M_B, W_M_BA, W_M_AB
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
        N_B=N_B,
        V_M_A=V_M_A,
        V_M_B=V_M_B,
        M_A=M_A,
        M_B=M_B,
        W_M_BA=W_M_BA,
        W_M_AB=W_M_AB
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


def effective_state_overlap_metric(
    D: float,
    n_alpha: int,
    n_beta: int,
    logabs_D: Optional[float] = None,
) -> float:
    """
    Return the per-occupied-orbital effective overlap |D|^(1 / N_occ).

    Here D is the many-electron state overlap determinant product and
    N_occ = n_alpha + n_beta. This rescales the determinant to a quantity that
    is easier to interpret across systems of different sizes.
    """
    n_occ = n_alpha + n_beta
    if n_occ <= 0:
        return float("nan")
    if logabs_D is not None and np.isfinite(logabs_D):
        return float(np.exp(logabs_D / n_occ))
    if not np.isfinite(D):
        return float("nan")
    return float(abs(D) ** (1.0 / n_occ))


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
import re
import shutil
import yaml
from dataclasses import field
from typing import List, Dict, Any, Tuple as _Tuple


_SCC_LINE_RE = re.compile(
    r"^\s*(\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s*$",
    re.M,
)


def _read_scc_from_detailed(detailed_path: Path) -> Optional[_Tuple[int, float]]:
    """Extract last (iSCC, SCC_error) from a DFTB+ detailed.out file.

    Returns None if file missing or no SCC convergence line found.
    """
    p = Path(detailed_path)
    if not p.exists():
        return None
    try:
        txt = p.read_text(errors="ignore")
    except Exception:
        return None
    matches = _SCC_LINE_RE.findall(txt)
    if not matches:
        return None
    iSCC, _, _, err = matches[-1]
    try:
        return int(iSCC), float(err)
    except ValueError:
        return None

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

from .phase_tracking import StatePhaseTracker, choose_phase_continuity_override
from .odin_overlap import compute_cross_overlap_odin
from .overlap_dump import save_frame_overlap_matrices_binary


# Default retry strategy used when no `scc.retry` is given in the YAML.
# Each entry is a dict consumed by `run_cdftb_with_retry`. Keys:
#   label                  : tag written to retry/SCC logs
#   use_initial_charges    : None -> follow per-frame policy
#                            (config.use_previous_charges + frame index);
#                            True/False -> override
#   disable_constraint     : drop ElectronicConstraints block
#   mixing_parameter       : override Broyden MixingParameter
#   max_scc_iterations     : override MaxSCCIterations
DEFAULT_RETRY_ATTEMPTS: List[Dict[str, Any]] = [
    {"label": "first", "use_initial_charges": None},
    {"label": "mix0.01", "use_initial_charges": True, "mixing_parameter": 0.01, "max_scc_iterations": 500},
    {"label": "no_init_charges", "use_initial_charges": False, "mixing_parameter": 0.01, "max_scc_iterations": 500},
    #{"label": "mix0.02", "use_initial_charges": False, "mixing_parameter": 0.02, "max_scc_iterations": 1500},
]


@dataclass
class CDFTBCIConfig(CDFTBConfig):
    """Configuration for CDFTB-CI analysis."""
    # CI settings
    ci_enabled: bool = True
    ci_N_A: float = 101.0  # Target charge population for constraint A
    ci_N_B: float = 101.0  # Target charge population for constraint B
    ci_M_A: float = 0.0  # Target spin population for constraint A (0 if no spin constraint)
    ci_M_B: float = 0.0  # Target spin population for constraint B (0 if no spin constraint)
    
    # Output mode
    output_mode: str = "online_ci"  # "online_ci" or "store_frames"
    work_directory: str = "work"
    
    # CI output files
    ci_output_file: Optional[Path] = None
    ci_sub_file: Optional[Path] = None
    overlap_matrices_binary_enabled: bool = False
    overlap_matrices_binary_dir: Optional[Path] = None
    
    # Spin output file
    spin_output_file: Optional[Path] = None
    
    # Retry log file
    retry_log_file: Optional[Path] = None
    
    # Use previous frame's charges as initial guess
    use_previous_charges: bool = False

    # SCC convergence logging (iSCC, SCC error per fragment per frame,
    # including every retry attempt).
    scc_log_enabled: bool = False
    scc_log_file: Optional[Path] = None

    # Retry strategy: ordered list of attempt dicts (see DEFAULT_RETRY_ATTEMPTS
    # for the schema). When None, the default sequence is used.
    retry_attempts: Optional[List[Dict[str, Any]]] = None

    # Phase (gauge) tracking of diabatic states across MD frames.
    # When enabled, the sign s_f(t_n) of each charge-localized diabatic
    # state is chosen so that <Phi_f^corr(t_{n-1}) | Phi_f(t_n)> > 0, and
    # the coupling and state overlap are corrected as
    #   H_AB^corr = s_A s_B H_AB,  S_AB^corr = s_A s_B S_AB.
    phase_tracking_enabled: bool = True
    # Approximation for the cross-geometry AO overlap S^{n-1,n}:
    #   "current"  -> S(t_n)                    (default, robust)
    #   "previous" -> S(t_{n-1})
    #   "midpoint" -> 0.5 [S(t_{n-1}) + S(t_n)]
    phase_tracking_cross_overlap_mode: str = "current"
    # Output file for the phase-tracking diagnostics.
    phase_output_file: Optional[Path] = None
    # Warn when |D_raw|^(1/N_occ) < this value.
    phase_warn_low_overlap: float = 0.5
    # Warn when the minimum occupied-space singular value gets too small.
    phase_warn_low_sigma_min: float = 0.9
    # If the latest occupied-overlap singular value falls below this threshold,
    # try older stored references before accepting the phase sign.
    phase_sigma_accept_threshold: float = 0.0
    # If True, remove rigid translation and rotation from QM geometries before
    # exact ODIN cross-overlap evaluation.
    phase_remove_translation_rotation_before_overlap: bool = False
    # If True, use singular-value thresholds to trigger lookback and frame
    # invalidation. If False, always accept the latest reference.
    phase_sigma_filtering_enabled: bool = True
    # Apply the SVD-based corresponding-orbital alignment before committing a
    # transported occupied-space reference.
    phase_corresponding_orbital_alignment: bool = True
    # If True, a frame is invalidated when none of the stored references reaches
    # the sigma threshold. Older references are still tried first.
    phase_invalidate_low_primary_sigma: bool = True
    # Number of recent gauge-fixed references kept temporarily for lookback.
    phase_reference_history: int = 5
    # Warn when |D_raw|^(1/N_occ) > this value.
    phase_warn_high_overlap: float = 1.5
    # Warn when sign(H_AB_corr) flips between consecutive frames after
    # gauge correction (possible genuine diabatic crossing OR tracker miss).
    phase_warn_post_correction_flip: bool = True
    # Number of transported references used in the wavefunction-only sign vote.
    phase_vote_history: int = 1
    # Geometric decay applied to older vote references.
    phase_vote_decay: float = 0.6
    # If |vote_score| / total_vote_weight is below this value, keep the
    # previous sign and mark the frame as ambiguous.
    phase_vote_ambiguity_ratio: float = 0.15
    # When the overlap-derived gauge is weak, optionally choose the sign that
    # keeps H_AB and S_AB continuous with the previous corrected frame.
    # Disabled by default: the primary sign decision should come from
    # wavefunction overlaps, not from Hamiltonian continuity.
    phase_continuity_override_enabled: bool = False
    phase_continuity_overlap_threshold: float = 0.995
    # ODIN cross-overlap settings (used when cross_overlap_mode == "odin")
    phase_odin_executable: Optional[str] = None
    phase_odin_sk_prefix: str = ""
    phase_odin_sk_separator: str = "-"
    phase_odin_sk_suffix: str = ".skf"
    phase_odin_lmax: Dict[str, int] = field(default_factory=dict)  # {"C": 2, "H": 1}
    phase_odin_work_subdir: str = "odin_work"
    phase_odin_keep_files: bool = False
    # Optional phase-tracking driven adaptive refinement. When enabled, a low
    # occupied-overlap singular value at a coarse frame triggers temporary
    # evaluation of the intermediate frames from a finer trajectory.
    phase_adaptive_refinement_enabled: bool = False
    phase_adaptive_refinement_traj_path: Optional[Path] = None
    phase_adaptive_refinement_sigma_min_threshold: float = 0.0
    phase_adaptive_refinement_fine_dt_fs: float = 1.0
    phase_adaptive_refinement_log_file: Optional[Path] = None


def _fine_frame_indices_between(
    prev_time_fs: float,
    curr_time_fs: float,
    fine_dt_fs: float,
    *,
    include_current: bool = False,
) -> List[int]:
    """Return fine-trajectory frame indices between two times."""
    if fine_dt_fs <= 0.0:
        raise ValueError("fine_dt_fs must be positive")
    start = int(np.floor(prev_time_fs / fine_dt_fs)) + 1
    stop = int(np.floor(curr_time_fs / fine_dt_fs)) + (
        1 if include_current else 0
    )
    return list(range(start, stop))


def _phase_adaptive_refinement_indices(
    *,
    enabled: bool,
    previous_time_fs: Optional[float],
    current_time_fs: float,
    sigma_min_a: float,
    sigma_min_b: float,
    threshold: float,
    fine_dt_fs: float,
) -> List[int]:
    """Return fine-frame indices requested by phase low-sigma diagnostics."""
    if not enabled or previous_time_fs is None:
        return []
    sigma_min = min(sigma_min_a, sigma_min_b)
    if not np.isfinite(sigma_min) or sigma_min > threshold:
        return []
    return _fine_frame_indices_between(
        previous_time_fs,
        current_time_fs,
        fine_dt_fs,
        include_current=True,
    )


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
    ci_M_A = ci_cfg.get('M_A', 0.0)  # Spin constraint target for A
    ci_M_B = ci_cfg.get('M_B', 0.0)  # Spin constraint target for B

    # Phase tracking settings
    phase_cfg = data.get('phase_tracking', {})
    phase_tracking_enabled = phase_cfg.get('enabled', True)
    phase_tracking_cross_overlap_mode = phase_cfg.get('cross_overlap_mode', 'current')
    phase_warn_low_overlap = float(phase_cfg.get('warn_low_overlap', 0.5))
    phase_warn_low_sigma_min = float(phase_cfg.get('warn_low_sigma_min', 0.9))
    phase_sigma_accept_threshold = float(
        phase_cfg.get('sigma_accept_threshold', phase_warn_low_sigma_min)
    )
    phase_remove_translation_rotation_before_overlap = bool(
        phase_cfg.get('remove_translation_rotation_before_overlap', False)
    )
    phase_sigma_filtering_enabled = bool(
        phase_cfg.get('sigma_filtering_enabled', True)
    )
    phase_corresponding_orbital_alignment = bool(
        phase_cfg.get('corresponding_orbital_alignment', True)
    )
    phase_invalidate_low_primary_sigma = bool(
        phase_cfg.get('invalidate_low_primary_sigma', True)
    )
    phase_reference_history = int(phase_cfg.get('reference_history', 5))
    phase_warn_high_overlap = float(phase_cfg.get('warn_high_overlap', 1.5))
    phase_warn_post_correction_flip = bool(phase_cfg.get('warn_post_correction_flip', True))
    phase_vote_history = int(phase_cfg.get('vote_history', 1))
    phase_vote_decay = float(phase_cfg.get('vote_decay', 0.6))
    phase_vote_ambiguity_ratio = float(phase_cfg.get('vote_ambiguity_ratio', 0.15))
    phase_continuity_override_enabled = bool(
        phase_cfg.get('continuity_override_enabled', False)
    )
    phase_continuity_overlap_threshold = float(
        phase_cfg.get('continuity_overlap_threshold', 0.995)
    )
    odin_cfg = phase_cfg.get('odin', {})
    phase_odin_executable = odin_cfg.get('executable', None)
    phase_odin_sk_prefix = odin_cfg.get('sk_prefix', '')
    phase_odin_sk_separator = odin_cfg.get('sk_separator', '-')
    phase_odin_sk_suffix = odin_cfg.get('sk_suffix', '.skf')
    phase_odin_lmax = {str(k): int(v) for k, v in odin_cfg.get('lmax', {}).items()}
    phase_odin_work_subdir = odin_cfg.get('work_subdir', 'odin_work')
    phase_odin_keep_files = bool(odin_cfg.get('keep_files', False))
    phase_adaptive_cfg = phase_cfg.get('adaptive_refinement', {})
    phase_adaptive_refinement_enabled = bool(phase_adaptive_cfg.get('enabled', False))
    phase_adaptive_refinement_traj_path = None
    if phase_adaptive_cfg.get('trajectory_1fs') is not None:
        phase_adaptive_refinement_traj_path = base_dir / phase_adaptive_cfg['trajectory_1fs']
    phase_adaptive_refinement_sigma_min_threshold = float(
        phase_adaptive_cfg.get('sigma_min_threshold', phase_warn_low_sigma_min)
    )
    phase_adaptive_refinement_fine_dt_fs = float(
        phase_adaptive_cfg.get('fine_dt_fs', 1.0)
    )
    phase_adaptive_refinement_log_file = output_dir / phase_adaptive_cfg.get(
        'log_file',
        'adaptive_refinement.dat',
    )
    
    # Parse SCC settings
    scc_cfg = data.get('scc', {})
    use_previous_charges = scc_cfg.get('use_previous_charges', False)
    scc_log_enabled = bool(scc_cfg.get('log_enabled', False))
    scc_log_filename = scc_cfg.get('logfile', 'scc_error.dat')

    # Parse retry strategy. Accepts either a list of dicts, or None/missing
    # (use DEFAULT_RETRY_ATTEMPTS).
    retry_cfg = scc_cfg.get('retry', None)
    if retry_cfg is None:
        retry_attempts = None
    else:
        if not isinstance(retry_cfg, list) or len(retry_cfg) == 0:
            raise ValueError("scc.retry must be a non-empty list of attempt dicts")
        retry_attempts = []
        for i, item in enumerate(retry_cfg):
            if not isinstance(item, dict):
                raise ValueError(f"scc.retry[{i}] must be a mapping")
            retry_attempts.append({
                "label": str(item.get("label", f"attempt{i}")),
                "use_initial_charges": item.get("use_initial_charges", None),
                "disable_constraint": bool(item.get("disable_constraint", False)),
                "mixing_parameter": item.get("mixing_parameter", None),
                "max_scc_iterations": item.get("max_scc_iterations", None),
                "mixer": item.get("mixer", None),
                "optimiser": item.get("optimiser", None),
                "max_constr_iterations": item.get("max_constr_iterations", None),
            })
    
    # CI output files
    ci_output_file = output_dir / output_cfg.get('ci_file', 'cdftbci.dat')
    ci_sub_file = output_dir / output_cfg.get('ci_sub_file', 'cdftbci_sub.dat')
    overlap_matrices_binary_enabled = bool(
        output_cfg.get('overlap_matrices_binary_enabled', False)
    )
    overlap_matrices_binary_dir = output_dir / output_cfg.get(
        'overlap_matrices_binary_dir',
        'overlap_matrices',
    )
    
    # Spin output file
    spin_output_file = output_dir / output_cfg.get('spin_file', 'spin.dat')
    
    # Retry log file
    retry_log_file = output_dir / output_cfg.get('retry_log_file', 'retry_log.dat')

    # SCC convergence log file
    scc_log_file = output_dir / scc_log_filename

    # Phase tracking output file
    phase_output_file = output_dir / output_cfg.get('phase_file', 'cdftbci_phase.dat')
    
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
        ci_M_A=ci_M_A,
        ci_M_B=ci_M_B,
        output_mode=output_mode,
        work_directory=work_directory,
        ci_output_file=ci_output_file,
        ci_sub_file=ci_sub_file,
        overlap_matrices_binary_enabled=overlap_matrices_binary_enabled,
        overlap_matrices_binary_dir=overlap_matrices_binary_dir,
        spin_output_file=spin_output_file,
        retry_log_file=retry_log_file,
        use_previous_charges=use_previous_charges,
        scc_log_enabled=scc_log_enabled,
        scc_log_file=scc_log_file,
        retry_attempts=retry_attempts,
        phase_tracking_enabled=phase_tracking_enabled,
        phase_tracking_cross_overlap_mode=phase_tracking_cross_overlap_mode,
        phase_output_file=phase_output_file,
        phase_warn_low_overlap=phase_warn_low_overlap,
        phase_warn_low_sigma_min=phase_warn_low_sigma_min,
        phase_sigma_accept_threshold=phase_sigma_accept_threshold,
        phase_remove_translation_rotation_before_overlap=phase_remove_translation_rotation_before_overlap,
        phase_sigma_filtering_enabled=phase_sigma_filtering_enabled,
        phase_corresponding_orbital_alignment=phase_corresponding_orbital_alignment,
        phase_invalidate_low_primary_sigma=phase_invalidate_low_primary_sigma,
        phase_reference_history=phase_reference_history,
        phase_warn_high_overlap=phase_warn_high_overlap,
        phase_warn_post_correction_flip=phase_warn_post_correction_flip,
        phase_vote_history=phase_vote_history,
        phase_vote_decay=phase_vote_decay,
        phase_vote_ambiguity_ratio=phase_vote_ambiguity_ratio,
        phase_continuity_override_enabled=phase_continuity_override_enabled,
        phase_continuity_overlap_threshold=phase_continuity_overlap_threshold,
        phase_odin_executable=phase_odin_executable,
        phase_odin_sk_prefix=phase_odin_sk_prefix,
        phase_odin_sk_separator=phase_odin_sk_separator,
        phase_odin_sk_suffix=phase_odin_sk_suffix,
        phase_odin_lmax=phase_odin_lmax,
        phase_odin_work_subdir=phase_odin_work_subdir,
        phase_odin_keep_files=phase_odin_keep_files,
        phase_adaptive_refinement_enabled=phase_adaptive_refinement_enabled,
        phase_adaptive_refinement_traj_path=phase_adaptive_refinement_traj_path,
        phase_adaptive_refinement_sigma_min_threshold=phase_adaptive_refinement_sigma_min_threshold,
        phase_adaptive_refinement_fine_dt_fs=phase_adaptive_refinement_fine_dt_fs,
        phase_adaptive_refinement_log_file=phase_adaptive_refinement_log_file,
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
    mixing_parameter: Optional[float] = None,
    max_scc_iterations: Optional[int] = None,
    mixer: Optional[Dict[str, Any]] = None,
    optimiser: Optional[Dict[str, Any]] = None,
    max_constr_iterations: Optional[int] = None,
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
    mixing_parameter : float, optional
        Override MixingParameter inside whatever Mixer block exists in the
        template. None means use template value.
    max_scc_iterations : int, optional
        Override MaxSCCIterations. None means use template value.
    mixer : dict, optional
        Replace the entire ``Hamiltonian.DFTB.Mixer`` block. Format is the
        same nested dict that ``hsd`` consumes, e.g.

            {"Anderson": {"MixingParameter": 0.05, "Generations": 4,
                          "InitMixingParameter": 0.01}}
            {"Broyden":  {"MixingParameter": 0.05}}
            {"Simple":   {"MixingParameter": 0.05}}

        ``mixing_parameter`` (if given) is applied AFTER the swap, into the
        new mixer's parameters.
    optimiser : dict, optional
        Replace the entire ``ElectronicConstraints.Optimiser`` block. Format
        is the same nested dict that ``hsd`` consumes, e.g.

            {"LBFGS": {"Memory": 20}}
            {"FIRE":  {}}
            {"SteepestDescent": {}}

        Has no effect when ``disable_constraint`` is True.
    max_constr_iterations : int, optional
        Override ``ElectronicConstraints.MaxConstrIterations`` (the outer
        Lagrange loop cap). None means use template value.
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
    
    # Override SCC parameters if specified
    # Replace the entire Mixer block first (e.g., Broyden -> Anderson) so
    # that a subsequent mixing_parameter override targets the new mixer.
    if mixer is not None:
        data['Hamiltonian']['DFTB']['Mixer'] = dict(mixer)
    if mixing_parameter is not None:
        if 'Mixer' in data['Hamiltonian']['DFTB']:
            mixer_block = data['Hamiltonian']['DFTB']['Mixer']
            # Update MixingParameter inside whichever scheme is active.
            for scheme in ('Broyden', 'Anderson', 'Simple', 'DIIS'):
                if scheme in mixer_block and isinstance(mixer_block[scheme], dict):
                    mixer_block[scheme]['MixingParameter'] = mixing_parameter
                    break
    if max_scc_iterations is not None:
        data['Hamiltonian']['DFTB']['MaxSCCIterations'] = max_scc_iterations

    # Override constraint-loop parameters (only meaningful when constraint
    # block is still present, i.e. disable_constraint is False).
    ec = data['Hamiltonian']['DFTB'].get('ElectronicConstraints')
    if isinstance(ec, dict):
        if optimiser is not None:
            ec['Optimiser'] = dict(optimiser)
        if max_constr_iterations is not None:
            ec['MaxConstrIterations'] = max_constr_iterations
    
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
    M_A: float = 0.0,
    M_B: float = 0.0,
) -> Tuple[
    Optional[UnrestrictedCDFTBCIHamiltonian],
    Optional[float],
    Optional[float],
    Optional[str],
    Optional[UnrestrictedOrbitalData],
    Optional[UnrestrictedOrbitalData],
    Optional[np.ndarray],
]:
    """
    Compute CDFTB-CI quantities for current frame using data in work directory.
    
    Parameters
    ----------
    work_dir : Path
        Working directory containing fragment calculations.
    frag1_name : str
        Name of first fragment directory.
    frag2_name : str
        Name of second fragment directory.
    E_A : float
        Total energy of state A.
    E_B : float
        Total energy of state B.
    N_A : float
        Target charge population for constraint A.
    N_B : float
        Target charge population for constraint B.
    M_A : float, optional
        Target spin population for constraint A. Default is 0.0.
    M_B : float, optional
        Target spin population for constraint B. Default is 0.0.
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
        
        # Get constraint potentials (charge and spin)
        V_A = data_A['Vc']  # charge constraint potential V_N
        V_B = data_B['Vc']
        
        # Get spin constraint potentials if available
        V_M_A = 0.0
        V_M_B = 0.0
        M_A_target = 0.0
        M_B_target = 0.0
        
        if data_A.get('Vc_potentials') is not None:
            potentials_A = data_A['Vc_potentials']
            if potentials_A.V_M is not None:
                V_M_A = potentials_A.V_M
                M_A_target = M_A  # Use passed M_A value
        
        if data_B.get('Vc_potentials') is not None:
            potentials_B = data_B['Vc_potentials']
            if potentials_B.V_M is not None:
                V_M_B = potentials_B.V_M
                M_B_target = M_B  # Use passed M_B value
        
        # Build CDFTB-CI Hamiltonian
        ham = build_cdftbci_hamiltonian_unrestricted(
            orb_A, orb_B, S_AO, w_A, w_B,
            E_A, E_B, V_A, V_B, N_A, N_B,
            V_M_A, V_M_B, M_A_target, M_B_target
        )
        
        # Compute transfer integrals
        J_direct = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
        J_lowdin = compute_transfer_integral_unrestricted(ham.H, ham.S, method="lowdin")

        return ham, J_direct, J_lowdin, None, orb_A, orb_B, S_AO

    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        return None, None, None, error_msg, None, None, None


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
    retry_attempts: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[float, List[float], bool, str, List[_Tuple[str, Optional[int], Optional[float]]]]:
    """
    Run CDFTB calculation with automatic retry on convergence failure.

    Retry strategy is configurable via ``retry_attempts``. If not given, the
    default sequence is used (kept for backward compatibility):

    1. ``first``           : Run with ``use_initial_charges`` as passed in
    2. ``no_init_charges`` : Force ``ReadInitialCharges = No``
    3. ``mix0.05``         : ``MixingParameter=0.05``, ``MaxSCCIterations=1000``
    4. ``mix0.02``         : ``MixingParameter=0.02``, ``MaxSCCIterations=1500``

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
        Default policy for ``ReadInitialCharges`` (used when an attempt sets
        ``use_initial_charges = None`` / "auto").
    ci_enabled : bool
        Whether CDFTB-CI is enabled (kept for API compatibility).
    retry_attempts : list of dict, optional
        Ordered list of attempts. Each entry supports the keys:

        - ``label`` (str)                 : tag written to the SCC/retry logs
        - ``use_initial_charges`` (bool or None) : None = follow the
          ``use_initial_charges`` argument; True/False overrides it
        - ``disable_constraint`` (bool)   : drop ``ElectronicConstraints``
        - ``mixing_parameter`` (float or None) : override MixingParameter
          inside whichever Mixer block is active
        - ``max_scc_iterations`` (int or None) : override
          ``Hamiltonian.DFTB.MaxSCCIterations``
        - ``mixer`` (dict or None)        : replace the whole
          ``Hamiltonian.DFTB.Mixer`` block, e.g.
          ``{"Anderson": {"MixingParameter": 0.05, "Generations": 4}}``
        - ``optimiser`` (dict or None)    : replace the whole
          ``ElectronicConstraints.Optimiser`` block, e.g.
          ``{"LBFGS": {"Memory": 20}}``, ``{"FIRE": {}}``,
          ``{"SteepestDescent": {}}``
        - ``max_constr_iterations`` (int or None) : override
          ``ElectronicConstraints.MaxConstrIterations`` (outer Lagrange loop)

    Returns
    -------
    energy : float
        Energy (nan if all attempts failed).
    all_frag_charges : list of float
        Mulliken charges for all fragments (nan if failed).
    success : bool
        Whether calculation succeeded.
    retry_info : str
        ``""`` if the very first attempt succeeded, otherwise the label of the
        successful attempt, or ``"all_failed:<error>"``.
    scc_attempts : list of tuple
        List of ``(attempt_label, iSCC, scc_error, status)`` for every
        attempt that produced a detailed.out. ``iSCC`` / ``scc_error`` may be
        ``None`` if detailed.out could not be parsed. ``status`` is ``"ok"``
        for a successful attempt, otherwise the error message returned by
        the DFTB+ subprocess (e.g. ``"Process crashed with exit code 1"``,
        ``"Timeout"``).
    """
    if retry_attempts is None:
        retry_attempts = DEFAULT_RETRY_ATTEMPTS

    frag_dir = work_dir / frag.name
    scc_attempts: list = []

    def _record(label: str, err_msg: Optional[str]) -> None:
        info = _read_scc_from_detailed(frag_dir / "detailed.out")
        status = "ok" if err_msg is None else err_msg
        if info is None:
            scc_attempts.append((label, None, None, status))
        else:
            scc_attempts.append((label, info[0], info[1], status))

    last_error = "no attempts configured"
    for idx, attempt in enumerate(retry_attempts):
        label = str(attempt.get("label", f"attempt{idx}"))
        a_use_init = attempt.get("use_initial_charges", None)
        if a_use_init is None:
            read_init = bool(use_initial_charges)
        else:
            read_init = bool(a_use_init)
        disable_constraint = bool(attempt.get("disable_constraint", False))
        mixing_parameter = attempt.get("mixing_parameter", None)
        max_scc_iterations = attempt.get("max_scc_iterations", None)
        mixer = attempt.get("mixer", None)
        optimiser = attempt.get("optimiser", None)
        max_constr_iterations = attempt.get("max_constr_iterations", None)

        if idx > 0:
            print(f"    [Retry] previous attempt failed, trying '{label}'")

        frag_dir = setup_fragment_work_directory(
            work_dir, frag.name, frag.atom_range, hsd_template,
            read_initial_charges=read_init,
            disable_constraint=disable_constraint,
            mixing_parameter=mixing_parameter,
            max_scc_iterations=max_scc_iterations,
            mixer=mixer,
            optimiser=optimiser,
            max_constr_iterations=max_constr_iterations,
        )

        energy, mcharge, error = run_dftb_in_subprocess(
            qm_coords_bohr, frag_dir, write_hs=False,
            dftb_library_path=dftb_library_path,
            num_threads=num_threads,
            timeout=timeout,
        )
        _record(label, error)

        if error is None:
            retry_info = "" if idx == 0 else label
            return energy, mcharge, True, retry_info, scc_attempts

        last_error = error

    # All attempts failed
    return float("nan"), None, False, f"all_failed:{last_error}", scc_attempts


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
        if config.phase_tracking_enabled:
            print(f"  Phase tracking: ENABLED "
                  f"(cross_overlap_mode = {config.phase_tracking_cross_overlap_mode})")
        else:
            print(f"  Phase tracking: DISABLED")
    
    if config.use_previous_charges:
        print(f"\nSCC settings:")
        print(f"  Use previous charges: Yes")

    active_retry = config.retry_attempts if config.retry_attempts is not None else DEFAULT_RETRY_ATTEMPTS
    print(f"\nRetry strategy ({len(active_retry)} attempt(s)):")
    for i, a in enumerate(active_retry):
        parts = []
        uic = a.get("use_initial_charges", None)
        parts.append(f"use_init={'auto' if uic is None else uic}")
        if a.get("disable_constraint"):
            parts.append("disable_constraint=True")
        if a.get("mixing_parameter") is not None:
            parts.append(f"mixing={a['mixing_parameter']}")
        if a.get("max_scc_iterations") is not None:
            parts.append(f"maxiter={a['max_scc_iterations']}")
        if a.get("mixer") is not None:
            try:
                scheme = next(iter(a["mixer"].keys()))
            except Exception:
                scheme = "?"
            parts.append(f"mixer={scheme}")
        if a.get("optimiser") is not None:
            try:
                opt = next(iter(a["optimiser"].keys()))
            except Exception:
                opt = "?"
            parts.append(f"optimiser={opt}")
        if a.get("max_constr_iterations") is not None:
            parts.append(f"maxconstr={a['max_constr_iterations']}")
        print(f"  {i+1}. {a.get('label', f'attempt{i}')}: " + ", ".join(parts))
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
    energy_file = open(config.energy_file, "w+")
    energy_file.write(energy_header)
    
    charge_file = open(config.charge_file, "w+")
    charge_file.write(charge_header)
    
    ci_file = None
    ci_sub_file = None
    spin_file = None
    retry_log_file = None
    phase_file = None
    phase_lookback_file = None
    phase_adaptive_file = None
    scc_log_file = None
    # Pre-initialize phase-tracking warning aggregators so that the `finally`
    # cleanup block can reference them even if an exception is raised before
    # the main per-frame loop starts.
    warn_counts = {"LOW_DA": 0, "LOW_DB": 0,
                   "LOW_SIG_A": 0, "LOW_SIG_B": 0,
                   "LOOKBACK_A": 0, "LOOKBACK_B": 0,
                   "INVALID_A": 0, "INVALID_B": 0,
                   "HIGH_DA": 0, "HIGH_DB": 0,
                   "AMBIG_A": 0, "AMBIG_B": 0,
                   "FLIP_HAB": 0,
                   "CONT_A": 0, "CONT_B": 0}
    warn_frames: list = []
    
    # Open retry log file
    retry_log_file = open(config.retry_log_file, "w+")
    retry_log_file.write("# Retry Log for CDFTB Calculations\n")
    retry_log_file.write("# Retry info: '-' = success on first attempt,\n")
    retry_log_file.write("#             'no_init_charges' = succeeded without initial charges, 'all_failed:...' = all attempts failed\n")
    frag_names = [f.name for f in config.fragments]
    retry_header = "# Frame  Time(fs)  " + "  ".join([f"retry_{name}" for name in frag_names]) + "\n"
    retry_log_file.write(retry_header)

    # Open SCC convergence log file (one row per fragment per attempt)
    if config.scc_log_enabled:
        scc_log_file = open(config.scc_log_file, "w")
        scc_log_file.write("# SCC Convergence Log for CDFTB Calculations\n")
        scc_log_file.write(
            "# Note: iSCC / SCC_error are read from the LAST inner-SCC summary in\n"
            "#       detailed.out. For constrained DFTB this is the inner SCC of\n"
            "#       the last Lagrange iteration -- it can look small even when the\n"
            "#       outer constraint loop failed. Check 'status' for the real result.\n"
        )
        scc_log_file.write(
            "# Frame  Time(fs)      fragment             attempt  iSCC  SCC_error(a.u.)  status\n"
        )
    
    if config.ci_enabled:
        ci_file = open(config.ci_output_file, "w")
        ci_file.write("# CDFTB-CI Results\n")
        if config.phase_tracking_enabled:
            ci_file.write("# Columns: H_AB and S_AB are gauge-corrected (s_A*s_B applied).\n")
            ci_file.write("# J_lowdin is written as nan for INVALID phase frames.\n")
        ci_file.write("# Frame  Time(fs)   J_lowdin(meV)        E1(Ha)        E2(Ha)        dE(eV)\n")
        
        ci_sub_file = open(config.ci_sub_file, "w")
        ci_sub_file.write("# CDFTB-CI Sub Values\n")
        if config.phase_tracking_enabled:
            ci_sub_file.write(
                "# Gauge-corrected: H_AB, S_AB, S_AB_alpha, S_AB_beta, W_BA, W_AB "
                "are multiplied by s_A*s_B. Raw values are in cdftbci_phase.dat.\n"
            )
            ci_sub_file.write("# J_direct is written as nan for INVALID phase frames.\n")
        ci_sub_file.write("# Frame  Time(fs)      E_A(Ha)      E_B(Ha)        H_AB(Ha)        J_direct(meV)        "
                         "V_A(Ha)        V_B(Ha)      N_A      N_B        S_AB        "
                         "S_AB_alpha      S_AB_beta        W_BA        W_BA_alpha      W_BA_beta        "
                         "W_AB        W_AB_alpha      W_AB_beta\n")

        if config.phase_tracking_enabled:
            phase_file = open(config.phase_output_file, "w")
            phase_lookback_path = config.phase_output_file.with_name(
                config.phase_output_file.stem
                + "_lookback"
                + config.phase_output_file.suffix
            )
            phase_lookback_file = open(phase_lookback_path, "w")
            phase_file.write("# CDFTB-CI Phase (Gauge) Tracking\n")
            phase_file.write(
                "# s_A, s_B in {-1,+1}: gauge factors making "
                "<Phi^corr(t_ref)|Phi(t_n)> > 0 for the selected reference.\n"
            )
            phase_file.write(
                "# D_A_raw = selected reference overlap used for sign determination "
                "before applying s_A(t_n).\n"
            )
            phase_file.write(
                "# D_A_corr = s_A(t_n) * D_A_raw. Same for B. Low-sigma modes are "
                "not truncated; older references are tried instead. NaN means tracker reset.\n"
            )
            phase_file.write(
                "# D_A_eff = |D_A_raw|^(1/N_occ) with N_occ = n_alpha + n_beta. "
                "Same for B.\n"
            )
            phase_file.write(
                "# H_AB_raw / S_AB_raw are pre-correction; *_corr = s_A*s_B * raw.\n"
            )
            phase_file.write(
                "# cross_overlap_mode = " + config.phase_tracking_cross_overlap_mode + "\n"
            )
            phase_file.write(
                "# remove_translation_rotation_before_overlap = "
                f"{config.phase_remove_translation_rotation_before_overlap}\n"
            )
            phase_file.write(
                "# sigma_accept_threshold = "
                f"{config.phase_sigma_accept_threshold:.3f}, "
                f"reference_history = {config.phase_reference_history}, "
                "invalidate_low_primary_sigma = "
                f"{config.phase_invalidate_low_primary_sigma}\n"
            )
            phase_file.write(
                "# vote_history = "
                f"{config.phase_vote_history}, vote_decay = {config.phase_vote_decay:.3f}, "
                f"vote_ambiguity_ratio = {config.phase_vote_ambiguity_ratio:.3f}\n"
            )
            if config.phase_tracking_cross_overlap_mode == "odin":
                phase_file.write(
                    "# In 'odin' mode, the cross-geometry AO overlap is computed "
                    "externally from each stored reference frame to the current "
                    "frame using the ODIN executable.\n"
                )
            phase_file.write(
                "# WARN flags:\n"
                "#   LOW_DA / LOW_DB : D_eff below warn_low_overlap "
                f"({config.phase_warn_low_overlap:.3f}) -> sign tracking may be unreliable.\n"
                "#   LOW_SIG_A / LOW_SIG_B : min(sigma_min_alpha, sigma_min_beta) below "
                f"warn_low_sigma_min ({config.phase_warn_low_sigma_min:.3f}).\n"
                "#   LOOKBACK_A / LOOKBACK_B : latest reference had low sigma, so an "
                "older stored reference was selected.\n"
                "#   INVALID_A / INVALID_B : no stored reference reached the sigma "
                "threshold, so the frame was not kept as a future phase reference.\n"
                "#   HIGH_DA / HIGH_DB : D_eff above warn_high_overlap "
                f"({config.phase_warn_high_overlap:.3f}) -> cross-overlap approximation breaking down.\n"
                "#   AMBIG_A / AMBIG_B : no stored reference passed the sigma threshold, "
                "so the previous sign was kept.\n"
                "#   CONT_A / CONT_B : continuity override flipped state A or B to keep "
                "H_AB / S_AB on the smooth gauge branch.\n"
                "#   FLIP_HAB : sign(H_AB_corr) flipped between consecutive frames AFTER "
                "gauge correction (possible diabatic crossing or tracker miss).\n"
                "#   OK if no flag.\n"
            )
            phase_file.write(
                "# Frame  Time(fs)  s_A  s_B  gauge        D_A_raw          D_A_corr        D_A_eff         "
                "A_sigmin_a       A_sigmin_b       A_cond_a         A_cond_b         A_vote_margin    A_vote_refs  A_ref_idx  "
                "D_B_raw          D_B_corr        D_B_eff         B_sigmin_a       B_sigmin_b       B_cond_a         "
                "B_cond_b         B_vote_margin    B_vote_refs  B_ref_idx  H_AB_raw(Ha)      H_AB_corr(Ha)      "
                "S_AB_raw          S_AB_corr     flags\n"
            )
            phase_lookback_file.write("# CDFTB-CI Phase Lookback Diagnostics\n")
            phase_lookback_file.write(
                "# One row per stored reference candidate when LOOKBACK, INVALID, "
                "or AMBIG is triggered for that state.\n"
            )
            phase_lookback_file.write(
                "# accepted = 1 if min(alpha_sigma_min,beta_sigma_min) >= "
                "sigma_accept_threshold.\n"
            )
            phase_lookback_file.write(
                "# selected = first accepted candidate in this diagnostic table. "
                "INVALID frames are still not committed as phase references.\n"
            )
            phase_lookback_file.write(
                "# Frame  Time(fs)  state  ref_idx  selected  accepted  sign  "
                "D_raw             D_eff             weight            sigma_min        "
                "alpha_det         beta_det          alpha_sign  beta_sign  "
                "alpha_sigmin      beta_sigmin       alpha_sigmax      beta_sigmax       "
                "alpha_cond        beta_cond         logabs_D\n"
            )
            if config.phase_adaptive_refinement_enabled:
                if config.phase_adaptive_refinement_traj_path is None:
                    raise ValueError(
                        "phase_tracking.adaptive_refinement.trajectory_1fs is required "
                        "when phase adaptive refinement is enabled"
                    )
                phase_adaptive_file = open(config.phase_adaptive_refinement_log_file, "w")
                phase_adaptive_file.write("# CDFTB-CI Phase Adaptive Refinement\n")
                phase_adaptive_file.write(
                    "# Frame  Time(fs)  sigma_min_A  sigma_min_B  refined  fine_frames\n"
                )
    
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
        coarse_processed_count = 0

        # Phase trackers for the two charge-localized diabatic states.
        # On CDFTB failure the previous-frame MOs are simply KEPT as the
        # phase reference; the next successful frame overwrites them via
        # tracker.update(). No automatic reset -- a single bad frame does
        # not break phase continuity. If the cross-frame overlap becomes
        # too low after many skipped frames, it is reported via
        # warn_low_overlap.
        tracker_A = StatePhaseTracker(
            name="A",
            vote_history=config.phase_vote_history,
            vote_decay=config.phase_vote_decay,
            vote_ambiguity_ratio=config.phase_vote_ambiguity_ratio,
            reference_history=config.phase_reference_history,
            sigma_accept_threshold=config.phase_sigma_accept_threshold,
            sigma_filtering_enabled=config.phase_sigma_filtering_enabled,
            corresponding_orbital_alignment=config.phase_corresponding_orbital_alignment,
            invalidate_low_primary_sigma=config.phase_invalidate_low_primary_sigma,
        )
        tracker_B = StatePhaseTracker(
            name="B",
            vote_history=config.phase_vote_history,
            vote_decay=config.phase_vote_decay,
            vote_ambiguity_ratio=config.phase_vote_ambiguity_ratio,
            reference_history=config.phase_reference_history,
            sigma_accept_threshold=config.phase_sigma_accept_threshold,
            sigma_filtering_enabled=config.phase_sigma_filtering_enabled,
            corresponding_orbital_alignment=config.phase_corresponding_orbital_alignment,
            invalidate_low_primary_sigma=config.phase_invalidate_low_primary_sigma,
        )
        # Previous frame's gauge-corrected H_AB for post-correction
        # sign-flip detection (NaN until the first valid frame).
        prev_H_AB_corr = float("nan")
        prev_S_AB_corr = float("nan")
        # QM coordinates at the last *successfully processed* frame (Angstrom).
        # Used only when cross_overlap_mode == "odin".
        phase_qm_coords_history: List[np.ndarray] = []
        phase_adaptive_previous_time_fs: Optional[float] = None
        pending_phase_refinement_frames: List[
            Tuple[int, Optional[float], np.ndarray, np.ndarray, bool]
        ] = []
        pending_phase_refinement_charge_snapshot: Optional[
            Dict[str, Optional[bytes]]
        ] = None

        def write_phase_lookback_rows(frame_id, time_val, state_name, tracker):
            if phase_lookback_file is None:
                return
            if not (
                tracker.last_lookback_used
                or tracker.last_ambiguous
                or tracker.last_low_sigma_refs > 0
            ):
                return
            threshold = config.phase_sigma_accept_threshold
            for ref_idx, vote in enumerate(tracker.last_reference_votes):
                selected = (
                    1
                    if ref_idx == tracker.last_diagnostic_selected_ref_index
                    else 0
                )
                accepted = int(
                    np.isfinite(vote.sigma_min)
                    and vote.sigma_min >= threshold
                )
                phase_lookback_file.write(
                    f"{frame_id:5d}  {time_val:8.2f}  {state_name:>5s}  "
                    f"{ref_idx:7d}  {selected:8d}  {accepted:8d}  "
                    f"{vote.sign:+5d}  {vote.det:16.10f}  {vote.eff:16.10f}  "
                    f"{vote.weight:16.10f}  {vote.sigma_min:16.10f}  "
                    f"{vote.alpha.det:16.10f}  {vote.beta.det:16.10f}  "
                    f"{vote.alpha.sign:+10d}  {vote.beta.sign:+9d}  "
                    f"{vote.alpha.sigma_min:16.10f}  {vote.beta.sigma_min:16.10f}  "
                    f"{vote.alpha.sigma_max:16.10f}  {vote.beta.sigma_max:16.10f}  "
                    f"{vote.alpha.condition_number:16.10f}  "
                    f"{vote.beta.condition_number:16.10f}  "
                    f"{vote.logabs:16.10f}\n"
                )
            phase_lookback_file.flush()

        def snapshot_fragment_charges() -> Dict[str, Optional[bytes]]:
            snapshot: Dict[str, Optional[bytes]] = {}
            for frag in config.fragments:
                charges_path = work_dir / frag.name / "charges.dat"
                snapshot[frag.name] = (
                    charges_path.read_bytes()
                    if charges_path.exists()
                    else None
                )
            return snapshot

        def restore_fragment_charges(snapshot: Dict[str, Optional[bytes]]) -> None:
            for frag in config.fragments:
                charges_path = work_dir / frag.name / "charges.dat"
                charges = snapshot.get(frag.name)
                if charges is None:
                    if charges_path.exists():
                        charges_path.unlink()
                    continue
                charges_path.parent.mkdir(parents=True, exist_ok=True)
                charges_path.write_bytes(charges)

        def truncate_last_output_row(handle) -> None:
            """Remove the last non-header row from an already-flushed text file."""
            handle.flush()
            handle.seek(0)
            lines = handle.readlines()
            header_count = 0
            for line in lines:
                if line.startswith("#"):
                    header_count += 1
                else:
                    break
            if len(lines) <= header_count:
                handle.seek(0, 2)
                return
            handle.seek(0)
            handle.truncate()
            handle.writelines(lines[:-1])
            handle.flush()

        def enqueue_phase_refinement_frames(frame_indices: List[int]) -> None:
            nonlocal pending_phase_refinement_charge_snapshot
            if not frame_indices:
                return
            pending_phase_refinement_charge_snapshot = current_charge_snapshot
            wanted = set(frame_indices)
            max_wanted = max(wanted)
            for fine_frame_id, fine_time_fs, fine_qm_bohr, fine_mm_ang in iter_qm_coordinates(
                config.phase_adaptive_refinement_traj_path,
                config.topology_path,
                qm_indices,
                mm_indices,
            ):
                if fine_frame_id in wanted:
                    pending_phase_refinement_frames.append(
                        (fine_frame_id, fine_time_fs, fine_qm_bohr, fine_mm_ang, True)
                    )
                if fine_frame_id >= max_wanted:
                    break

        coarse_frame_iter = iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices, mm_indices
        )

        while True:
            if pending_phase_refinement_frames:
                frame_id, time_fs, qm_coords_bohr, mm_coords_ang, frame_is_fine = (
                    pending_phase_refinement_frames.pop(0)
                )
            else:
                while True:
                    try:
                        frame_id, time_fs, qm_coords_bohr, mm_coords_ang = next(coarse_frame_iter)
                    except StopIteration:
                        return
                    frame_is_fine = False
                    # Skip frames before start_frame
                    if frame_id < config.start_frame:
                        continue
                    # Stop if we've processed enough coarse frames
                    if (
                        config.n_frames is not None
                        and coarse_processed_count >= config.n_frames
                    ):
                        return
                    break
            
            # Convert QM coords to Angstrom
            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            
            # Calculate time using user-defined t0 and dt
            frame_dt_fs = (
                config.phase_adaptive_refinement_fine_dt_fs
                if frame_is_fine
                else config.dt_fs
            )
            time_val = config.t0_fs + frame_id * frame_dt_fs
            time_str = f"t = {time_val:.3f} fs"
            
            source_label = "fine" if frame_is_fine else "coarse"
            print(f"Frame {frame_id:05d} ({time_str}, {source_label})")

            if frame_is_fine and pending_phase_refinement_charge_snapshot is not None:
                restore_fragment_charges(pending_phase_refinement_charge_snapshot)
                pending_phase_refinement_charge_snapshot = None

            current_charge_snapshot = (
                snapshot_fragment_charges()
                if not frame_is_fine
                else {}
            )
            
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
            phase_refinement_triggered = False
            
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
                energy, mcharge, success, retry_info, scc_attempts = run_cdftb_with_retry(
                    work_dir, frag, qm_coords_bohr, config.hsd_template,
                    config.dftb_library_path, config.num_threads, config.timeout,
                    use_initial_charges, config.ci_enabled,
                    retry_attempts=config.retry_attempts,
                )
                retry_infos.append(retry_info)

                # Log SCC convergence info for every attempt (including retries)
                if scc_log_file is not None:
                    for entry in scc_attempts:
                        # Backward compat: entry may be 3- or 4-tuple
                        if len(entry) == 4:
                            attempt_label, iSCC, scc_err, status = entry
                        else:
                            attempt_label, iSCC, scc_err = entry
                            status = "?"
                        iSCC_str = f"{iSCC:5d}" if iSCC is not None else "  nan"
                        err_str = f"{scc_err:.12e}" if scc_err is not None else "nan"
                        # Compress whitespace in status so the column stays single-token.
                        status_str = " ".join(str(status).split()) or "?"
                        scc_log_file.write(
                            f"{frame_id:5d}  {time_val:12.3f}  "
                            f"{frag.name:>12s}  {attempt_label:>18s}  "
                            f"{iSCC_str}  {err_str}  {status_str}\n"
                        )
                    scc_log_file.flush()
                
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
            
            # Compute CDFTB-CI if enabled
            if config.ci_enabled and len(config.fragments) == 2:
                if not cdftb_success:
                    # CDFTB failed - write NaN values
                    print(f"  CDFTB-CI: SKIPPED (CDFTB convergence failed)")
                    ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {'nan':>12s}  "
                                  f"{'nan':>16s}  {'nan':>16s}  {'nan':>10s}\n")
                    ci_file.flush()
                    ci_sub_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 17) + "\n")
                    ci_sub_file.flush()
                    # Keep previous-frame MOs in the trackers as the phase
                    # reference; the next successful frame will overwrite
                    # them. Just record a SKIP marker in the phase log.
                    if config.phase_tracking_enabled and phase_file is not None:
                        phase_file.write(
                            f"{frame_id:5d}  {time_val:8.2f}  "
                            + "  ".join(["nan"] * 27) + "  SKIP\n"
                        )
                        phase_file.flush()
                else:
                    E_A = energies[0]
                    E_B = energies[1]
                    
                    ham, J_direct, J_lowdin, ci_error, orb_A_data, orb_B_data, S_AO_data = compute_cdftbci_for_frame(
                        work_dir,
                        config.fragments[0].name,
                        config.fragments[1].name,
                        E_A, E_B,
                        config.ci_N_A, config.ci_N_B,
                        config.ci_M_A, config.ci_M_B
                    )
                    
                    if ci_error:
                        print(f"  CDFTB-CI: ERROR - {ci_error}")
                        # Write NaN values
                        ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {'nan':>12s}  "
                                      f"{'nan':>16s}  {'nan':>16s}  {'nan':>10s}\n")
                        ci_sub_file.write(f"{frame_id:5d}  {time_val:8.2f}  " + "  ".join(["nan"] * 17) + "\n")
                        if config.phase_tracking_enabled and phase_file is not None:
                            phase_file.write(
                                f"{frame_id:5d}  {time_val:8.2f}  "
                                + "  ".join(["nan"] * 27) + "  SKIP\n"
                            )
                            phase_file.flush()
                    else:
                        # ---- Phase (gauge) tracking --------------------------
                        H_AB_raw = ham.H_AB
                        S_AB_raw = ham.S_AB
                        S_AB_alpha_raw = ham.S_AB_alpha
                        S_AB_beta_raw = ham.S_AB_beta
                        W_BA_raw = ham.W_BA
                        W_AB_raw = ham.W_AB
                        W_BA_alpha_raw = ham.W_BA_alpha
                        W_BA_beta_raw = ham.W_BA_beta
                        W_AB_alpha_raw = ham.W_AB_alpha
                        W_AB_beta_raw = ham.W_AB_beta
                        W_M_BA_raw = ham.W_M_BA
                        W_M_AB_raw = ham.W_M_AB

                        s_A = 1
                        s_B = 1
                        D_A_raw = float("nan")
                        D_A_corr = float("nan")
                        D_A_eff = float("nan")
                        D_B_raw = float("nan")
                        D_B_corr = float("nan")
                        D_B_eff = float("nan")
                        continuity_override: Optional[str] = None
                        phase_frame_invalid = False
                        phase_refinement_triggered = False
                        requested_fine_indices: List[int] = []
                        S_ao_cross_history: Optional[List[Optional[np.ndarray]]] = None

                        if config.phase_tracking_enabled and orb_A_data is not None and orb_B_data is not None:
                            # Compute exact cross-geometry AO overlaps via ODIN
                            # when requested. One matrix is generated per
                            # stored history reference.
                            if (
                                config.phase_tracking_cross_overlap_mode == "odin"
                                and config.phase_odin_executable is not None
                                and phase_qm_coords_history
                            ):
                                n_history_refs = min(
                                    max(len(tracker_A.ref_history), len(tracker_B.ref_history)),
                                    config.phase_reference_history,
                                    len(phase_qm_coords_history),
                                )
                                S_ao_cross_history = []
                                for hist_idx, hist_coords_ang in enumerate(
                                    phase_qm_coords_history[:n_history_refs]
                                ):
                                    try:
                                        odin_subdir = None
                                        if config.phase_odin_work_subdir:
                                            odin_subdir = config.output_dir / config.phase_odin_work_subdir
                                            if config.phase_odin_keep_files:
                                                odin_subdir = odin_subdir / (
                                                    f"frame_{frame_id:05d}_ref_{hist_idx:02d}"
                                                )
                                        S_ao_cross_history.append(
                                            compute_cross_overlap_odin(
                                                hist_coords_ang,
                                                qm_coords_ang,
                                                atom_types,
                                                config.phase_odin_lmax,
                                                config.phase_odin_sk_prefix,
                                                config.phase_odin_sk_separator,
                                                config.phase_odin_sk_suffix,
                                                config.phase_odin_executable,
                                                remove_translation_rotation=(
                                                    config.phase_remove_translation_rotation_before_overlap
                                                ),
                                                work_dir=odin_subdir,
                                                keep_files=config.phase_odin_keep_files,
                                            )
                                        )
                                    except Exception as _odin_exc:
                                        print(
                                            f"  [ODIN] Warning: cross-overlap computation failed "
                                            f"for history ref {hist_idx}: {_odin_exc}. "
                                            "Falling back to current-frame S_AO for that reference."
                                        )
                                        S_ao_cross_history.append(None)

                        if config.phase_tracking_enabled and orb_A_data is not None and orb_B_data is not None:
                            _cross_mode = config.phase_tracking_cross_overlap_mode
                            _S_prev_for_tracker = (
                                S_ao_cross_history[0]
                                if (
                                    _cross_mode == "odin"
                                    and S_ao_cross_history
                                    and S_ao_cross_history[0] is not None
                                )
                                else S_AO_data
                            )
                            has_phase_reference = (
                                tracker_A.prev_occ_alpha is not None
                                and tracker_B.prev_occ_alpha is not None
                            )
                            if has_phase_reference:
                                s_A, D_A_raw = tracker_A.update(
                                    orb_A_data.C_alpha, orb_A_data.C_beta,
                                    orb_A_data.n_alpha, orb_A_data.n_beta,
                                    S_AO_data,
                                    S_ao_previous=_S_prev_for_tracker,
                                    cross_overlap_mode=_cross_mode,
                                    S_ao_cross_history=S_ao_cross_history,
                                    commit=False,
                                )
                                s_B, D_B_raw = tracker_B.update(
                                    orb_B_data.C_alpha, orb_B_data.C_beta,
                                    orb_B_data.n_alpha, orb_B_data.n_beta,
                                    S_AO_data,
                                    S_ao_previous=_S_prev_for_tracker,
                                    cross_overlap_mode=_cross_mode,
                                    S_ao_cross_history=S_ao_cross_history,
                                    commit=False,
                                )
                                phase_frame_invalid = (
                                    tracker_A.last_invalid
                                    or tracker_B.last_invalid
                                )
                                if not frame_is_fine:
                                    sigma_min_A_probe = min(
                                        tracker_A.last_sigma_min_alpha,
                                        tracker_A.last_sigma_min_beta,
                                    )
                                    sigma_min_B_probe = min(
                                        tracker_B.last_sigma_min_alpha,
                                        tracker_B.last_sigma_min_beta,
                                    )
                                    requested_fine_indices = _phase_adaptive_refinement_indices(
                                        enabled=config.phase_adaptive_refinement_enabled,
                                        previous_time_fs=phase_adaptive_previous_time_fs,
                                        current_time_fs=time_val,
                                        sigma_min_a=sigma_min_A_probe,
                                        sigma_min_b=sigma_min_B_probe,
                                        threshold=config.phase_adaptive_refinement_sigma_min_threshold,
                                        fine_dt_fs=config.phase_adaptive_refinement_fine_dt_fs,
                                    )
                                    phase_refinement_triggered = bool(requested_fine_indices)
                                if not phase_frame_invalid and not phase_refinement_triggered:
                                    s_A, D_A_raw = tracker_A.update(
                                        orb_A_data.C_alpha, orb_A_data.C_beta,
                                        orb_A_data.n_alpha, orb_A_data.n_beta,
                                        S_AO_data,
                                        S_ao_previous=_S_prev_for_tracker,
                                        cross_overlap_mode=_cross_mode,
                                        S_ao_cross_history=S_ao_cross_history,
                                    )
                                    s_B, D_B_raw = tracker_B.update(
                                        orb_B_data.C_alpha, orb_B_data.C_beta,
                                        orb_B_data.n_alpha, orb_B_data.n_beta,
                                        S_AO_data,
                                        S_ao_previous=_S_prev_for_tracker,
                                        cross_overlap_mode=_cross_mode,
                                        S_ao_cross_history=S_ao_cross_history,
                                    )
                                else:
                                    s_A = tracker_A.s
                                    s_B = tracker_B.s
                            else:
                                s_A, D_A_raw = tracker_A.update(
                                    orb_A_data.C_alpha, orb_A_data.C_beta,
                                    orb_A_data.n_alpha, orb_A_data.n_beta,
                                    S_AO_data,
                                    S_ao_previous=_S_prev_for_tracker,
                                    cross_overlap_mode=_cross_mode,
                                    S_ao_cross_history=S_ao_cross_history,
                                )
                                s_B, D_B_raw = tracker_B.update(
                                    orb_B_data.C_alpha, orb_B_data.C_beta,
                                    orb_B_data.n_alpha, orb_B_data.n_beta,
                                    S_AO_data,
                                    S_ao_previous=_S_prev_for_tracker,
                                    cross_overlap_mode=_cross_mode,
                                    S_ao_cross_history=S_ao_cross_history,
                                )
                            D_A_corr = tracker_A.last_D_corr
                            D_B_corr = tracker_B.last_D_corr
                            D_A_eff = effective_state_overlap_metric(
                                D_A_raw,
                                orb_A_data.n_alpha,
                                orb_A_data.n_beta,
                                logabs_D=tracker_A.last_logabs_D,
                            )
                            D_B_eff = effective_state_overlap_metric(
                                D_B_raw,
                                orb_B_data.n_alpha,
                                orb_B_data.n_beta,
                                logabs_D=tracker_B.last_logabs_D,
                            )
                            if phase_frame_invalid:
                                D_A_corr = float("nan")
                                D_B_corr = float("nan")

                            gauge = s_A * s_B
                            if (
                                config.phase_continuity_override_enabled
                                and not phase_frame_invalid
                            ):
                                continuity_override = choose_phase_continuity_override(
                                    gauge=gauge,
                                    H_raw=H_AB_raw,
                                    S_raw=S_AB_raw,
                                    prev_H_corr=prev_H_AB_corr,
                                    prev_S_corr=prev_S_AB_corr,
                                    D_A_eff=D_A_eff,
                                    D_B_eff=D_B_eff,
                                    overlap_threshold=config.phase_continuity_overlap_threshold,
                                )
                                if continuity_override == "A":
                                    tracker_A.flip_current_sign()
                                    s_A = tracker_A.s
                                    D_A_corr = tracker_A.last_D_corr
                                elif continuity_override == "B":
                                    tracker_B.flip_current_sign()
                                    s_B = tracker_B.s
                                    D_B_corr = tracker_B.last_D_corr
                                gauge = s_A * s_B

                            # Apply gauge correction: H_AB and all <Phi_B|...|Phi_A>
                            # quantities transform with s_A * s_B.
                            ham.H_AB = gauge * H_AB_raw
                            ham.S_AB = gauge * S_AB_raw
                            ham.S_AB_alpha = gauge * S_AB_alpha_raw
                            ham.S_AB_beta = gauge * S_AB_beta_raw
                            ham.W_BA = gauge * W_BA_raw
                            ham.W_AB = gauge * W_AB_raw
                            ham.W_BA_alpha = gauge * W_BA_alpha_raw
                            ham.W_BA_beta = gauge * W_BA_beta_raw
                            ham.W_AB_alpha = gauge * W_AB_alpha_raw
                            ham.W_AB_beta = gauge * W_AB_beta_raw
                            ham.W_M_BA = gauge * W_M_BA_raw
                            ham.W_M_AB = gauge * W_M_AB_raw
                            ham.H = np.array([[ham.E_A, ham.H_AB],
                                              [ham.H_AB, ham.E_B]])
                            ham.S = np.array([[1.0, ham.S_AB],
                                              [ham.S_AB, 1.0]])
                            # Transfer integrals must be recomputed from the
                            # gauge-corrected matrices.
                            J_direct = compute_transfer_integral_unrestricted(
                                ham.H, ham.S, method="direct")
                            J_lowdin = compute_transfer_integral_unrestricted(
                                ham.H, ham.S, method="lowdin")

                        if (
                            config.overlap_matrices_binary_enabled
                            and not phase_refinement_triggered
                        ):
                            save_frame_overlap_matrices_binary(
                                config.overlap_matrices_binary_dir,
                                frame_id,
                                time_val,
                                tracker_A_votes=(
                                    tracker_A.last_reference_votes
                                    if config.phase_tracking_enabled
                                    else None
                                ),
                                tracker_B_votes=(
                                    tracker_B.last_reference_votes
                                    if config.phase_tracking_enabled
                                    else None
                                ),
                                selected_ref_index_A=(
                                    tracker_A.last_selected_ref_index
                                    if config.phase_tracking_enabled
                                    else -1
                                ),
                                selected_ref_index_B=(
                                    tracker_B.last_selected_ref_index
                                    if config.phase_tracking_enabled
                                    else -1
                                ),
                            )

                        gauge = s_A * s_B
                        # ------------------------------------------------------

                        # Solve eigenvalue problem
                        eigenvalues, _ = solve_cdftbci_unrestricted(ham.H, ham.S)
                        
                        J_direct_meV = J_direct * 27211.386
                        J_lowdin_meV = J_lowdin * 27211.386
                        J_direct_meV_out = (
                            float("nan") if phase_frame_invalid else J_direct_meV
                        )
                        J_lowdin_meV_out = (
                            float("nan") if phase_frame_invalid else J_lowdin_meV
                        )
                        dE_eV = (eigenvalues[1] - eigenvalues[0]) * 27.211386
                        
                        if config.phase_tracking_enabled:
                            invalid_note = " INVALID" if phase_frame_invalid else ""
                            print(f"  CDFTB-CI: S_AB={ham.S_AB:+.6f}, J={J_lowdin_meV_out:+.2f} meV, ΔE={dE_eV:.4f} eV "
                                  f"(s_A={s_A:+d}, s_B={s_B:+d}){invalid_note}")
                            if phase_refinement_triggered:
                                print(
                                    "  [PHASE REFINE] coarse probe is not written "
                                    "to main output; rerunning previous interval "
                                    f"with fine frames {requested_fine_indices}"
                                )
                        else:
                            print(f"  CDFTB-CI: S_AB={ham.S_AB:.6f}, J={J_lowdin_meV:.2f} meV, ΔE={dE_eV:.4f} eV")
                        
                        # Write CI results (gauge-corrected if enabled)
                        if not phase_refinement_triggered:
                            ci_file.write(f"{frame_id:5d}  {time_val:8.2f}  {J_lowdin_meV_out:12.4f}  "
                                          f"{eigenvalues[0]:16.10f}  {eigenvalues[1]:16.10f}  {dE_eV:10.6f}\n")
                            ci_file.flush()
                        
                        # Write CI sub values
                        if not phase_refinement_triggered:
                            ci_sub_file.write(
                                f"{frame_id:5d}  {time_val:8.2f}  {E_A:14.10f}  {E_B:14.10f}  "
                                f"{ham.H_AB:14.10f}  {J_direct_meV_out:12.4f}  "
                                f"{ham.V_A:14.10f}  {ham.V_B:14.10f}  "
                                f"{ham.N_A:6.1f}  {ham.N_B:6.1f}  {ham.S_AB:14.10f}  "
                                f"{ham.S_AB_alpha:14.10f}  {ham.S_AB_beta:14.10f}  "
                                f"{ham.W_BA:14.10f}  {ham.W_BA_alpha:14.10f}  {ham.W_BA_beta:14.10f}  "
                                f"{ham.W_AB:14.10f}  {ham.W_AB_alpha:14.10f}  {ham.W_AB_beta:14.10f}\n"
                            )
                            ci_sub_file.flush()

                        if config.phase_tracking_enabled:
                            # ------- Phase-tracking sanity diagnostics -------
                            flags: list = []
                            warn_flags: list = []
                            sigma_min_A = min(
                                tracker_A.last_sigma_min_alpha,
                                tracker_A.last_sigma_min_beta,
                            )
                            sigma_min_B = min(
                                tracker_B.last_sigma_min_alpha,
                                tracker_B.last_sigma_min_beta,
                            )
                            if not frame_is_fine:
                                if phase_adaptive_file is not None:
                                    fine_frames_text = (
                                        ",".join(str(i) for i in requested_fine_indices)
                                        if requested_fine_indices
                                        else "-"
                                    )
                                    phase_adaptive_file.write(
                                        f"{frame_id:5d}  {time_val:8.2f}  "
                                        f"{sigma_min_A:16.10f}  {sigma_min_B:16.10f}  "
                                        f"{int(bool(requested_fine_indices)):7d}  "
                                        f"{fine_frames_text}\n"
                                    )
                                    phase_adaptive_file.flush()
                                enqueue_phase_refinement_frames(requested_fine_indices)
                                if phase_refinement_triggered:
                                    truncate_last_output_row(energy_file)
                                    truncate_last_output_row(charge_file)
                                    truncate_last_output_row(retry_log_file)
                            if tracker_A.last_invalid:
                                flags.append("INVALID_A")
                                warn_flags.append("INVALID_A")
                                warn_counts["INVALID_A"] += 1
                            elif tracker_A.last_lookback_used:
                                flags.append("LOOKBACK_A")
                                warn_flags.append("LOOKBACK_A")
                                warn_counts["LOOKBACK_A"] += 1
                            if tracker_B.last_invalid:
                                flags.append("INVALID_B")
                                warn_flags.append("INVALID_B")
                                warn_counts["INVALID_B"] += 1
                            elif tracker_B.last_lookback_used:
                                flags.append("LOOKBACK_B")
                                warn_flags.append("LOOKBACK_B")
                                warn_counts["LOOKBACK_B"] += 1
                            if tracker_A.last_ambiguous and not tracker_A.last_invalid:
                                flags.append("AMBIG_A")
                                warn_flags.append("AMBIG_A")
                                warn_counts["AMBIG_A"] += 1
                            if tracker_B.last_ambiguous and not tracker_B.last_invalid:
                                flags.append("AMBIG_B")
                                warn_flags.append("AMBIG_B")
                                warn_counts["AMBIG_B"] += 1
                            if phase_refinement_triggered:
                                flags.append("REFINE_PREV")
                            if np.isfinite(sigma_min_A) and sigma_min_A < config.phase_warn_low_sigma_min:
                                flags.append("LOW_SIG_A")
                                warn_flags.append("LOW_SIG_A")
                                warn_counts["LOW_SIG_A"] += 1
                            if np.isfinite(sigma_min_B) and sigma_min_B < config.phase_warn_low_sigma_min:
                                flags.append("LOW_SIG_B")
                                warn_flags.append("LOW_SIG_B")
                                warn_counts["LOW_SIG_B"] += 1
                            if np.isfinite(D_A_eff):
                                if D_A_eff < config.phase_warn_low_overlap:
                                    flags.append("LOW_DA")
                                    warn_flags.append("LOW_DA")
                                    warn_counts["LOW_DA"] += 1
                                if D_A_eff > config.phase_warn_high_overlap:
                                    flags.append("HIGH_DA")
                                    warn_flags.append("HIGH_DA")
                                    warn_counts["HIGH_DA"] += 1
                            if np.isfinite(D_B_eff):
                                if D_B_eff < config.phase_warn_low_overlap:
                                    flags.append("LOW_DB")
                                    warn_flags.append("LOW_DB")
                                    warn_counts["LOW_DB"] += 1
                                if D_B_eff > config.phase_warn_high_overlap:
                                    flags.append("HIGH_DB")
                                    warn_flags.append("HIGH_DB")
                                    warn_counts["HIGH_DB"] += 1
                            if continuity_override is not None:
                                cont_flag = f"CONT_{continuity_override}"
                                flags.append(cont_flag)
                                warn_counts[cont_flag] += 1
                            # Post-correction sign flip of H_AB between
                            # consecutive successful frames.
                            if (config.phase_warn_post_correction_flip
                                    and not phase_frame_invalid
                                    and not phase_refinement_triggered
                                    and np.isfinite(prev_H_AB_corr)
                                    and np.isfinite(ham.H_AB)
                                    and prev_H_AB_corr != 0.0 and ham.H_AB != 0.0
                                    and (prev_H_AB_corr * ham.H_AB) < 0.0):
                                flags.append("FLIP_HAB")
                                warn_flags.append("FLIP_HAB")
                                warn_counts["FLIP_HAB"] += 1
                            flags_str = ",".join(flags) if flags else "OK"
                            if warn_flags:
                                warn_frames.append((frame_id, ",".join(warn_flags)))
                                print(f"  [PHASE WARN] frame {frame_id}: {flags_str} "
                                      f"(D_A_eff={D_A_eff:.3f}, D_B_eff={D_B_eff:.3f}, "
                                      f"sigma_min_A={sigma_min_A:.3f}, sigma_min_B={sigma_min_B:.3f}, "
                                      f"ref_A={tracker_A.last_selected_ref_index}, "
                                      f"ref_B={tracker_B.last_selected_ref_index}, "
                                      f"H_AB_corr: {prev_H_AB_corr:+.4f} -> {ham.H_AB:+.4f} Ha)")
                            elif continuity_override is not None:
                                print(f"  [PHASE FIX] frame {frame_id}: {flags_str} "
                                      f"(H_AB_corr: {prev_H_AB_corr:+.4f} -> {ham.H_AB:+.4f} Ha, "
                                      f"S_AB_corr: {prev_S_AB_corr:+.4e} -> {ham.S_AB:+.4e})")
                            # Invalid frames are not phase references. The next
                            # valid frame is compared against the previous
                            # valid frame.
                            if not phase_frame_invalid and not phase_refinement_triggered:
                                prev_H_AB_corr = float(ham.H_AB)
                                prev_S_AB_corr = float(ham.S_AB)

                            if phase_file is not None and not phase_refinement_triggered:
                                phase_file.write(
                                    f"{frame_id:5d}  {time_val:8.2f}  "
                                    f"{s_A:+3d}  {s_B:+3d}  {gauge:+3d}  "
                                    f"{D_A_raw:16.10f}  {D_A_corr:16.10f}  {D_A_eff:16.10f}  "
                                    f"{tracker_A.last_sigma_min_alpha:16.10f}  {tracker_A.last_sigma_min_beta:16.10f}  "
                                    f"{tracker_A.last_cond_alpha:16.10f}  {tracker_A.last_cond_beta:16.10f}  "
                                    f"{tracker_A.last_vote_margin:16.10f}  {tracker_A.last_vote_refs:6d}  "
                                    f"{tracker_A.last_selected_ref_index:9d}  "
                                    f"{D_B_raw:16.10f}  {D_B_corr:16.10f}  {D_B_eff:16.10f}  "
                                    f"{tracker_B.last_sigma_min_alpha:16.10f}  {tracker_B.last_sigma_min_beta:16.10f}  "
                                    f"{tracker_B.last_cond_alpha:16.10f}  {tracker_B.last_cond_beta:16.10f}  "
                                    f"{tracker_B.last_vote_margin:16.10f}  {tracker_B.last_vote_refs:6d}  "
                                    f"{tracker_B.last_selected_ref_index:9d}  "
                                    f"{H_AB_raw:16.10f}  {ham.H_AB:16.10f}  "
                                    f"{S_AB_raw:16.10f}  {ham.S_AB:16.10f}  "
                                    f"{flags_str}\n"
                                )
                                phase_file.flush()

                            write_phase_lookback_rows(frame_id, time_val, "A", tracker_A)
                            write_phase_lookback_rows(frame_id, time_val, "B", tracker_B)

                            # Keep recent successful QM geometries aligned with
                            # the phase reference history. Failed/SKIP frames do
                            # not advance ODIN references.
                            if not phase_frame_invalid and not phase_refinement_triggered:
                                history_limit = max(1, int(config.phase_reference_history))
                                phase_qm_coords_history = (
                                    [qm_coords_ang.copy()]
                                    + phase_qm_coords_history[:history_limit - 1]
                                )
                                phase_adaptive_previous_time_fs = time_val
            
            # Compute spin populations for both constraint states
            if (
                cdftb_success
                and len(config.fragments) == 2
                and not phase_refinement_triggered
            ):
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
            if not frame_is_fine:
                coarse_processed_count += 1
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
        if phase_file:
            phase_file.close()
        if phase_lookback_file:
            phase_lookback_file.close()
        if phase_adaptive_file:
            phase_adaptive_file.close()
        if scc_log_file:
            scc_log_file.close()
        
        print("=" * 70)
        print(f"Energies saved to {config.energy_file}")
        print(f"Charges saved to {config.charge_file}")
        print(f"Spin populations saved to {config.spin_output_file}")
        print(f"Retry log saved to {config.retry_log_file}")
        if config.scc_log_enabled:
            print(f"SCC convergence log saved to {config.scc_log_file}")
        if config.ci_enabled:
            print(f"CDFTB-CI results saved to {config.ci_output_file}")
            print(f"CDFTB-CI sub values saved to {config.ci_sub_file}")
            if config.overlap_matrices_binary_enabled:
                print(f"CDFTB-CI overlap matrices saved to {config.overlap_matrices_binary_dir}")
            if config.phase_tracking_enabled:
                print(f"CDFTB-CI phase tracking saved to {config.phase_output_file}")
                if phase_lookback_file is not None:
                    print(f"CDFTB-CI phase lookback saved to {phase_lookback_file.name}")
                warning_keys = (
                    "LOW_DA", "LOW_DB",
                    "LOW_SIG_A", "LOW_SIG_B",
                    "LOOKBACK_A", "LOOKBACK_B",
                    "INVALID_A", "INVALID_B",
                    "HIGH_DA", "HIGH_DB",
                    "AMBIG_A", "AMBIG_B",
                    "FLIP_HAB",
                )
                total_warn = sum(warn_counts[k] for k in warning_keys)
                if total_warn == 0:
                    print("Phase tracking: no warnings (all frames passed sanity checks).")
                else:
                    print(f"Phase tracking: {total_warn} warning(s) raised:")
                    for k in warning_keys:
                        v = warn_counts[k]
                        if v > 0:
                            print(f"  {k}: {v}")
                    if warn_frames:
                        preview = ", ".join(f"{fid}:{flg}" for fid, flg in warn_frames[:10])
                        more = "" if len(warn_frames) <= 10 else f", ... (+{len(warn_frames)-10} more)"
                        print(f"  Frames flagged: {preview}{more}")
                for k in ("CONT_A", "CONT_B"):
                    if warn_counts[k] > 0:
                        print(f"  {k}: {warn_counts[k]} continuity override(s)")
        print("=" * 70)
