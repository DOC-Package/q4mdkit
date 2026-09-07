"""CDFT-CI core mathematics (dependency-light, NumPy only).

These routines implement the Wu--Van Voorhis constrained-DFT configuration
interaction (CDFT-CI) coupling for spin-polarized (unrestricted) states.  They
are shared by the DFTB+ (``cdftbci.py``) and NWChem (``nwchem_cdftci.py``)
backends, so only the data-extraction layer differs between backends.

References:
    Wu and Van Voorhis, J. Chem. Phys. 125, 164105 (2006)
    Wu et al., J. Chem. Phys. 127, 164119 (2007)
    Oberhofer and Blumberger, J. Chem. Phys. 133, 244105 (2010)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np


@dataclass
class UnrestrictedOrbitalData:
    """Container for unrestricted orbital data (alpha and beta separately)."""

    C_alpha: np.ndarray  # alpha MO coefficients (n_AO x n_MO)
    C_beta: np.ndarray   # beta MO coefficients (n_AO x n_MO)
    occ_alpha: np.ndarray  # alpha occupation numbers
    occ_beta: np.ndarray   # beta occupation numbers
    n_alpha: int  # number of alpha electrons
    n_beta: int   # number of beta electrons


@dataclass
class UnrestrictedCDFTCIHamiltonian:
    """Container for the 2x2 CDFT-CI Hamiltonian and related quantities."""

    H: np.ndarray  # 2x2 Hamiltonian matrix
    S: np.ndarray  # 2x2 overlap matrix
    E_A: float
    E_B: float
    H_AB: float
    S_AB: float
    S_AB_alpha: float
    S_AB_beta: float
    W_BA: float
    W_AB: float
    W_BA_alpha: float
    W_BA_beta: float
    W_AB_alpha: float
    W_AB_beta: float
    V_A: float = 0.0
    V_B: float = 0.0
    N_A: float = 0.0
    N_B: float = 0.0
    V_M_A: float = 0.0
    V_M_B: float = 0.0
    M_A: float = 0.0
    M_B: float = 0.0
    W_M_BA: float = 0.0
    W_M_AB: float = 0.0

    def save(
        self, filepath: str | Path, frame_id: int = 0, time_fs: float = 0.0
    ) -> None:
        """Append the Hamiltonian parameters to ``filepath``."""
        filepath = Path(filepath)
        write_header = not filepath.exists()
        has_spin = (
            self.V_M_A != 0.0
            or self.V_M_B != 0.0
            or self.M_A != 0.0
            or self.M_B != 0.0
        )

        with filepath.open("a", encoding="utf-8") as handle:
            if write_header:
                handle.write("# CDFT-CI Hamiltonian Parameters\n")
                if has_spin:
                    handle.write(
                        "# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_N_A(Ha)  V_N_B(Ha)  "
                        "N_A  N_B  V_M_A(Ha)  V_M_B(Ha)  M_A  M_B  "
                        "S_AB  S_AB_alpha  S_AB_beta  W_BA  W_AB  W_M_BA  W_M_AB  "
                        "H_AB(Ha)  J(meV)\n"
                    )
                else:
                    handle.write(
                        "# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  V_A(Ha)  V_B(Ha)  "
                        "N_A  N_B  S_AB  S_AB_alpha  S_AB_beta  W_BA  W_AB  "
                        "H_AB(Ha)  J(meV)\n"
                    )

            J = self.H_AB - self.S_AB * (self.E_A + self.E_B) / 2
            J_meV = J * 27211.386

            if has_spin:
                handle.write(
                    f"{frame_id:6d}  {time_fs:8.2f}  {self.E_A:14.10f}  {self.E_B:14.10f}  "
                    f"{self.V_A:12.8f}  {self.V_B:12.8f}  {self.N_A:6.1f}  {self.N_B:6.1f}  "
                    f"{self.V_M_A:12.8f}  {self.V_M_B:12.8f}  {self.M_A:6.1f}  {self.M_B:6.1f}  "
                    f"{self.S_AB:12.8f}  {self.S_AB_alpha:12.8f}  {self.S_AB_beta:12.8f}  "
                    f"{self.W_BA:12.8f}  {self.W_AB:12.8f}  {self.W_M_BA:12.8f}  {self.W_M_AB:12.8f}  "
                    f"{self.H_AB:14.10f}  {J_meV:10.4f}\n"
                )
            else:
                handle.write(
                    f"{frame_id:6d}  {time_fs:8.2f}  {self.E_A:14.10f}  {self.E_B:14.10f}  "
                    f"{self.V_A:12.8f}  {self.V_B:12.8f}  {self.N_A:6.1f}  {self.N_B:6.1f}  "
                    f"{self.S_AB:12.8f}  {self.S_AB_alpha:12.8f}  {self.S_AB_beta:12.8f}  "
                    f"{self.W_BA:12.8f}  {self.W_AB:12.8f}  {self.H_AB:14.10f}  {J_meV:10.4f}\n"
                )


def compute_mo_overlap_matrix_spin(
    C_A: np.ndarray,
    C_B: np.ndarray,
    S_AO: np.ndarray,
    n_occ: int,
) -> np.ndarray:
    """Compute the occupied MO overlap O^BA = C_B^T S_AO C_A."""
    if n_occ == 0:
        return np.array([[1.0]])
    C_A_occ = C_A[:, :n_occ]
    C_B_occ = C_B[:, :n_occ]
    return C_B_occ.T @ S_AO @ C_A_occ


def compute_omega_matrix_spin(
    C_A: np.ndarray,
    C_B: np.ndarray,
    w: np.ndarray,
    n_occ: int,
) -> np.ndarray:
    """Compute the occupied constraint-weight overlap Omega^BA = C_B^T w C_A."""
    if n_occ == 0:
        return np.array([[0.0]])
    C_A_occ = C_A[:, :n_occ]
    C_B_occ = C_B[:, :n_occ]
    return C_B_occ.T @ w @ C_A_occ


def compute_state_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray,
) -> Tuple[float, float, float]:
    """Return S_AB = det(O_alpha) * det(O_beta) and both determinants."""
    S_AB_alpha = float(np.linalg.det(O_BA_alpha))
    S_AB_beta = float(np.linalg.det(O_BA_beta))
    return S_AB_alpha * S_AB_beta, S_AB_alpha, S_AB_beta


def compute_weight_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray,
    Omega_BA_alpha: np.ndarray,
    Omega_BA_beta: np.ndarray,
) -> Tuple[float, float, float]:
    """Return W_BA = S_AB * (Tr(O_alpha^-1 Omega_alpha) + Tr(O_beta^-1 Omega_beta))."""
    det_alpha = float(np.linalg.det(O_BA_alpha))
    det_beta = float(np.linalg.det(O_BA_beta))
    S_AB = det_alpha * det_beta

    trace_alpha = (
        float(np.trace(np.linalg.inv(O_BA_alpha) @ Omega_BA_alpha))
        if O_BA_alpha.shape[0] > 0
        else 0.0
    )
    trace_beta = (
        float(np.trace(np.linalg.inv(O_BA_beta) @ Omega_BA_beta))
        if O_BA_beta.shape[0] > 0
        else 0.0
    )

    W_BA_alpha = S_AB * trace_alpha
    W_BA_beta = S_AB * trace_beta
    return W_BA_alpha + W_BA_beta, W_BA_alpha, W_BA_beta


def compute_spin_weight_overlap_unrestricted(
    O_BA_alpha: np.ndarray,
    O_BA_beta: np.ndarray,
    Omega_BA_alpha: np.ndarray,
    Omega_BA_beta: np.ndarray,
) -> float:
    """Return spin weight overlap W_M = S_AB * (trace_alpha - trace_beta)."""
    S_AB, _, _ = compute_state_overlap_unrestricted(O_BA_alpha, O_BA_beta)

    trace_alpha = (
        float(np.trace(np.linalg.inv(O_BA_alpha) @ Omega_BA_alpha))
        if O_BA_alpha.shape[0] > 0
        else 0.0
    )
    trace_beta = (
        float(np.trace(np.linalg.inv(O_BA_beta) @ Omega_BA_beta))
        if O_BA_beta.shape[0] > 0
        else 0.0
    )
    return S_AB * (trace_alpha - trace_beta)


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
    W_M_AB: float = 0.0,
) -> float:
    """Compute the CDFT-CI coupling element H_AB (charge and/or spin)."""
    term1 = 0.5 * (E_A + E_B + N_A * V_A + N_B * V_B) * S_AB
    term2 = 0.5 * (V_A * W_BA + V_B * W_AB)
    term1 += 0.5 * (M_A * V_M_A + M_B * V_M_B) * S_AB
    term2 += 0.5 * (V_M_A * W_M_BA + V_M_B * W_M_AB)
    return term1 - term2


def build_cdft_ci_hamiltonian_unrestricted(
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
    M_B: float = 0.0,
) -> UnrestrictedCDFTCIHamiltonian:
    """Build the 2x2 CDFT-CI Hamiltonian for unrestricted wavefunctions."""
    O_BA_alpha = compute_mo_overlap_matrix_spin(
        orb_A.C_alpha, orb_B.C_alpha, S_AO, orb_A.n_alpha
    )
    O_BA_beta = compute_mo_overlap_matrix_spin(
        orb_A.C_beta, orb_B.C_beta, S_AO, orb_A.n_beta
    )

    O_AB_alpha = O_BA_alpha.T
    O_AB_beta = O_BA_beta.T

    S_AB, S_AB_alpha, S_AB_beta = compute_state_overlap_unrestricted(
        O_BA_alpha, O_BA_beta
    )

    Omega_BA_alpha = compute_omega_matrix_spin(
        orb_A.C_alpha, orb_B.C_alpha, w_A, orb_A.n_alpha
    )
    Omega_BA_beta = compute_omega_matrix_spin(
        orb_A.C_beta, orb_B.C_beta, w_A, orb_A.n_beta
    )
    Omega_AB_alpha = compute_omega_matrix_spin(
        orb_B.C_alpha, orb_A.C_alpha, w_B, orb_B.n_alpha
    )
    Omega_AB_beta = compute_omega_matrix_spin(
        orb_B.C_beta, orb_A.C_beta, w_B, orb_B.n_beta
    )

    W_BA, W_BA_alpha, W_BA_beta = compute_weight_overlap_unrestricted(
        O_BA_alpha, O_BA_beta, Omega_BA_alpha, Omega_BA_beta
    )
    W_AB, W_AB_alpha, W_AB_beta = compute_weight_overlap_unrestricted(
        O_AB_alpha, O_AB_beta, Omega_AB_alpha, Omega_AB_beta
    )

    W_M_BA = 0.0
    W_M_AB = 0.0
    if V_M_A != 0.0 or V_M_B != 0.0 or M_A != 0.0 or M_B != 0.0:
        W_M_BA = compute_spin_weight_overlap_unrestricted(
            O_BA_alpha, O_BA_beta, Omega_BA_alpha, Omega_BA_beta
        )
        W_M_AB = compute_spin_weight_overlap_unrestricted(
            O_AB_alpha, O_AB_beta, Omega_AB_alpha, Omega_AB_beta
        )

    H_AB = compute_coupling_element_unrestricted(
        E_A, E_B, V_A, V_B, N_A, N_B, S_AB, W_BA, W_AB,
        V_M_A, V_M_B, M_A, M_B, W_M_BA, W_M_AB,
    )

    H = np.array([[E_A, H_AB], [H_AB, E_B]])
    S = np.array([[1.0, S_AB], [S_AB, 1.0]])

    return UnrestrictedCDFTCIHamiltonian(
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
        W_M_AB=W_M_AB,
    )


def solve_cdft_ci_unrestricted(
    H: np.ndarray, S: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve the generalized eigenvalue problem HC = SCE."""
    from scipy.linalg import eigh

    return eigh(H, S)


def compute_transfer_integral_unrestricted(
    H: np.ndarray,
    S: np.ndarray,
    method: str = "direct",
) -> float:
    """Return the effective transfer integral J (direct or Lowdin)."""
    E_A = H[0, 0]
    E_B = H[1, 1]
    H_AB = H[0, 1]
    S_AB = S[0, 1]

    J_direct = H_AB - S_AB * (E_A + E_B) / 2
    if method == "direct":
        return J_direct
    if method == "lowdin":
        return J_direct / (1 - S_AB**2)
    raise ValueError(f"Unknown transfer-integral method: {method}")


def build_mulliken_weight_matrix(
    S_AO: np.ndarray,
    fragment_atoms: list[int],
    n_atoms: int,
    atom_to_orbitals: Optional[dict[int, list[int]]] = None,
) -> np.ndarray:
    """Build a Mulliken fragment weight matrix (w_uv = delta_{u in frag} S_uv)."""
    n_orbitals = S_AO.shape[0]
    w = np.zeros((n_orbitals, n_orbitals), dtype=np.float64)

    if atom_to_orbitals is None:
        orbitals_per_atom = n_orbitals // n_atoms
        atom_to_orbitals = {
            atom_idx: list(
                range((atom_idx - 1) * orbitals_per_atom, atom_idx * orbitals_per_atom)
            )
            for atom_idx in range(1, n_atoms + 1)
        }

    for atom_idx in fragment_atoms:
        for orbital_idx in atom_to_orbitals.get(atom_idx, []):
            w[orbital_idx, :] = S_AO[orbital_idx, :]
    return w


def build_lowdin_weight_matrix(
    S_AO: np.ndarray,
    fragment_atoms: list[int],
    n_atoms: int,
    atom_to_orbitals: Optional[dict[int, list[int]]] = None,
) -> np.ndarray:
    """Build a Lowdin fragment weight matrix ``w = S^{1/2} P_frag S^{1/2}``."""
    n_orbitals = S_AO.shape[0]
    if atom_to_orbitals is None:
        orbitals_per_atom = n_orbitals // n_atoms
        atom_to_orbitals = {
            atom_idx: list(
                range((atom_idx - 1) * orbitals_per_atom, atom_idx * orbitals_per_atom)
            )
            for atom_idx in range(1, n_atoms + 1)
        }

    P = np.zeros((n_orbitals, n_orbitals), dtype=np.float64)
    for atom_idx in fragment_atoms:
        for orbital_idx in atom_to_orbitals.get(atom_idx, []):
            P[orbital_idx, orbital_idx] = 1.0

    eigvals, eigvecs = np.linalg.eigh(S_AO)
    S_half = (eigvecs * np.sqrt(eigvals)) @ eigvecs.T
    return S_half @ P @ S_half
