import numpy as np

from q4mdkit.analysis.cdft_ci_core import (
    UnrestrictedOrbitalData,
    build_cdft_ci_hamiltonian_unrestricted,
    build_lowdin_weight_matrix,
    build_mulliken_weight_matrix,
    compute_mo_overlap_matrix_spin,
    compute_state_overlap_unrestricted,
    compute_transfer_integral_unrestricted,
)


def _two_state_problem():
    S_AO = np.array([[1.0, 0.5], [0.5, 1.0]])

    orb_A = UnrestrictedOrbitalData(
        C_alpha=np.array([[1.0], [0.0]]),
        C_beta=np.array([[1.0], [0.0]]),
        occ_alpha=np.array([1.0]),
        occ_beta=np.array([1.0]),
        n_alpha=1,
        n_beta=1,
    )
    orb_B = UnrestrictedOrbitalData(
        C_alpha=np.array([[0.0], [1.0]]),
        C_beta=np.array([[0.0], [1.0]]),
        occ_alpha=np.array([1.0]),
        occ_beta=np.array([1.0]),
        n_alpha=1,
        n_beta=1,
    )
    atom_to_orbitals = {1: [0], 2: [1]}
    w_A = build_mulliken_weight_matrix(S_AO, [1], 2, atom_to_orbitals)
    w_B = build_mulliken_weight_matrix(S_AO, [2], 2, atom_to_orbitals)
    return S_AO, orb_A, orb_B, w_A, w_B


def test_compute_mo_overlap_matrix_spin():
    S_AO, orb_A, orb_B, _, _ = _two_state_problem()
    O = compute_mo_overlap_matrix_spin(orb_A.C_alpha, orb_B.C_alpha, S_AO, 1)
    assert O.shape == (1, 1)
    assert O[0, 0] == 0.5


def test_compute_state_overlap_unrestricted():
    S_AO, orb_A, orb_B, _, _ = _two_state_problem()
    O_alpha = compute_mo_overlap_matrix_spin(orb_A.C_alpha, orb_B.C_alpha, S_AO, 1)
    O_beta = compute_mo_overlap_matrix_spin(orb_A.C_beta, orb_B.C_beta, S_AO, 1)
    S_AB, S_alpha, S_beta = compute_state_overlap_unrestricted(O_alpha, O_beta)
    assert S_alpha == 0.5
    assert S_beta == 0.5
    assert S_AB == 0.25


def test_build_mulliken_weight_matrix():
    S_AO, _, _, _, _ = _two_state_problem()
    w_A = build_mulliken_weight_matrix(S_AO, [1], 2, {1: [0], 2: [1]})
    w_B = build_mulliken_weight_matrix(S_AO, [2], 2, {1: [0], 2: [1]})
    assert np.allclose(w_A, [[1.0, 0.5], [0.0, 0.0]])
    assert np.allclose(w_B, [[0.0, 0.0], [0.5, 1.0]])


def test_build_lowdin_weight_matrix():
    S_AO = np.array([[1.0, 0.5], [0.5, 1.0]])
    w_A = build_lowdin_weight_matrix(S_AO, [1], 2, {1: [0], 2: [1]})
    w_B = build_lowdin_weight_matrix(S_AO, [2], 2, {1: [0], 2: [1]})

    # w = S^{1/2} P S^{1/2}; check symmetry, PSD, and P_A + P_B = S.
    assert np.allclose(w_A, w_A.T)
    assert np.allclose(w_B, w_B.T)
    assert np.allclose(w_A + w_B, S_AO)


def test_build_hamiltonian_and_transfer_integral():
    S_AO, orb_A, orb_B, w_A, w_B = _two_state_problem()

    ham = build_cdft_ci_hamiltonian_unrestricted(
        orb_A,
        orb_B,
        S_AO,
        w_A,
        w_B,
        E_A=0.5,
        E_B=1.5,
        V_A=1.0,
        V_B=2.0,
        N_A=1.0,
        N_B=1.0,
    )

    assert ham.S_AB == 0.25
    # H_AB = 0.5 * (0.5 + 1.5 + 1.0 + 2.0) * 0.25 = 0.625
    assert np.isclose(ham.H_AB, 0.625)
    assert np.isclose(ham.H[0, 0], 0.5)
    assert np.isclose(ham.H[1, 1], 1.5)

    J = compute_transfer_integral_unrestricted(ham.H, ham.S, method="direct")
    # J = 0.625 - 0.25 * (0.5 + 1.5) / 2 = 0.375
    assert np.isclose(J, 0.375)
