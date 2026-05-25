"""
Phase (gauge) tracking for diabatic states along an MD trajectory.

In CDFTB-CI calculations, each charge-localized diabatic state |Phi_f(t_n)>
on fragment f at time t_n is defined only up to an arbitrary sign s_f(t_n)=±1
(the wavefunction is real-valued in DFTB). The coupling transforms as
    H_AB(t_n) -> s_A(t_n) s_B(t_n) H_AB(t_n),
so when each snapshot is optimized independently the raw H_AB(t) trajectory
can exhibit artificial sign flips that come purely from this gauge freedom.

This module fixes the gauge by overlap-based phase tracking. At each step
the sign of |Phi_f(t_n)> is chosen so that its overlap with the sign-
corrected state at the previous snapshot is positive:

    D_f(t_{n-1}, t_n) = <Phi_f^corr(t_{n-1}) | Phi_f(t_n)>
                      = prod_{sigma=alpha,beta} det O_f^sigma(t_{n-1}, t_n),

    [O_f^sigma]_{ij} = sum_{mu,nu} C^*_{mu i}(t_{n-1}) S^{n-1,n}_{mu nu} C_{nu j}(t_n),

    s_f(t_n) = s_f(t_{n-1}) * sgn[D_f(t_{n-1}, t_n)].

After multiplying H_AB(t_n) and S_AB(t_n) by s_A(t_n) s_B(t_n), the gauge-
corrected trajectory is continuous and is the quantity required when
constructing spectral densities from time-correlation functions of the
transfer-integral fluctuations.

Notes on the cross-geometry AO overlap S^{n-1,n}
------------------------------------------------
The exact quantity is

    S^{n-1,n}_{mu nu} = <chi_mu(t_{n-1}) | chi_nu(t_n)>,

i.e. the AO overlap between basis functions evaluated at the two different
geometries. For DFTB+ this is not directly available without a dedicated
"doubled-system" calculation. For *sign tracking* (the only thing this
module needs) it is sufficient to approximate

    S^{n-1,n}_{mu nu} ~= S_{mu nu}(t_n)         (default)

or, optionally, the midpoint average 1/2 [S(t_{n-1}) + S(t_n)]. For typical
MD time steps (<~ few fs) the displacement is small compared to atomic
orbital extents and this approximation reliably preserves the sign of the
determinant det O_f^sigma even when its magnitude is somewhat reduced.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np


def _occ_block(C: np.ndarray, n_occ: int) -> np.ndarray:
    """Return the occupied block of the MO coefficient matrix."""
    if n_occ <= 0:
        return np.zeros((C.shape[0], 0), dtype=C.dtype)
    return C[:, :n_occ]


def compute_state_overlap_signed(
    C_prev_alpha: np.ndarray, C_curr_alpha: np.ndarray,
    C_prev_beta: np.ndarray, C_curr_beta: np.ndarray,
    n_alpha: int, n_beta: int,
    S_ao_cross: np.ndarray,
) -> Tuple[float, float, float]:
    """
    Compute <Phi(t_{n-1}) | Phi(t_n)> = det O_alpha * det O_beta.

    Parameters
    ----------
    C_prev_alpha, C_prev_beta : np.ndarray
        Raw (uncorrected) MO coefficients at t_{n-1}, shape (n_AO, n_MO).
    C_curr_alpha, C_curr_beta : np.ndarray
        MO coefficients at t_n.
    n_alpha, n_beta : int
        Number of occupied alpha / beta orbitals.
    S_ao_cross : np.ndarray
        Cross-geometry AO overlap matrix S^{n-1,n} (n_AO x n_AO).
        For sign tracking, S(t_n) is a usable approximation.

    Returns
    -------
    D : float
        det(O_alpha) * det(O_beta).
    det_alpha : float
    det_beta : float
    """
    Ca_prev = _occ_block(C_prev_alpha, n_alpha)
    Ca_curr = _occ_block(C_curr_alpha, n_alpha)
    Cb_prev = _occ_block(C_prev_beta, n_beta)
    Cb_curr = _occ_block(C_curr_beta, n_beta)

    if Ca_prev.shape[1] == 0:
        det_alpha = 1.0
    else:
        O_alpha = Ca_prev.T @ S_ao_cross @ Ca_curr
        det_alpha = float(np.linalg.det(O_alpha))

    if Cb_prev.shape[1] == 0:
        det_beta = 1.0
    else:
        O_beta = Cb_prev.T @ S_ao_cross @ Cb_curr
        det_beta = float(np.linalg.det(O_beta))

    return det_alpha * det_beta, det_alpha, det_beta


@dataclass
class StatePhaseTracker:
    """
    Tracks the gauge factor s_f for a single diabatic state along an MD
    trajectory.

    Usage
    -----
        tracker = StatePhaseTracker(name="A")
        for frame in trajectory:
            s, D = tracker.update(C_alpha, C_beta, n_alpha, n_beta, S_ao)
            # use s to gauge-correct couplings involving this state
    """
    name: str = ""
    s: int = 1                               # current gauge factor s_f(t_n)
    prev_C_alpha: Optional[np.ndarray] = field(default=None, repr=False)
    prev_C_beta: Optional[np.ndarray] = field(default=None, repr=False)
    prev_n_alpha: int = 0
    prev_n_beta: int = 0
    last_D: float = 1.0                      # raw <Phi(t_{n-1})|Phi(t_n)>
    last_D_corr: float = 1.0                 # s_{n-1} * raw D (gauge-corrected)
    last_det_alpha: float = 1.0
    last_det_beta: float = 1.0

    def reset(self) -> None:
        self.s = 1
        self.prev_C_alpha = None
        self.prev_C_beta = None
        self.prev_n_alpha = 0
        self.prev_n_beta = 0
        self.last_D = 1.0
        self.last_D_corr = 1.0
        self.last_det_alpha = 1.0
        self.last_det_beta = 1.0

    def update(
        self,
        C_alpha: np.ndarray,
        C_beta: np.ndarray,
        n_alpha: int,
        n_beta: int,
        S_ao_current: np.ndarray,
        S_ao_previous: Optional[np.ndarray] = None,
        cross_overlap_mode: str = "current",
    ) -> Tuple[int, float]:
        """
        Update gauge factor using overlap with the previous frame.

        Parameters
        ----------
        C_alpha, C_beta : np.ndarray
            Raw MO coefficients at the current snapshot t_n.
        n_alpha, n_beta : int
            Number of occupied alpha / beta orbitals at t_n.
        S_ao_current : np.ndarray
            AO overlap matrix at t_n.
        S_ao_previous : np.ndarray, optional
            AO overlap matrix at t_{n-1}. Only used when
            cross_overlap_mode == "midpoint".
        cross_overlap_mode : {"current", "previous", "midpoint", "odin"}
            Approximation for the cross-geometry AO overlap S^{n-1,n}:
              "current"  -> S(t_n)                    (default, robust)
              "previous" -> S(t_{n-1})
              "midpoint" -> 0.5 [S(t_{n-1}) + S(t_n)]
              "odin"     -> use the pre-computed exact cross-overlap passed
                            in ``S_ao_previous`` (computed externally via ODIN).

        Returns
        -------
        s : int
            Updated gauge factor s_f(t_n) in {-1, +1}.
        D : float
            Raw overlap <Phi_f(t_{n-1}) | Phi_f(t_n)> (before applying
            s_f(t_{n-1})). On the first call this is 1.0 by convention.
        """
        if self.prev_C_alpha is None:
            # First frame: fix s = +1 by convention.
            self.s = 1
            self.last_D = 1.0
            self.last_D_corr = 1.0
            self.last_det_alpha = 1.0
            self.last_det_beta = 1.0
            self.prev_C_alpha = C_alpha.copy()
            self.prev_C_beta = C_beta.copy()
            self.prev_n_alpha = n_alpha
            self.prev_n_beta = n_beta
            return self.s, 1.0

        # Build cross-geometry AO overlap approximation.
        if cross_overlap_mode == "current":
            S_cross = S_ao_current
        elif cross_overlap_mode == "previous":
            if S_ao_previous is None:
                raise ValueError("S_ao_previous required for 'previous' mode")
            S_cross = S_ao_previous
        elif cross_overlap_mode == "midpoint":
            if S_ao_previous is None:
                raise ValueError("S_ao_previous required for 'midpoint' mode")
            if S_ao_previous.shape != S_ao_current.shape:
                # Number of AOs changed (shouldn't happen) -> fall back.
                S_cross = S_ao_current
            else:
                S_cross = 0.5 * (S_ao_previous + S_ao_current)
        elif cross_overlap_mode == "odin":
            # S_ao_previous must be the pre-computed exact cross-geometry
            # overlap S^{n-1,n} provided by the caller (e.g., from ODIN).
            if S_ao_previous is None:
                raise ValueError(
                    "S_ao_previous must be the ODIN-computed cross-overlap "
                    "for cross_overlap_mode='odin'."
                )
            S_cross = S_ao_previous
        else:
            raise ValueError(f"Unknown cross_overlap_mode: {cross_overlap_mode}")

        # If the AO dimensions changed (e.g. different number of basis
        # functions across geometries), we cannot track -> keep previous s.
        if (S_cross.shape[0] != self.prev_C_alpha.shape[0] or
                S_cross.shape[1] != C_alpha.shape[0]):
            self.last_D = float("nan")
            self.last_D_corr = float("nan")
            self.prev_C_alpha = C_alpha.copy()
            self.prev_C_beta = C_beta.copy()
            self.prev_n_alpha = n_alpha
            self.prev_n_beta = n_beta
            return self.s, float("nan")

        # If the number of occupied orbitals changed (different multiplicity
        # or constraint state), keep previous s and reset reference.
        if n_alpha != self.prev_n_alpha or n_beta != self.prev_n_beta:
            self.last_D = float("nan")
            self.last_D_corr = float("nan")
            self.prev_C_alpha = C_alpha.copy()
            self.prev_C_beta = C_beta.copy()
            self.prev_n_alpha = n_alpha
            self.prev_n_beta = n_beta
            return self.s, float("nan")

        D_raw, det_a, det_b = compute_state_overlap_signed(
            self.prev_C_alpha, C_alpha,
            self.prev_C_beta, C_beta,
            n_alpha, n_beta,
            S_cross,
        )

        # Standard overlap-based phase tracking recursion:
        #   s_n = s_{n-1} * sgn[<Phi(t_{n-1}) | Phi(t_n)>_raw],
        # which guarantees the gauge-corrected overlap
        #   <Phi^corr(t_{n-1}) | Phi^corr(t_n)> = s_{n-1} * s_n * D_raw
        # = s_{n-1}^2 * |D_raw| = |D_raw| > 0,
        # i.e. the corrected diabatic-state trajectory is continuous.
        s_prev = self.s

        if D_raw == 0.0 or not np.isfinite(D_raw):
            # Degenerate / unreliable -> keep previous sign.
            sign_update = 1
        else:
            sign_update = 1 if D_raw > 0.0 else -1

        self.s = s_prev * sign_update

        # Bookkeeping: D_corr = <Phi^corr(t_{n-1}) | Phi(t_n)> = s_prev * D_raw.
        self.last_D = D_raw
        self.last_D_corr = s_prev * D_raw
        self.last_det_alpha = det_a
        self.last_det_beta = det_b

        # Store current MOs as the new reference for the next frame.
        self.prev_C_alpha = C_alpha.copy()
        self.prev_C_beta = C_beta.copy()
        self.prev_n_alpha = n_alpha
        self.prev_n_beta = n_beta

        return self.s, D_raw
