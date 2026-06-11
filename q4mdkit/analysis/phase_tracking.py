"""
Phase (gauge) tracking for diabatic states along an MD trajectory.

In CDFTB-CI calculations, each charge-localized diabatic state |Phi_f(t_n)>
on fragment f at time t_n is defined only up to an arbitrary sign s_f(t_n)=+-1
(the wavefunction is real-valued in DFTB). The coupling transforms as
    H_AB(t_n) -> s_A(t_n) s_B(t_n) H_AB(t_n),
so when each snapshot is optimized independently the raw H_AB(t) trajectory
can exhibit artificial sign flips that come purely from this gauge freedom.

This module fixes the gauge by overlap-based phase tracking. At each step the
sign of |Phi_f(t_n)> is chosen from occupied-space overlaps with recent
gauge-fixed references. The most recent frame is tried first. If its occupied
overlap has a small singular value, older references are tried until a
well-conditioned overlap is found. This suppresses isolated spurious sign
flips without deleting singular modes or consulting H_AB or S_AB.

The occupied orbitals at t_n are parallel-transported within their occupied
subspace so that the next frame is compared against smoothly varying
gauge-fixed references rather than against the raw canonical MOs.

Notes on the cross-geometry AO overlap S^{n-1,n}
------------------------------------------------
The exact quantity is

    S^{n-1,n}_{mu nu} = <chi_mu(t_{n-1}) | chi_nu(t_n)>,

i.e. the AO overlap between basis functions evaluated at the two different
geometries. For DFTB+ this is not directly available without a dedicated
"doubled-system" calculation. For sign tracking it is sufficient to
approximate

    S^{n-1,n}_{mu nu} ~= S_{mu nu}(t_n)         (default)

or, optionally, the midpoint average 1/2 [S(t_{n-1}) + S(t_n)]. For typical
MD time steps (<~ few fs) the displacement is small compared to atomic
orbital extents and this approximation reliably preserves the sign of the
determinant det O_f^sigma even when its magnitude is somewhat reduced.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np


def _occ_block(C: np.ndarray, n_occ: int) -> np.ndarray:
    """Return the occupied block of the MO coefficient matrix."""
    if n_occ <= 0:
        return np.zeros((C.shape[0], 0), dtype=C.dtype)
    return C[:, :n_occ]


def _signed_value_from_logabs(sign: float, logabs: float) -> float:
    """Convert ``slogdet`` output to a floating-point determinant value."""
    if sign == 0.0:
        return 0.0
    if not np.isfinite(sign) or not np.isfinite(logabs):
        return float("nan")
    tiny_log = np.log(np.finfo(float).tiny)
    max_log = np.log(np.finfo(float).max)
    if logabs < tiny_log:
        return 0.0
    if logabs > max_log:
        return float(np.sign(sign) * np.inf)
    return float(np.sign(sign) * np.exp(logabs))


def _effective_overlap_from_logabs(logabs: float, n_occ: int) -> float:
    """Return |D|^(1/N_occ) from log |D|."""
    if n_occ <= 0 or not np.isfinite(logabs):
        return float("nan")
    return float(np.exp(logabs / n_occ))


@dataclass
class _SpinChannelOverlap:
    """Cross-frame overlap diagnostics and transported occupied block."""
    det: float
    logabs: float
    sign: int
    full_det: float
    full_logabs: float
    full_sign: int
    rotation_det: int
    sigma_min: float
    sigma_max: float
    condition_number: float
    retained_sigma_min: float
    retained_count: int
    truncated_count: int
    overlap_matrix: np.ndarray
    transported_block: np.ndarray


@dataclass
class _ReferenceState:
    """Gauge-fixed occupied blocks expressed in a common AO basis."""
    occ_alpha: np.ndarray
    occ_beta: np.ndarray


@dataclass
class _ReferenceVote:
    """Overlap diagnostics for one history reference."""
    sign: int
    det: float
    logabs: float
    eff: float
    weight: float
    sigma_min: float
    alpha: _SpinChannelOverlap
    beta: _SpinChannelOverlap


def _compute_spin_channel_overlap(
    C_prev_occ: np.ndarray,
    C_curr_occ: np.ndarray,
    S_ao_cross: np.ndarray,
    corresponding_orbital_alignment: bool = True,
) -> _SpinChannelOverlap:
    """
    Compute a spin-channel occupied-space overlap and its parallel transport.

    When corresponding-orbital alignment is enabled, the current occupied
    block is rotated by the orthogonal Procrustes solution so that it is
    maximally aligned with the previous occupied block.
    """
    if C_prev_occ.shape[1] == 0:
        return _SpinChannelOverlap(
            det=1.0,
            logabs=0.0,
            sign=1,
            full_det=1.0,
            full_logabs=0.0,
            full_sign=1,
            rotation_det=1,
            sigma_min=1.0,
            sigma_max=1.0,
            condition_number=1.0,
            retained_sigma_min=1.0,
            retained_count=0,
            truncated_count=0,
            overlap_matrix=np.array([[1.0]]),
            transported_block=C_curr_occ.copy(),
        )

    overlap = C_prev_occ.T @ S_ao_cross @ C_curr_occ
    full_sign, full_logabs = np.linalg.slogdet(overlap)
    full_det = _signed_value_from_logabs(full_sign, full_logabs)

    U, singular_values, Vt = np.linalg.svd(overlap, full_matrices=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    if sigma_min <= 0.0:
        condition_number = float("inf")
    else:
        condition_number = float(sigma_max / sigma_min)

    if corresponding_orbital_alignment:
        rotation = Vt.T @ U.T
        rot_sign, _ = np.linalg.slogdet(rotation)
        rotation_det = 1 if rot_sign >= 0.0 else -1
        transported_block = C_curr_occ @ rotation
    else:
        rotation_det = 1
        transported_block = C_curr_occ.copy()

    sign = int(np.sign(full_sign))
    logabs = float(full_logabs)
    det = full_det
    retained_sigma_min = sigma_min
    retained_count = int(singular_values.size)

    return _SpinChannelOverlap(
        det=det,
        logabs=logabs,
        sign=sign,
        full_det=full_det,
        full_logabs=float(full_logabs),
        full_sign=int(np.sign(full_sign)),
        rotation_det=rotation_det,
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        condition_number=condition_number,
        retained_sigma_min=retained_sigma_min,
        retained_count=retained_count,
        truncated_count=0,
        overlap_matrix=overlap.copy(),
        transported_block=transported_block,
    )


def _overlap_vote_weight(
    alpha: _SpinChannelOverlap,
    beta: _SpinChannelOverlap,
) -> float:
    """
    Return a conservative confidence weight for sign voting.

    The sign of the many-electron determinant becomes unreliable when even one
    singular value gets small, so the minimum occupied-space singular value is
    used together with the determinant-based D_eff metric.
    """
    retained_occ = alpha.retained_count + beta.retained_count
    if retained_occ <= 0:
        return 0.0

    eff = _effective_overlap_from_logabs(alpha.logabs + beta.logabs, retained_occ)
    if not np.isfinite(eff):
        return 0.0

    sigma_candidates = []
    if alpha.retained_count > 0 and np.isfinite(alpha.retained_sigma_min):
        sigma_candidates.append(alpha.retained_sigma_min)
    if beta.retained_count > 0 and np.isfinite(beta.retained_sigma_min):
        sigma_candidates.append(beta.retained_sigma_min)
    if not sigma_candidates:
        return 0.0
    sigma_min = min(sigma_candidates)

    return float(max(0.0, min(eff, sigma_min)))


def _build_reference_vote(
    ref: _ReferenceState,
    C_curr_alpha: np.ndarray,
    C_curr_beta: np.ndarray,
    n_alpha: int,
    n_beta: int,
    S_ao_cross: np.ndarray,
    corresponding_orbital_alignment: bool = True,
) -> _ReferenceVote:
    """Compute the sign proposal from one transported history reference."""
    alpha = _compute_spin_channel_overlap(
        ref.occ_alpha,
        C_curr_alpha,
        S_ao_cross,
        corresponding_orbital_alignment=corresponding_orbital_alignment,
    )
    beta = _compute_spin_channel_overlap(
        ref.occ_beta,
        C_curr_beta,
        S_ao_cross,
        corresponding_orbital_alignment=corresponding_orbital_alignment,
    )

    logabs = alpha.logabs + beta.logabs
    sign = alpha.sign * beta.sign
    det = _signed_value_from_logabs(sign, logabs)
    eff = _effective_overlap_from_logabs(
        logabs,
        alpha.retained_count + beta.retained_count,
    )
    weight = _overlap_vote_weight(alpha, beta)
    sigma_min = min(alpha.sigma_min, beta.sigma_min)

    return _ReferenceVote(
        sign=sign,
        det=det,
        logabs=logabs,
        eff=eff,
        weight=weight,
        sigma_min=sigma_min,
        alpha=alpha,
        beta=beta,
    )


def _transport_reference_blocks(
    alpha: _SpinChannelOverlap,
    beta: _SpinChannelOverlap,
    target_sign: int,
) -> _ReferenceState:
    """
    Return the transported occupied blocks with the requested many-electron sign.
    """
    alpha_ref = alpha.transported_block
    beta_ref = beta.transported_block
    total_ref_sign = alpha.rotation_det * beta.rotation_det

    if target_sign not in (-1, 1):
        target_sign = 1

    if total_ref_sign != target_sign:
        if alpha_ref.shape[1] > 0:
            alpha_ref = alpha_ref.copy()
            alpha_ref[:, 0] *= -1.0
        elif beta_ref.shape[1] > 0:
            beta_ref = beta_ref.copy()
            beta_ref[:, 0] *= -1.0

    return _ReferenceState(alpha_ref.copy(), beta_ref.copy())


def choose_phase_sign_from_weighted_votes(
    signs: List[int],
    weights: List[float],
    fallback_sign: int = 1,
    ambiguity_ratio: float = 0.15,
    min_total_weight: float = 0.0,
) -> Tuple[int, float, float, float, bool]:
    """
    Choose the phase sign from weighted overlap votes.

    Parameters
    ----------
    signs, weights
        Parallel arrays of sign proposals (+1 or -1) and their non-negative
        confidence weights.
    fallback_sign
        Sign to keep when the vote is ambiguous.
    ambiguity_ratio
        Require |score| / total_weight >= ambiguity_ratio to accept the vote.
    min_total_weight
        If the total vote weight is below this value, the vote is treated as
        ambiguous regardless of the margin.

    Returns
    -------
    chosen_sign, score, total_weight, margin, ambiguous
    """
    valid_votes = []
    for sign, weight in zip(signs, weights):
        if sign not in (-1, 1):
            continue
        if not np.isfinite(weight) or weight <= 0.0:
            continue
        valid_votes.append((sign, float(weight)))

    if fallback_sign not in (-1, 1):
        fallback_sign = 1

    if not valid_votes:
        return fallback_sign, 0.0, 0.0, 0.0, True

    score = float(sum(sign * weight for sign, weight in valid_votes))
    total_weight = float(sum(weight for _, weight in valid_votes))
    if total_weight <= 0.0:
        return fallback_sign, score, total_weight, 0.0, True

    margin = float(abs(score) / total_weight)
    if total_weight < min_total_weight or margin < ambiguity_ratio:
        return fallback_sign, score, total_weight, margin, True

    chosen_sign = 1 if score > 0.0 else -1
    return chosen_sign, score, total_weight, margin, False


def compute_state_overlap_signed(
    C_prev_alpha: np.ndarray,
    C_curr_alpha: np.ndarray,
    C_prev_beta: np.ndarray,
    C_curr_beta: np.ndarray,
    n_alpha: int,
    n_beta: int,
    S_ao_cross: np.ndarray,
) -> Tuple[float, float, float, float]:
    """
    Compute <Phi(t_{n-1}) | Phi(t_n)> = det O_alpha * det O_beta.

    Parameters
    ----------
    C_prev_alpha, C_prev_beta : np.ndarray
        MO coefficients at t_{n-1}, shape (n_AO, n_MO).
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
    logabs_D : float
        Natural logarithm of |D|. This remains usable even when D underflows.
    """
    Ca_prev = _occ_block(C_prev_alpha, n_alpha)
    Ca_curr = _occ_block(C_curr_alpha, n_alpha)
    Cb_prev = _occ_block(C_prev_beta, n_beta)
    Cb_curr = _occ_block(C_curr_beta, n_beta)

    alpha = _compute_spin_channel_overlap(Ca_prev, Ca_curr, S_ao_cross)
    beta = _compute_spin_channel_overlap(Cb_prev, Cb_curr, S_ao_cross)

    total_sign = alpha.sign * beta.sign
    total_logabs = alpha.logabs + beta.logabs
    total_det = _signed_value_from_logabs(total_sign, total_logabs)

    return total_det, alpha.det, beta.det, total_logabs


@dataclass
class StatePhaseTracker:
    """
    Tracks the gauge factor s_f for a single diabatic state along an MD
    trajectory.

    Parameters
    ----------
    vote_history
        Maximum number of references used for backup vote diagnostics when no
        reference passes ``sigma_accept_threshold``. The sign itself is kept
        unchanged in that case and the frame is marked ambiguous.
    vote_decay
        Geometric decay applied to older references.
    vote_ambiguity_ratio
        If |vote_score| / total_vote_weight is below this threshold, the sign
        is treated as ambiguous and the previous sign is kept.
    reference_history
        Number of recent gauge-fixed references kept in memory. The latest
        reference is tried first; older references are tried only when the
        latest overlap has a small singular value.
    sigma_accept_threshold
        Minimum combined spin-channel singular value required to accept a
        reference for deterministic sign tracking. Set to 0.0 to always accept
        the latest reference.
    sigma_filtering_enabled
        If True, use ``sigma_accept_threshold`` to decide whether the latest
        reference is accepted, whether to look back to older references, and
        whether low-sigma frames can be invalidated. If False, always use the
        latest reference regardless of singular values.
    corresponding_orbital_alignment
        If True, align the current occupied orbitals to the previous occupied
        subspace by the corresponding-orbital SVD rotation before committing a
        new reference. If False, keep the raw occupied orbitals.
    invalidate_low_primary_sigma
        If True, the current frame is invalidated when none of the stored
        references reaches ``sigma_accept_threshold``. Older references are
        still tried before invalidating the frame.
    """

    name: str = ""
    vote_history: int = 1
    vote_decay: float = 0.6
    vote_ambiguity_ratio: float = 0.15
    reference_history: int = 5
    sigma_accept_threshold: float = 0.0
    sigma_filtering_enabled: bool = True
    corresponding_orbital_alignment: bool = True
    invalidate_low_primary_sigma: bool = False
    s: int = 1
    prev_occ_alpha: Optional[np.ndarray] = field(default=None, repr=False)
    prev_occ_beta: Optional[np.ndarray] = field(default=None, repr=False)
    prev_n_alpha: int = 0
    prev_n_beta: int = 0
    ref_history: List[_ReferenceState] = field(default_factory=list, repr=False)
    last_D: float = 1.0
    last_D_corr: float = 1.0
    last_logabs_D: float = 0.0
    last_det_alpha: float = 1.0
    last_det_beta: float = 1.0
    last_sigma_min_alpha: float = 1.0
    last_sigma_min_beta: float = 1.0
    last_cond_alpha: float = 1.0
    last_cond_beta: float = 1.0
    last_vote_score: float = 1.0
    last_vote_total_weight: float = 1.0
    last_vote_margin: float = 1.0
    last_vote_refs: int = 1
    last_ambiguous: bool = False
    last_retained_alpha: int = 0
    last_retained_beta: int = 0
    last_truncated_alpha: int = 0
    last_truncated_beta: int = 0
    last_selected_ref_index: int = 0
    last_diagnostic_selected_ref_index: int = 0
    last_lookback_used: bool = False
    last_low_sigma_refs: int = 0
    last_reference_refs: int = 1
    last_reference_votes: List[_ReferenceVote] = field(default_factory=list, repr=False)
    last_invalid: bool = False

    def reset(self) -> None:
        self.s = 1
        self.prev_occ_alpha = None
        self.prev_occ_beta = None
        self.prev_n_alpha = 0
        self.prev_n_beta = 0
        self.ref_history = []
        self.last_D = 1.0
        self.last_D_corr = 1.0
        self.last_logabs_D = 0.0
        self.last_det_alpha = 1.0
        self.last_det_beta = 1.0
        self.last_sigma_min_alpha = 1.0
        self.last_sigma_min_beta = 1.0
        self.last_cond_alpha = 1.0
        self.last_cond_beta = 1.0
        self.last_vote_score = 1.0
        self.last_vote_total_weight = 1.0
        self.last_vote_margin = 1.0
        self.last_vote_refs = 1
        self.last_ambiguous = False
        self.last_retained_alpha = 0
        self.last_retained_beta = 0
        self.last_truncated_alpha = 0
        self.last_truncated_beta = 0
        self.last_selected_ref_index = 0
        self.last_diagnostic_selected_ref_index = 0
        self.last_lookback_used = False
        self.last_low_sigma_refs = 0
        self.last_reference_refs = 1
        self.last_reference_votes = []
        self.last_invalid = False

    def _initialize_reference_history(
        self,
        C_alpha: np.ndarray,
        C_beta: np.ndarray,
        n_alpha: int,
        n_beta: int,
    ) -> None:
        occ_alpha = _occ_block(C_alpha, n_alpha).copy()
        occ_beta = _occ_block(C_beta, n_beta).copy()
        ref = _ReferenceState(occ_alpha, occ_beta)

        self.prev_occ_alpha = occ_alpha.copy()
        self.prev_occ_beta = occ_beta.copy()
        self.prev_n_alpha = n_alpha
        self.prev_n_beta = n_beta
        self.ref_history = [ref]

    def _build_cross_overlap(
        self,
        S_ao_current: np.ndarray,
        S_ao_previous: Optional[np.ndarray],
        cross_overlap_mode: str,
    ) -> np.ndarray:
        """Build the cross-geometry AO overlap approximation."""
        if cross_overlap_mode == "current":
            return S_ao_current
        if cross_overlap_mode == "previous":
            if S_ao_previous is None:
                raise ValueError("S_ao_previous required for 'previous' mode")
            return S_ao_previous
        if cross_overlap_mode == "midpoint":
            if S_ao_previous is None:
                raise ValueError("S_ao_previous required for 'midpoint' mode")
            if S_ao_previous.shape != S_ao_current.shape:
                return S_ao_current
            return 0.5 * (S_ao_previous + S_ao_current)
        if cross_overlap_mode == "odin":
            if S_ao_previous is None:
                raise ValueError(
                    "S_ao_previous must be the ODIN-computed cross-overlap "
                    "for cross_overlap_mode='odin'."
                )
            return S_ao_previous
        raise ValueError(f"Unknown cross_overlap_mode: {cross_overlap_mode}")

    def _set_reset_diagnostics(self) -> None:
        """Mark the latest diagnostics as invalid after a tracker reset."""
        self.last_D = float("nan")
        self.last_D_corr = float("nan")
        self.last_logabs_D = float("nan")
        self.last_det_alpha = float("nan")
        self.last_det_beta = float("nan")
        self.last_sigma_min_alpha = float("nan")
        self.last_sigma_min_beta = float("nan")
        self.last_cond_alpha = float("nan")
        self.last_cond_beta = float("nan")
        self.last_vote_score = float("nan")
        self.last_vote_total_weight = float("nan")
        self.last_vote_margin = float("nan")
        self.last_vote_refs = 0
        self.last_ambiguous = False
        self.last_retained_alpha = 0
        self.last_retained_beta = 0
        self.last_truncated_alpha = 0
        self.last_truncated_beta = 0
        self.last_selected_ref_index = 0
        self.last_diagnostic_selected_ref_index = 0
        self.last_lookback_used = False
        self.last_low_sigma_refs = 0
        self.last_reference_refs = 0
        self.last_reference_votes = []
        self.last_invalid = False

    def _build_cross_overlap_history(
        self,
        S_ao_current: np.ndarray,
        S_ao_previous: Optional[np.ndarray],
        cross_overlap_mode: str,
        n_refs: int,
        S_ao_cross_history: Optional[List[Optional[np.ndarray]]],
    ) -> List[np.ndarray]:
        """
        Build one AO cross-overlap matrix per stored reference.

        ``S_ao_cross_history`` is used by exact cross-overlap providers such as
        ODIN, where each history reference has a different geometry.
        """
        if S_ao_cross_history is None:
            S_cross = self._build_cross_overlap(
                S_ao_current,
                S_ao_previous,
                cross_overlap_mode,
            )
            return [S_cross for _ in range(n_refs)]

        cross_history: List[np.ndarray] = []
        for idx in range(n_refs):
            S_cross = None
            if idx < len(S_ao_cross_history):
                S_cross = S_ao_cross_history[idx]
            if S_cross is None:
                S_cross = S_ao_current
            cross_history.append(S_cross)
        return cross_history

    def update(
        self,
        C_alpha: np.ndarray,
        C_beta: np.ndarray,
        n_alpha: int,
        n_beta: int,
        S_ao_current: np.ndarray,
        S_ao_previous: Optional[np.ndarray] = None,
        cross_overlap_mode: str = "current",
        S_ao_cross_history: Optional[List[Optional[np.ndarray]]] = None,
        commit: bool = True,
    ) -> Tuple[int, float]:
        """
        Update gauge factor using overlap with transported history references.

        Returns
        -------
        s : int
            Updated gauge factor s_f(t_n) in {-1, +1}.
        D : float
            Selected reference overlap <Phi_f^corr(t_ref) | Phi_f(t_n)> before
            applying the current gauge factor s_f(t_n). On the first call this
            is 1.0 by convention.
        """
        if self.prev_occ_alpha is None:
            self.s = 1
            self.last_D = 1.0
            self.last_D_corr = 1.0
            self.last_logabs_D = 0.0
            self.last_det_alpha = 1.0
            self.last_det_beta = 1.0
            self.last_sigma_min_alpha = 1.0
            self.last_sigma_min_beta = 1.0
            self.last_cond_alpha = 1.0
            self.last_cond_beta = 1.0
            self.last_vote_score = 1.0
            self.last_vote_total_weight = 1.0
            self.last_vote_margin = 1.0
            self.last_vote_refs = 1
            self.last_ambiguous = False
            self.last_retained_alpha = n_alpha
            self.last_retained_beta = n_beta
            self.last_truncated_alpha = 0
            self.last_truncated_beta = 0
            self.last_selected_ref_index = 0
            self.last_diagnostic_selected_ref_index = 0
            self.last_lookback_used = False
            self.last_low_sigma_refs = 0
            self.last_reference_refs = 1
            self.last_reference_votes = []
            self.last_invalid = False
            self._initialize_reference_history(C_alpha, C_beta, n_alpha, n_beta)
            return self.s, 1.0

        if n_alpha != self.prev_n_alpha or n_beta != self.prev_n_beta:
            self._set_reset_diagnostics()
            self._initialize_reference_history(C_alpha, C_beta, n_alpha, n_beta)
            return self.s, float("nan")

        C_curr_alpha = _occ_block(C_alpha, n_alpha)
        C_curr_beta = _occ_block(C_beta, n_beta)
        max_refs = max(1, int(self.reference_history))
        refs = self.ref_history[:max_refs]
        if not refs:
            refs = [_ReferenceState(self.prev_occ_alpha.copy(), self.prev_occ_beta.copy())]

        S_crosses = self._build_cross_overlap_history(
            S_ao_current,
            S_ao_previous,
            cross_overlap_mode,
            len(refs),
            S_ao_cross_history,
        )
        if not S_crosses:
            self._set_reset_diagnostics()
            self._initialize_reference_history(C_alpha, C_beta, n_alpha, n_beta)
            return self.s, float("nan")

        first_cross = S_crosses[0]
        if (first_cross.shape[0] != refs[0].occ_alpha.shape[0] or
                first_cross.shape[1] != C_alpha.shape[0]):
            self._set_reset_diagnostics()
            self._initialize_reference_history(C_alpha, C_beta, n_alpha, n_beta)
            return self.s, float("nan")

        votes = [
            _build_reference_vote(
                ref,
                C_curr_alpha,
                C_curr_beta,
                n_alpha,
                n_beta,
                S_cross,
                corresponding_orbital_alignment=self.corresponding_orbital_alignment,
            )
            for ref, S_cross in zip(refs, S_crosses)
        ]

        threshold = (
            max(0.0, float(self.sigma_accept_threshold))
            if self.sigma_filtering_enabled
            else 0.0
        )
        selected_idx: Optional[int] = None
        if threshold <= 0.0:
            selected_idx = 0
        else:
            for idx, vote in enumerate(votes):
                if np.isfinite(vote.sigma_min) and vote.sigma_min >= threshold:
                    selected_idx = idx
                    break

        diagnostic_selected_idx = selected_idx if selected_idx is not None else -1
        no_accepted_reference = selected_idx is None
        invalid_frame = bool(
            self.sigma_filtering_enabled
            and self.invalidate_low_primary_sigma
            and no_accepted_reference
        )
        if selected_idx is None:
            weights = [
                vote.weight if np.isfinite(vote.weight) else -1.0
                for vote in votes
            ]
            selected_idx = int(np.argmax(weights))

        selected_vote = votes[selected_idx]

        s_prev = self.s
        D_raw = selected_vote.det
        last_logabs_D = selected_vote.logabs

        max_vote_refs = max(1, int(self.vote_history))
        vote_signs: List[int] = []
        vote_weights: List[float] = []
        for idx, vote in enumerate(votes[:max_vote_refs]):
            if vote.sign == 0:
                continue
            weight = float(vote.weight * (self.vote_decay ** idx))
            vote_signs.append(vote.sign)
            vote_weights.append(weight)

        if no_accepted_reference:
            _, vote_score, vote_total_weight, vote_margin, _ = (
                choose_phase_sign_from_weighted_votes(
                    vote_signs,
                    vote_weights,
                    fallback_sign=s_prev,
                    ambiguity_ratio=self.vote_ambiguity_ratio,
                )
            )
            chosen_sign = s_prev
            ambiguous = True
        else:
            chosen_sign = selected_vote.sign if selected_vote.sign in (-1, 1) else s_prev
            vote_score = float(chosen_sign * selected_vote.weight)
            vote_total_weight = float(selected_vote.weight)
            vote_margin = 1.0 if vote_total_weight > 0.0 else 0.0
            ambiguous = False

        if invalid_frame:
            D_raw = float("nan")
            last_logabs_D = float("nan")
            D_corr = float("nan")
        else:
            D_corr = (
                chosen_sign * D_raw
                if np.isfinite(D_raw)
                else float("nan")
            )

        self.last_D = D_raw
        self.last_D_corr = D_corr
        self.last_logabs_D = last_logabs_D
        self.last_det_alpha = selected_vote.alpha.det
        self.last_det_beta = selected_vote.beta.det
        self.last_sigma_min_alpha = selected_vote.alpha.sigma_min
        self.last_sigma_min_beta = selected_vote.beta.sigma_min
        self.last_cond_alpha = selected_vote.alpha.condition_number
        self.last_cond_beta = selected_vote.beta.condition_number
        self.last_vote_score = vote_score
        self.last_vote_total_weight = vote_total_weight
        self.last_vote_margin = vote_margin
        self.last_vote_refs = min(len(votes), max_vote_refs)
        self.last_ambiguous = ambiguous
        self.last_retained_alpha = selected_vote.alpha.retained_count
        self.last_retained_beta = selected_vote.beta.retained_count
        self.last_truncated_alpha = 0
        self.last_truncated_beta = 0
        self.last_selected_ref_index = -1 if invalid_frame else selected_idx
        self.last_diagnostic_selected_ref_index = diagnostic_selected_idx
        self.last_lookback_used = (not invalid_frame) and selected_idx > 0
        self.last_low_sigma_refs = 1 if invalid_frame else (
            selected_idx if not no_accepted_reference else len(votes)
        )
        self.last_reference_refs = len(votes)
        self.last_reference_votes = votes
        self.last_invalid = invalid_frame

        if invalid_frame or not commit:
            return self.s if invalid_frame else chosen_sign, D_raw

        self.s = chosen_sign
        current_ref = _transport_reference_blocks(
            selected_vote.alpha,
            selected_vote.beta,
            target_sign=self.s,
        )
        self.ref_history = [current_ref] + refs[:max(0, max_refs - 1)]
        self.prev_occ_alpha = self.ref_history[0].occ_alpha.copy()
        self.prev_occ_beta = self.ref_history[0].occ_beta.copy()
        self.prev_n_alpha = n_alpha
        self.prev_n_beta = n_beta

        return self.s, D_raw

    def flip_current_sign(self) -> None:
        """
        Flip the sign of the current corrected-state reference in place.

        This is used only by optional higher-level continuity heuristics.
        """
        self.s *= -1
        self.last_D_corr *= -1.0
        self.last_vote_score *= -1.0
        if self.ref_history:
            ref0 = self.ref_history[0]
            if ref0.occ_alpha.shape[1] > 0:
                ref0 = _ReferenceState(ref0.occ_alpha.copy(), ref0.occ_beta.copy())
                ref0.occ_alpha[:, 0] *= -1.0
            elif ref0.occ_beta.shape[1] > 0:
                ref0 = _ReferenceState(ref0.occ_alpha.copy(), ref0.occ_beta.copy())
                ref0.occ_beta[:, 0] *= -1.0
            self.ref_history[0] = ref0
            self.prev_occ_alpha = ref0.occ_alpha.copy()
            self.prev_occ_beta = ref0.occ_beta.copy()
            return
        if self.prev_occ_alpha is not None and self.prev_occ_alpha.shape[1] > 0:
            self.prev_occ_alpha = self.prev_occ_alpha.copy()
            self.prev_occ_alpha[:, 0] *= -1.0
        elif self.prev_occ_beta is not None and self.prev_occ_beta.shape[1] > 0:
            self.prev_occ_beta = self.prev_occ_beta.copy()
            self.prev_occ_beta[:, 0] *= -1.0


def choose_phase_continuity_override(
    gauge: int,
    H_raw: float,
    S_raw: float,
    prev_H_corr: float,
    prev_S_corr: float,
    D_A_eff: float,
    D_B_eff: float,
    overlap_threshold: float = 0.995,
    score_margin: float = 1.0e-8,
) -> Optional[str]:
    """
    Return which state sign to flip ("A" or "B") to preserve gauge continuity.

    The override is disabled by default at the caller level. It remains
    available as an optional post-processing heuristic when overlap evidence is
    weak and a smoother gauge branch is desired.
    """
    if gauge not in (-1, 1):
        raise ValueError(f"Gauge must be +/-1, got {gauge}")
    vals = (H_raw, S_raw, prev_H_corr, prev_S_corr)
    if not all(np.isfinite(v) for v in vals):
        return None

    conf_a = D_A_eff if np.isfinite(D_A_eff) else float("inf")
    conf_b = D_B_eff if np.isfinite(D_B_eff) else float("inf")
    if min(conf_a, conf_b) > overlap_threshold:
        return None

    def continuity_score(test_gauge: int) -> float:
        H_corr = test_gauge * H_raw
        S_corr = test_gauge * S_raw
        h_scale = max(abs(prev_H_corr), abs(H_corr), 1.0e-12)
        s_scale = max(abs(prev_S_corr), abs(S_corr), 1.0e-12)
        return abs(H_corr - prev_H_corr) / h_scale + abs(S_corr - prev_S_corr) / s_scale

    score_same = continuity_score(gauge)
    score_flip = continuity_score(-gauge)
    if score_flip + score_margin >= score_same:
        return None

    return "A" if conf_a <= conf_b else "B"
