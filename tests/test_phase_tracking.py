import numpy as np

from q4mdkit.analysis.phase_tracking import (
    StatePhaseTracker,
    choose_phase_continuity_override,
    choose_phase_sign_from_weighted_votes,
)


def test_state_phase_tracker_corrects_single_spin_sign_flip():
    tracker = StatePhaseTracker(name="A")
    C_alpha = np.array([[1.0], [0.0]])
    C_beta = np.zeros((2, 0))
    S_ao = np.eye(2)

    s0, d0 = tracker.update(C_alpha, C_beta, 1, 0, S_ao)
    assert s0 == 1
    assert d0 == 1.0

    s1, d1 = tracker.update(-C_alpha, C_beta, 1, 0, S_ao)
    assert s1 == -1
    assert np.isclose(d1, -1.0)
    assert np.isclose(tracker.last_D_corr, 1.0)
    assert np.isclose(tracker.last_sigma_min_alpha, 1.0)
    assert np.isclose(tracker.last_cond_alpha, 1.0)


def test_state_phase_tracker_uses_corrected_reference_for_next_frame():
    tracker = StatePhaseTracker(name="A")
    C_beta = np.zeros((2, 0))
    S_ao = np.eye(2)
    C_ref = np.eye(2)
    C_swap = C_ref[:, [1, 0]]

    tracker.update(C_ref, C_beta, 2, 0, S_ao)

    s1, d1 = tracker.update(C_swap, C_beta, 2, 0, S_ao)
    assert s1 == -1
    assert np.isclose(d1, -1.0)

    s2, d2 = tracker.update(C_ref, C_beta, 2, 0, S_ao)
    assert s2 == 1
    assert np.isclose(d2, 1.0)
    assert np.isclose(tracker.last_D_corr, 1.0)


def test_state_phase_tracker_uses_global_many_electron_sign():
    tracker = StatePhaseTracker(name="A")
    C_alpha = np.array([[1.0], [0.0]])
    C_beta = np.array([[0.0], [1.0]])
    S_ao = np.eye(2)

    tracker.update(C_alpha, C_beta, 1, 1, S_ao)
    s1, d1 = tracker.update(-C_alpha, -C_beta, 1, 1, S_ao)

    assert s1 == 1
    assert np.isclose(d1, 1.0)
    assert np.isclose(tracker.last_D_corr, 1.0)


def test_state_phase_tracker_looks_back_when_latest_sigma_is_low():
    tracker = StatePhaseTracker(
        name="A",
        vote_history=1,
        reference_history=5,
        sigma_accept_threshold=0.5,
    )
    C_beta = np.zeros((3, 0))
    S_ao = np.eye(3)
    C_ref = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]
    )
    C_low_sigma = np.array(
        [
            [1.0, 0.0],
            [0.0, -0.1],
            [0.0, np.sqrt(1.0 - 0.1**2)],
        ]
    )

    tracker.update(C_ref, C_beta, 2, 0, S_ao)
    s1, d1 = tracker.update(C_low_sigma, C_beta, 2, 0, S_ao)
    assert s1 == 1
    assert np.isclose(d1, -0.1)
    assert tracker.last_ambiguous is True

    s2, d2 = tracker.update(C_ref, C_beta, 2, 0, S_ao)

    assert s2 == 1
    assert np.isclose(d2, 1.0)
    assert tracker.last_selected_ref_index == 1
    assert tracker.last_diagnostic_selected_ref_index == 1
    assert tracker.last_lookback_used is True
    assert tracker.last_low_sigma_refs == 1
    assert len(tracker.last_reference_votes) == 2
    assert np.isclose(tracker.last_reference_votes[0].sigma_min, 0.1)
    assert np.isclose(tracker.last_reference_votes[1].sigma_min, 1.0)
    assert np.isclose(tracker.last_sigma_min_alpha, 1.0)
    assert tracker.last_truncated_alpha == 0
    assert tracker.last_ambiguous is False


def test_state_phase_tracker_invalid_low_primary_does_not_advance_reference():
    tracker = StatePhaseTracker(
        name="A",
        reference_history=5,
        sigma_accept_threshold=0.5,
        invalidate_low_primary_sigma=True,
    )
    C_beta = np.zeros((3, 0))
    S_ao = np.eye(3)
    C_ref = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]
    )
    C_low_sigma = np.array(
        [
            [1.0, 0.0],
            [0.0, -0.1],
            [0.0, np.sqrt(1.0 - 0.1**2)],
        ]
    )

    tracker.update(C_ref, C_beta, 2, 0, S_ao)
    s1, d1 = tracker.update(C_low_sigma, C_beta, 2, 0, S_ao)

    assert s1 == 1
    assert np.isnan(d1)
    assert tracker.last_invalid is True
    assert tracker.last_selected_ref_index == -1
    assert tracker.last_diagnostic_selected_ref_index == -1
    assert len(tracker.ref_history) == 1
    assert np.allclose(tracker.ref_history[0].occ_alpha, C_ref)

    s2, d2 = tracker.update(C_ref, C_beta, 2, 0, S_ao)

    assert s2 == 1
    assert np.isclose(d2, 1.0)
    assert tracker.last_invalid is False
    assert tracker.last_selected_ref_index == 0
    assert tracker.last_diagnostic_selected_ref_index == 0
    assert tracker.last_lookback_used is False


def test_state_phase_tracker_low_primary_uses_accepted_lookback_reference():
    tracker = StatePhaseTracker(
        name="A",
        reference_history=5,
        sigma_accept_threshold=0.5,
        invalidate_low_primary_sigma=False,
    )
    C_beta = np.zeros((3, 0))
    S_ao = np.eye(3)
    C_ref = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]
    )
    C_low_sigma = np.array(
        [
            [1.0, 0.0],
            [0.0, -0.1],
            [0.0, np.sqrt(1.0 - 0.1**2)],
        ]
    )

    tracker.update(C_ref, C_beta, 2, 0, S_ao)
    tracker.update(C_low_sigma, C_beta, 2, 0, S_ao)

    tracker.invalidate_low_primary_sigma = True
    _, d2 = tracker.update(C_ref, C_beta, 2, 0, S_ao)

    assert np.isclose(d2, 1.0)
    assert tracker.last_invalid is False
    assert tracker.last_selected_ref_index == 1
    assert tracker.last_diagnostic_selected_ref_index == 1
    assert tracker.last_lookback_used is True
    assert tracker.last_reference_votes[0].sigma_min < 0.5
    assert tracker.last_reference_votes[1].sigma_min >= 0.5


def test_weighted_phase_vote_marks_small_margin_as_ambiguous():
    chosen_sign, score, total_weight, margin, ambiguous = (
        choose_phase_sign_from_weighted_votes(
            signs=[-1, 1, 1],
            weights=[0.40, 0.30, 0.20],
            fallback_sign=1,
            ambiguity_ratio=0.15,
        )
    )

    assert chosen_sign == 1
    assert ambiguous is True
    assert np.isclose(score, 0.10)
    assert np.isclose(total_weight, 0.90)
    assert np.isclose(margin, score / total_weight)


def test_weighted_phase_vote_follows_strong_majority():
    chosen_sign, score, total_weight, margin, ambiguous = (
        choose_phase_sign_from_weighted_votes(
            signs=[-1, -1, 1],
            weights=[0.60, 0.25, 0.05],
            fallback_sign=1,
            ambiguity_ratio=0.15,
        )
    )

    assert chosen_sign == -1
    assert ambiguous is False
    assert score < 0.0
    assert np.isclose(total_weight, 0.90)
    assert margin > 0.15


def test_state_phase_tracker_flip_current_sign_updates_reference():
    tracker = StatePhaseTracker(name="A")
    C_alpha = np.array([[1.0], [0.0]])
    C_beta = np.zeros((2, 0))
    S_ao = np.eye(2)

    tracker.update(C_alpha, C_beta, 1, 0, S_ao)
    prev_occ = tracker.prev_occ_alpha.copy()

    tracker.flip_current_sign()

    assert tracker.s == -1
    assert np.isclose(tracker.last_D_corr, -1.0)
    assert np.allclose(tracker.prev_occ_alpha[:, 0], -prev_occ[:, 0])


def test_choose_phase_continuity_override_prefers_smooth_branch():
    choice = choose_phase_continuity_override(
        gauge=1,
        H_raw=-0.3834212921,
        S_raw=0.0044804521,
        prev_H_corr=0.6676052604,
        prev_S_corr=-0.0078009037,
        D_A_eff=0.9941756093,
        D_B_eff=0.9854903362,
        overlap_threshold=0.995,
    )

    assert choice == "B"


def test_choose_phase_continuity_override_skips_high_confidence_frames():
    choice = choose_phase_continuity_override(
        gauge=1,
        H_raw=-0.3834212921,
        S_raw=0.0044804521,
        prev_H_corr=0.6676052604,
        prev_S_corr=-0.0078009037,
        D_A_eff=0.999,
        D_B_eff=0.998,
        overlap_threshold=0.995,
    )

    assert choice is None
