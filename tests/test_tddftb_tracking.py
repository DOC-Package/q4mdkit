import numpy as np

from q4mdkit.analysis.tddftb import (
    Excitation,
    FrontierOrbitalTracker,
    apply_orbital_tracking_to_excitations,
)


def test_frontier_orbital_tracker_preserves_labels_across_permutation():
    tracker = FrontierOrbitalTracker(n_excitations=2, window=2, overlap_threshold=0.1)

    eigenvalues = np.array([-2.0, -1.0, 1.0, 2.0])
    occupations = np.array([1.0, 1.0, 0.0, 0.0])
    coeff_ref = np.eye(4)
    coeff_swap = coeff_ref[:, [0, 1, 3, 2]]

    first = tracker.update(eigenvalues, occupations, coeff_ref)
    second = tracker.update(eigenvalues, occupations, coeff_swap)

    assert first.orbital_label_map[3] == 3
    assert first.orbital_label_map[4] == 4
    assert second.orbital_label_map[3] == 4
    assert second.orbital_label_map[4] == 3
    assert second.valid is True
    assert np.isclose(second.sigma_min, 1.0)


def test_apply_orbital_tracking_to_excitations_reorders_by_tracked_transition():
    excitations = [
        Excitation(
            energy_ev=1.0,
            oscillator_strength=0.1,
            transition_from=2,
            transition_to=4,
        ),
        Excitation(
            energy_ev=1.1,
            oscillator_strength=0.2,
            transition_from=2,
            transition_to=3,
        ),
    ]
    tracker = FrontierOrbitalTracker(n_excitations=2, window=2, overlap_threshold=0.1)
    eigenvalues = np.array([-2.0, -1.0, 1.0, 2.0])
    occupations = np.array([1.0, 1.0, 0.0, 0.0])
    coeff_ref = np.eye(4)
    coeff_swap = coeff_ref[:, [0, 1, 3, 2]]

    tracker.update(eigenvalues, occupations, coeff_ref)
    result = tracker.update(eigenvalues, occupations, coeff_swap)
    reordered = apply_orbital_tracking_to_excitations(excitations, result)

    assert reordered[0].tracked_transition_to == 3
    assert reordered[1].tracked_transition_to == 4


def test_frontier_orbital_tracker_marks_low_overlap_as_invalid():
    tracker = FrontierOrbitalTracker(n_excitations=2, window=2, overlap_threshold=0.95)

    eigenvalues = np.array([-2.0, -1.0, 1.0, 2.0])
    occupations = np.array([1.0, 1.0, 0.0, 0.0])
    coeff_ref = np.eye(4)
    coeff_same = np.eye(4)
    cross_overlap = np.diag([1.0, 1.0, 1.0, 0.4])

    tracker.update(eigenvalues, occupations, coeff_ref)
    result = tracker.update(
        eigenvalues,
        occupations,
        coeff_same,
        cross_overlap=cross_overlap,
    )

    assert result.valid is False
    assert np.isclose(result.sigma_min, 0.4)


def test_frontier_orbital_tracker_uses_cross_overlap_matrix_for_assignment():
    tracker = FrontierOrbitalTracker(n_excitations=2, window=2, overlap_threshold=0.1)

    eigenvalues = np.array([-2.0, -1.0, 1.0, 2.0])
    occupations = np.array([1.0, 1.0, 0.0, 0.0])
    coeff_ref = np.eye(4)
    coeff_same = np.eye(4)
    cross_overlap = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )

    tracker.update(eigenvalues, occupations, coeff_ref)
    result = tracker.update(
        eigenvalues,
        occupations,
        coeff_same,
        cross_overlap=cross_overlap,
    )

    assert result.valid is True
    assert result.orbital_label_map[3] == 4
    assert result.orbital_label_map[4] == 3
    assert np.isclose(result.sigma_min, 1.0)