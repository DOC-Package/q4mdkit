import numpy as np

from q4mdkit.analysis.overlap_dump import save_frame_overlap_matrices_binary


class _FakeSpinOverlap:
    def __init__(self, matrix):
        self.overlap_matrix = np.array(matrix)


class _FakeVote:
    def __init__(self, alpha, beta):
        self.alpha = _FakeSpinOverlap(alpha)
        self.beta = _FakeSpinOverlap(beta)


def test_save_frame_overlap_matrices_binary_saves_mo_overlaps(tmp_path):
    output_path = save_frame_overlap_matrices_binary(
        tmp_path,
        frame_id=12,
        time_fs=48.0,
        tracker_A_votes=[
            _FakeVote([[0.9, 0.0], [0.0, 0.8]], [[0.7]]),
            _FakeVote([[0.6, 0.1], [0.1, 0.5]], [[0.4]]),
        ],
        tracker_B_votes=[_FakeVote([[0.3]], [[0.2, 0.0], [0.0, 0.1]])],
        selected_ref_index_A=1,
        selected_ref_index_B=0,
    )

    assert output_path.name == "frame_00012.npz"
    data = np.load(output_path)

    assert int(data["frame_id"]) == 12
    assert np.isclose(float(data["time_fs"]), 48.0)
    assert int(data["selected_ref_index_A"]) == 1
    assert int(data["selected_ref_index_B"]) == 0
    assert np.array_equal(data["state_A_reference_indices"], np.array([0, 1]))
    assert np.array_equal(data["state_B_reference_indices"], np.array([0]))
    assert np.allclose(data["state_A_alpha_ref_00"], np.array([[0.9, 0.0], [0.0, 0.8]]))
    assert np.allclose(data["state_A_beta_ref_00"], np.array([[0.7]]))
    assert np.allclose(data["state_A_alpha_ref_01"], np.array([[0.6, 0.1], [0.1, 0.5]]))
    assert np.allclose(data["state_B_beta_ref_00"], np.array([[0.2, 0.0], [0.0, 0.1]]))


def test_save_frame_overlap_matrices_binary_handles_missing_votes(tmp_path):
    output_path = save_frame_overlap_matrices_binary(
        tmp_path,
        frame_id=3,
        time_fs=12.0,
        tracker_A_votes=None,
        tracker_B_votes=None,
    )

    data = np.load(output_path)

    assert np.array_equal(data["state_A_reference_indices"], np.array([], dtype=int))
    assert np.array_equal(data["state_B_reference_indices"], np.array([], dtype=int))
    assert "state_A_alpha_ref_00" not in data.files