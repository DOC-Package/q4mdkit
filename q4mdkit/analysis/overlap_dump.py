from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np


def save_frame_mo_overlap_matrices_binary(
    output_dir: Path,
    frame_id: int,
    time_fs: float,
    tracker_A_votes: Optional[List[Any]] = None,
    tracker_B_votes: Optional[List[Any]] = None,
    selected_ref_index_A: int = -1,
    selected_ref_index_B: int = -1,
) -> Path:
    """Save per-frame occupied-space MO overlap matrices as a compressed archive."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    archive_data: Dict[str, Any] = {
        "frame_id": np.array(frame_id, dtype=int),
        "time_fs": np.array(time_fs, dtype=float),
        "selected_ref_index_A": np.array(selected_ref_index_A, dtype=int),
        "selected_ref_index_B": np.array(selected_ref_index_B, dtype=int),
    }

    _append_tracker_overlap_matrices(
        archive_data,
        state_label="A",
        reference_votes=tracker_A_votes,
    )
    _append_tracker_overlap_matrices(
        archive_data,
        state_label="B",
        reference_votes=tracker_B_votes,
    )

    archive_path = output_dir / f"frame_{frame_id:05d}.npz"
    np.savez_compressed(archive_path, **archive_data)
    return archive_path


def _append_tracker_overlap_matrices(
    archive_data: Dict[str, Any],
    state_label: str,
    reference_votes: Optional[List[Any]],
) -> None:
    valid_ref_indices: List[int] = []
    if reference_votes is None:
        archive_data[f"state_{state_label}_reference_indices"] = np.array([], dtype=int)
        return

    for ref_idx, vote in enumerate(reference_votes):
        alpha_overlap = getattr(vote.alpha, "overlap_matrix", None)
        beta_overlap = getattr(vote.beta, "overlap_matrix", None)
        if alpha_overlap is not None:
            archive_data[
                f"state_{state_label}_alpha_ref_{ref_idx:02d}"
            ] = np.asarray(alpha_overlap)
        if beta_overlap is not None:
            archive_data[
                f"state_{state_label}_beta_ref_{ref_idx:02d}"
            ] = np.asarray(beta_overlap)
        if alpha_overlap is not None or beta_overlap is not None:
            valid_ref_indices.append(ref_idx)
    archive_data[
        f"state_{state_label}_reference_indices"
    ] = np.asarray(valid_ref_indices, dtype=int)


def save_frame_overlap_matrices_binary(
    output_dir: Path,
    frame_id: int,
    time_fs: float,
    tracker_A_votes: Optional[List[Any]] = None,
    tracker_B_votes: Optional[List[Any]] = None,
    selected_ref_index_A: int = -1,
    selected_ref_index_B: int = -1,
) -> Path:
    """Backward-compatible wrapper for per-frame MO overlap archive saving."""
    return save_frame_mo_overlap_matrices_binary(
        output_dir=output_dir,
        frame_id=frame_id,
        time_fs=time_fs,
        tracker_A_votes=tracker_A_votes,
        tracker_B_votes=tracker_B_votes,
        selected_ref_index_A=selected_ref_index_A,
        selected_ref_index_B=selected_ref_index_B,
    )