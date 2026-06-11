#!/usr/bin/env python3
"""
Analyze per-frame MO overlap matrices saved as compressed NumPy archives.

This script is intentionally kept outside the q4mdkit package so it can be
copied and customized independently for one-off analyses.

Usage:
    python analyze_overlap_matrices.py
    python analyze_overlap_matrices.py output_phase_test/overlap_matrices
    python analyze_overlap_matrices.py output_phase_test/overlap_matrices \
        --output overlap_summary.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def iter_overlap_archives(input_dir: Path) -> List[Path]:
    """Return overlap archive files sorted by frame index in the filename."""
    return sorted(input_dir.glob("frame_*.npz"))


def iter_archive_matrices(data: np.lib.npyio.NpzFile) -> Iterable[tuple[str, np.ndarray, int]]:
    """Yield matrix label, matrix data, and reference index for stored overlaps."""
    for matrix_key in sorted(data.files):
        if matrix_key in {
            "frame_id",
            "time_fs",
            "selected_ref_index_A",
            "selected_ref_index_B",
            "state_A_reference_indices",
            "state_B_reference_indices",
        }:
            continue
        matrix = data[matrix_key]
        if matrix.ndim != 2:
            continue
        ref_idx = -1
        if "_ref_" in matrix_key:
            try:
                ref_idx = int(matrix_key.rsplit("_ref_", 1)[1])
            except ValueError:
                ref_idx = -1
        yield matrix_key, matrix, ref_idx


def matrix_metrics(matrix: np.ndarray) -> dict[str, float]:
    """Compute compact diagnostics for one overlap matrix."""
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(np.max(singular_values))
    sigma_min = float(np.min(singular_values))
    condition_number = float(np.inf) if sigma_min == 0.0 else sigma_max / sigma_min
    trace_value = np.trace(matrix)
    symmetry_error = float(np.linalg.norm(matrix - matrix.T))
    det_sign, logabs_det = np.linalg.slogdet(matrix)
    if det_sign == 0.0:
        determinant = 0.0
    else:
        determinant = float(det_sign * np.exp(logabs_det))
    return {
        "determinant": determinant,
        "det_sign": float(det_sign),
        "logabs_det": float(logabs_det),
        "fro_norm": float(np.linalg.norm(matrix)),
        "max_abs": float(np.max(np.abs(matrix))),
        "trace_real": float(np.real(trace_value)),
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "condition_number": condition_number,
        "symmetry_error": symmetry_error,
    }


def smallest_singular_mode_metrics(
    matrix: np.ndarray,
    n_components: int = 5,
) -> dict[str, float | int | str]:
    """Return dominant components of the smallest-singular-value mode."""
    U, singular_values, Vt = np.linalg.svd(matrix, full_matrices=False)
    left_vec = U[:, -1]
    right_vec = Vt[-1, :]

    left_order = np.argsort(np.abs(left_vec))[::-1][:n_components]
    right_order = np.argsort(np.abs(right_vec))[::-1][:n_components]

    def _format_components(indices: np.ndarray, vector: np.ndarray) -> str:
        return ",".join(
            f"{int(idx)}:{float(vector[idx]):+.6f}"
            for idx in indices
        )

    return {
        "sigma_min": float(singular_values[-1]),
        "left_top_index": int(left_order[0]),
        "left_top_abs": float(np.abs(left_vec[left_order[0]])),
        "right_top_index": int(right_order[0]),
        "right_top_abs": float(np.abs(right_vec[right_order[0]])),
        "left_top_components": _format_components(left_order, left_vec),
        "right_top_components": _format_components(right_order, right_vec),
    }


def analyze_archive(archive_path: Path) -> List[dict[str, float | int | str]]:
    """Extract one summary row for each saved overlap matrix in the archive."""
    rows: List[dict[str, float | int | str]] = []
    with np.load(archive_path) as data:
        frame_id = int(data["frame_id"])
        time_fs = float(data["time_fs"])

        for matrix_key, matrix, ref_idx in iter_archive_matrices(data):
            rows.append(
                {
                    "frame_id": frame_id,
                    "time_fs": time_fs,
                    "matrix_label": matrix_key,
                    "reference_index": ref_idx,
                    "n_rows": int(matrix.shape[0]),
                    "n_cols": int(matrix.shape[1]),
                    **matrix_metrics(matrix),
                }
            )
    return rows


def save_archive_heatmaps(archive_path: Path, output_dir: Path) -> int:
    """Save one heatmap PNG per overlap matrix stored in an archive."""
    output_dir.mkdir(parents=True, exist_ok=True)
    saved = 0
    with np.load(archive_path) as data:
        frame_id = int(data["frame_id"])
        time_fs = float(data["time_fs"])
        for matrix_key, matrix, _ in iter_archive_matrices(data):
            vmax = float(np.max(np.abs(matrix)))
            if vmax == 0.0:
                vmax = 1.0
            fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
            image = ax.imshow(
                matrix,
                cmap="coolwarm",
                origin="lower",
                interpolation="nearest",
                vmin=-vmax,
                vmax=vmax,
            )
            ax.set_title(f"frame {frame_id}  t={time_fs:.1f} fs\n{matrix_key}")
            ax.set_xlabel("Current occupied MO index")
            ax.set_ylabel("Reference occupied MO index")
            fig.colorbar(image, ax=ax, shrink=0.85, label="Overlap")
            png_path = output_dir / f"frame_{frame_id:05d}_{matrix_key}.png"
            fig.savefig(png_path, dpi=180)
            plt.close(fig)
            saved += 1
    return saved


def write_summary(rows: Iterable[dict[str, float | int | str]], output_path: Path) -> None:
    """Write a whitespace-delimited summary table."""
    header = (
        "# frame_id time_fs matrix_label reference_index n_rows n_cols "
        "determinant det_sign logabs_det fro_norm max_abs trace_real sigma_min sigma_max "
        "condition_number symmetry_error"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(
                f"{row['frame_id']:8d} "
                f"{row['time_fs']:12.6f} "
                f"{str(row['matrix_label']):>18s} "
                f"{row['reference_index']:8d} "
                f"{row['n_rows']:6d} "
                f"{row['n_cols']:6d} "
                f"{row['determinant']:16.10f} "
                f"{row['det_sign']:8.1f} "
                f"{row['logabs_det']:16.10f} "
                f"{row['fro_norm']:16.10f} "
                f"{row['max_abs']:16.10f} "
                f"{row['trace_real']:16.10f} "
                f"{row['sigma_min']:16.10f} "
                f"{row['sigma_max']:16.10f} "
                f"{row['condition_number']:16.10f} "
                f"{row['symmetry_error']:16.10f}\n"
            )


def write_smallest_singular_modes(
    rows: Iterable[dict[str, float | int | str]],
    output_path: Path,
) -> None:
    """Write dominant-component diagnostics for the smallest singular mode."""
    header = (
        "# frame_id time_fs matrix_label reference_index sigma_min "
        "left_top_index left_top_abs right_top_index right_top_abs "
        "left_top_components right_top_components"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(
                f"{row['frame_id']:8d} "
                f"{row['time_fs']:12.6f} "
                f"{str(row['matrix_label']):>18s} "
                f"{row['reference_index']:8d} "
                f"{row['sigma_min']:16.10f} "
                f"{row['left_top_index']:8d} "
                f"{row['left_top_abs']:16.10f} "
                f"{row['right_top_index']:8d} "
                f"{row['right_top_abs']:16.10f} "
                f"{row['left_top_components']:>64s} "
                f"{row['right_top_components']:>64s}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze MO overlap matrices saved by run_cdftbci_analysis"
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="overlap_matrices",
        help="Directory containing frame_*.npz overlap archives",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="overlap_matrix_summary.dat",
        help="Output summary file path",
    )
    parser.add_argument(
        "--heatmap-dir",
        default="overlap_heatmaps",
        help="Directory for per-matrix heatmap PNG files",
    )
    parser.add_argument(
        "--smin-output",
        default="overlap_smin_components.dat",
        help="Output file path for smallest-singular-mode component diagnostics",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise SystemExit(f"Input directory not found: {input_dir}")

    archive_paths = iter_overlap_archives(input_dir)
    if not archive_paths:
        raise SystemExit(f"No frame_*.npz files found in: {input_dir}")

    summary_rows: List[dict[str, float | int | str]] = []
    smin_rows: List[dict[str, float | int | str]] = []
    for archive_path in archive_paths:
        archive_rows = analyze_archive(archive_path)
        summary_rows.extend(archive_rows)
        with np.load(archive_path) as data:
            frame_id = int(data["frame_id"])
            time_fs = float(data["time_fs"])
            for matrix_key, matrix, ref_idx in iter_archive_matrices(data):
                smin_rows.append(
                    {
                        "frame_id": frame_id,
                        "time_fs": time_fs,
                        "matrix_label": matrix_key,
                        "reference_index": ref_idx,
                        **smallest_singular_mode_metrics(matrix),
                    }
                )

    output_path = Path(args.output).resolve()
    write_summary(summary_rows, output_path)
    smin_output_path = Path(args.smin_output).resolve()
    write_smallest_singular_modes(smin_rows, smin_output_path)

    heatmap_dir = Path(args.heatmap_dir).resolve()
    heatmap_count = 0
    for archive_path in archive_paths:
        heatmap_count += save_archive_heatmaps(archive_path, heatmap_dir)

    print(f"Analyzed {len(archive_paths)} overlap archives from {input_dir}")
    print(f"Wrote {len(summary_rows)} summary rows to {output_path}")
    print(f"Wrote {len(smin_rows)} smallest-singular-mode rows to {smin_output_path}")
    print(f"Wrote {heatmap_count} heatmaps to {heatmap_dir}")


if __name__ == "__main__":
    main()