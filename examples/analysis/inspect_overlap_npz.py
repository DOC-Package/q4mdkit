#!/usr/bin/env python3
"""
Load overlap matrices from an NPZ file, print determinant/singular values,
and save heatmaps.

Usage:
    python inspect_overlap_npz.py overlap_matrices/frame_09094.npz
    python inspect_overlap_npz.py overlap_matrices/frame_09094.npz --key state_B_beta_ref_00
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def iter_matrix_items(data: np.lib.npyio.NpzFile, key: str | None) -> list[tuple[str, np.ndarray]]:
    """Return the requested matrix items from the NPZ archive."""
    if key is not None:
        matrix = data[key]
        if matrix.ndim != 2:
            raise ValueError(f"{key} is not a 2D matrix")
        return [(key, matrix)]

    items: list[tuple[str, np.ndarray]] = []
    for name in sorted(data.files):
        matrix = data[name]
        if matrix.ndim == 2:
            items.append((name, matrix))
    return items


def save_heatmap(matrix: np.ndarray, title: str, output_path: Path) -> None:
    """Save a single heatmap PNG for a matrix."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
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
    ax.set_title(title)
    ax.set_xlabel("Column index")
    ax.set_ylabel("Row index")
    fig.colorbar(image, ax=ax, shrink=0.85, label="Overlap")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect overlap matrices stored in an NPZ file")
    parser.add_argument("input_npz", help="Path to the NPZ file")
    parser.add_argument("--key", help="Specific matrix key to inspect")
    parser.add_argument(
        "--heatmap-dir",
        default="overlap_heatmaps",
        help="Directory for output heatmap PNG files",
    )
    args = parser.parse_args()

    input_path = Path(args.input_npz).resolve()
    if not input_path.exists():
        raise SystemExit(f"NPZ file not found: {input_path}")

    heatmap_dir = Path(args.heatmap_dir).resolve()

    with np.load(input_path) as data:
        frame_id = data["frame_id"] if "frame_id" in data else None
        time_fs = data["time_fs"] if "time_fs" in data else None
        matrix_items = iter_matrix_items(data, args.key)

        if not matrix_items:
            raise SystemExit("No 2D matrices found in the NPZ file")

        for name, matrix in matrix_items:
            determinant = np.linalg.det(matrix)
            singular_values = np.linalg.svd(matrix, compute_uv=False)
            sigma_max = float(np.max(singular_values))
            sigma_min = float(np.min(singular_values))
            print(f"matrix: {name}")
            if frame_id is not None:
                print(f"  frame_id: {int(frame_id)}")
            if time_fs is not None:
                print(f"  time_fs: {float(time_fs):.6f}")
            print(f"  shape: {matrix.shape}")
            print(f"  determinant: {determinant:+.16e}")
            print(f"  sigma_max: {sigma_max:.16e}")
            print(f"  sigma_min: {sigma_min:.16e}")

            png_name = f"{input_path.stem}_{name}.png"
            title = name
            if frame_id is not None and time_fs is not None:
                title = f"frame {int(frame_id)}  t={float(time_fs):.1f} fs\n{name}"
            save_heatmap(matrix, title, heatmap_dir / png_name)
            print(f"  heatmap: {heatmap_dir / png_name}")
            print()


if __name__ == "__main__":
    main()