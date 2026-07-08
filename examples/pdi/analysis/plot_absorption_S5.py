#!/usr/bin/env python3
"""Plot cumulant absorption output file (absorption_S5.dat by default)."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, nargs="?", default=Path("absorption_S5.dat"))
    parser.add_argument("--out", type=Path, default=Path("absorption_S5.png"))
    parser.add_argument("--energy-min", type=float, default=0.0)
    parser.add_argument("--energy-max", type=float, default=6.0)
    parser.add_argument("--wavelength-min", type=float, default=200.0)
    parser.add_argument("--wavelength-max", type=float, default=800.0)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    path = args.file
    if not path.exists():
        raise SystemExit(f"Input file not found: {path}")

    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] < 3:
        raise SystemExit("Absorption file must have at least three columns: energy, wavelength, intensity")
    energy = data[:, 0]
    wavelength = data[:, 1]
    intensity = data[:, 2]
    norm = data[:, 3] if data.shape[1] >= 4 else None

    fig, (ax_e, ax_w) = plt.subplots(1, 2, figsize=(14, 5))

    mask_e = (energy >= args.energy_min) & (energy <= args.energy_max) & np.isfinite(energy)
    ax_e.plot(energy[mask_e], intensity[mask_e], lw=1.5, color="C0", label="Intensity")
    if norm is not None:
        ax_e.plot(energy[mask_e], norm[mask_e], lw=1.0, color="C1", ls="--", label="Normalized")
    ax_e.set_xlabel("Energy (eV)")
    ax_e.set_ylabel("Intensity (arb. units)")
    ax_e.set_title(f"Absorption: {path.name}")
    ax_e.set_xlim(args.energy_min, args.energy_max)
    ax_e.grid(True, alpha=0.3)
    ax_e.legend()

    mask_w = (wavelength >= args.wavelength_min) & (wavelength <= args.wavelength_max) & np.isfinite(wavelength)
    ax_w.plot(wavelength[mask_w], intensity[mask_w], lw=1.5, color="C3", label="Intensity")
    if norm is not None:
        ax_w.plot(wavelength[mask_w], norm[mask_w], lw=1.0, color="C4", ls="--", label="Normalized")
    ax_w.set_xlabel("Wavelength (nm)")
    ax_w.set_ylabel("Intensity (arb. units)")
    ax_w.set_xlim(args.wavelength_max, args.wavelength_min)
    ax_w.set_title("Wavelength domain")
    ax_w.grid(True, alpha=0.3)
    ax_w.legend()

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
