"""CLI for building a spectral density from an excitation trajectory."""

from __future__ import annotations

import argparse
from pathlib import Path

from .energy_gap_spectral_density import (
    build_spectral_density,
    read_excitation_state,
    write_spectral_density,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--state", type=int, default=5)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--n-fft", type=int)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    output = args.output or Path(f"spectral_density_S{args.state}.dat")
    series = read_excitation_state(args.input, args.state)
    result = build_spectral_density(series, args.temperature, args.n_fft)
    write_spectral_density(output, result)
    print(
        f"Saved S{result.state} spectral density to {output} "
        f"({result.sample_count} samples, mean energy {result.mean_energy_ev:.6f} eV)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
