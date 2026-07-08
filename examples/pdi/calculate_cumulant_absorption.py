"""CLI for second-order-cumulant absorption from a spectral density."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .absorption import absorption_spectrum, oscillator_strength_to_mu2
from .energy_gap_spectral_density import HBAR_EV_FS, KB_EV_K, read_spectral_density

NM_PER_EV = 1239.841984


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectral_density", type=Path)
    parser.add_argument("--temperature", type=float)
    parser.add_argument("--t-max-fs", type=float, default=500.0)
    parser.add_argument("--dt-fs", type=float, default=0.1)
    parser.add_argument("--damping-ev", type=float, default=0.0)
    parser.add_argument("--energy-min", type=float, default=0.0)
    parser.add_argument("--energy-max", type=float, default=4.0)
    parser.add_argument("--n-energy", type=int, default=2001)
    parser.add_argument("--output", type=Path, default=Path("absorption.dat"))
    parser.add_argument("--correlation-output", type=Path)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    sd = read_spectral_density(args.spectral_density)
    temperature = sd.temperature if args.temperature is None else args.temperature
    if temperature <= 0.0:
        raise ValueError("temperature must be positive")
    if args.energy_min < 0.0 or args.energy_max <= args.energy_min:
        raise ValueError("invalid output energy range")
    if args.n_energy < 2:
        raise ValueError("n-energy must be at least 2")
    if args.damping_ev < 0.0:
        raise ValueError("damping must be nonnegative")
    mu2 = oscillator_strength_to_mu2(
        sd.mean_oscillator_strength, sd.mean_energy_ev
    )
    beta = 1.0 / (KB_EV_K * temperature)
    omega_out = np.linspace(args.energy_min, args.energy_max, args.n_energy)
    omega, intensity, chi = absorption_spectrum(
        [mu2], [sd.mean_energy_ev], sd.omega_ev, sd.J_ev, beta,
        args.t_max_fs / HBAR_EV_FS, args.dt_fs / HBAR_EV_FS,
        damping=args.damping_ev, omega_out=omega_out,
    )
    positive_max = float(np.max(np.maximum(intensity, 0.0)))
    normalized = intensity / positive_max if positive_max > 0.0 else np.zeros_like(intensity)
    wavelength = np.full_like(omega, np.inf)
    np.divide(NM_PER_EV, omega, out=wavelength, where=omega > 0.0)
    np.savetxt(
        args.output,
        np.column_stack([omega, wavelength, intensity, normalized]),
        header=(
            f"state={sd.state}\n"
            f"temperature={temperature:.16g}\n"
            f"mean_energy_ev={sd.mean_energy_ev:.16g}\n"
            f"mean_oscillator_strength={sd.mean_oscillator_strength:.16g}\n"
            "columns=energy_eV wavelength_nm intensity normalized_intensity"
        ),
        comments="# ", fmt="%.12e",
    )
    if args.correlation_output is not None:
        time_natural = np.arange(chi.size) * (args.dt_fs / HBAR_EV_FS)
        np.savetxt(
            args.correlation_output,
            np.column_stack([
                time_natural * HBAR_EV_FS, np.real(chi), np.imag(chi), np.abs(chi)
            ]),
            header="columns=time_fs chi_real chi_imag chi_abs",
            comments="# ", fmt="%.12e",
        )
    print(f"Saved S{sd.state} cumulant absorption spectrum to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
