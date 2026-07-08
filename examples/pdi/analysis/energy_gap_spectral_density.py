"""Energy-gap spectral densities from TD-DFTB excitation trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

HBAR_EV_FS = 0.6582119569
KB_EV_K = 8.617333262145e-5
EV_TO_WAVENUMBER = 8065.543937


@dataclass(frozen=True)
class ExcitationStateSeries:
    state: int
    frame: np.ndarray
    time_fs: np.ndarray
    energy_ev: np.ndarray
    oscillator_strength: np.ndarray


@dataclass(frozen=True)
class SpectralDensityResult:
    state: int
    temperature: float
    sample_count: int
    dt_fs: float
    mean_energy_ev: float
    mean_oscillator_strength: float
    std_energy_ev: float
    omega_ev: np.ndarray
    wavenumber_cm: np.ndarray
    J_ev: np.ndarray


def _read_columns(path: Path) -> tuple[list[str], np.ndarray]:
    header = None
    with Path(path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("# Frame"):
                header = stripped[1:].split()
                break
    if header is None:
        raise ValueError(f"excitation column header not found in {path}")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(header):
        raise ValueError(
            f"data has {data.shape[1]} columns but header has {len(header)}"
        )
    return header, data


def read_excitation_state(path: Path, state: int) -> ExcitationStateSeries:
    """Read one one-based excited-state trajectory from horizontal output."""
    if state < 1:
        raise ValueError("state index must be at least 1")
    columns, data = _read_columns(Path(path))
    energy_label = f"E{state}(eV)"
    strength_patterns = (f"f{state}", f"Osc{state}")
    if energy_label not in columns:
        raise ValueError(f"state {state} energy column {energy_label} not found")
    strength_label = next((x for x in strength_patterns if x in columns), None)
    if strength_label is None:
        raise ValueError(f"state {state} oscillator-strength column not found")
    required = ["Frame", "Time(fs)", energy_label, strength_label]
    indices = [columns.index(label) for label in required]
    selected = data[:, indices]
    if selected.shape[0] < 4:
        raise ValueError("at least four excitation samples are required")
    if not np.all(np.isfinite(selected)):
        raise ValueError("excitation data contain nonfinite values")
    time_fs = selected[:, 1]
    differences = np.diff(time_fs)
    if np.any(differences <= 0.0):
        raise ValueError("time samples must be strictly increasing")
    if not np.allclose(differences, differences[0], rtol=1e-6, atol=1e-10):
        raise ValueError("time samples must be uniformly spaced")
    return ExcitationStateSeries(
        state=state,
        frame=selected[:, 0].astype(int),
        time_fs=time_fs,
        energy_ev=selected[:, 2],
        oscillator_strength=selected[:, 3],
    )


def build_spectral_density(
    series: ExcitationStateSeries,
    temperature: float = 300.0,
    n_fft: int | None = None,
) -> SpectralDensityResult:
    """Build J(omega)=beta*omega*S_deltaE(omega)/2 in eV units."""
    if temperature <= 0.0:
        raise ValueError("temperature must be positive")
    energy = np.asarray(series.energy_ev, dtype=float)
    if energy.size < 4 or not np.all(np.isfinite(energy)):
        raise ValueError("at least four finite energy samples are required")
    dt_fs = float(np.diff(series.time_fs)[0])
    if dt_fs <= 0.0 or not np.allclose(np.diff(series.time_fs), dt_fs):
        raise ValueError("time samples must be uniformly spaced")
    n = energy.size
    if n_fft is None:
        n_fft = n
    if n_fft < n:
        raise ValueError("n_fft must be at least the sample count")
    delta = energy - np.mean(energy)
    window = np.hanning(n)
    window_power = float(np.mean(window ** 2))
    transformed = np.fft.rfft(delta * window, n=n_fft)
    dt_natural = dt_fs / HBAR_EV_FS
    psd_two_sided = (
        dt_natural / (n * window_power) * np.abs(transformed) ** 2
    )
    omega_ev = 2.0 * np.pi * np.fft.rfftfreq(n_fft, d=dt_natural)
    beta = 1.0 / (KB_EV_K * temperature)
    J_ev = 0.5 * beta * omega_ev * psd_two_sided
    J_ev[0] = 0.0
    return SpectralDensityResult(
        state=series.state,
        temperature=float(temperature),
        sample_count=n,
        dt_fs=dt_fs,
        mean_energy_ev=float(np.mean(energy)),
        mean_oscillator_strength=float(np.mean(series.oscillator_strength)),
        std_energy_ev=float(np.std(energy)),
        omega_ev=omega_ev,
        wavenumber_cm=omega_ev * EV_TO_WAVENUMBER,
        J_ev=J_ev,
    )


def write_spectral_density(path: Path, result: SpectralDensityResult) -> None:
    """Write spectral-density columns and scalar metadata."""
    metadata = [
        f"state={result.state}",
        f"temperature={result.temperature:.16g}",
        f"sample_count={result.sample_count}",
        f"dt_fs={result.dt_fs:.16g}",
        f"mean_energy_ev={result.mean_energy_ev:.16g}",
        f"mean_oscillator_strength={result.mean_oscillator_strength:.16g}",
        f"std_energy_ev={result.std_energy_ev:.16g}",
        "columns=omega_eV wavenumber_cm-1 J_eV",
    ]
    np.savetxt(
        path,
        np.column_stack(
            [result.omega_ev, result.wavenumber_cm, result.J_ev]
        ),
        header="\n".join(metadata),
        comments="# ",
        fmt="%.12e",
    )


def read_spectral_density(path: Path) -> SpectralDensityResult:
    """Read a spectral-density file written by :func:`write_spectral_density`."""
    metadata: dict[str, str] = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if "=" in text:
                key, value = text.split("=", 1)
                metadata[key.strip()] = value.strip()
    required = {
        "state", "temperature", "sample_count", "dt_fs",
        "mean_energy_ev", "mean_oscillator_strength", "std_energy_ev",
    }
    missing = required - metadata.keys()
    if missing:
        raise ValueError(f"missing spectral-density metadata: {sorted(missing)}")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != 3 or data.shape[0] < 2 or not np.all(np.isfinite(data)):
        raise ValueError("spectral-density file must contain three finite columns")
    return SpectralDensityResult(
        state=int(metadata["state"]),
        temperature=float(metadata["temperature"]),
        sample_count=int(metadata["sample_count"]),
        dt_fs=float(metadata["dt_fs"]),
        mean_energy_ev=float(metadata["mean_energy_ev"]),
        mean_oscillator_strength=float(metadata["mean_oscillator_strength"]),
        std_energy_ev=float(metadata["std_energy_ev"]),
        omega_ev=data[:, 0],
        wavenumber_cm=data[:, 1],
        J_ev=data[:, 2],
    )
