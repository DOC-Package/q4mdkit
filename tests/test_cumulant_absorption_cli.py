import numpy as np

import q4mdkit.analysis.absorption as absorption_module
from q4mdkit.analysis.absorption import (
    dipole_correlation,
    lineshape_function,
    oscillator_strength_to_mu2,
)
from q4mdkit.analysis.calculate_cumulant_absorption import main as absorption_main
from q4mdkit.analysis.energy_gap_spectral_density import (
    SpectralDensityResult,
    write_spectral_density,
)


def test_oscillator_strength_converts_to_squared_transition_dipole():
    expected = 3.0 * 0.5 / (2.0 * (2.0 / 27.211386245988))
    assert np.isclose(oscillator_strength_to_mu2(0.5, 2.0), expected)


def test_zero_spectral_density_gives_undamped_phase_correlation():
    t = np.linspace(0.0, 5.0, 11)
    omega = np.linspace(0.0, 1.0, 32)
    chi = dipole_correlation(t, [2.0], [0.4], omega, np.zeros_like(omega), 10.0)
    assert np.allclose(chi, 2.0 * np.exp(-1j * 0.4 * t))


def test_lineshape_bounds_temporary_kernel_size(monkeypatch):
    t = np.linspace(0.0, 5.0, 130)
    omega = np.linspace(0.0, 1.0, 301)
    spectral_density = 0.002 * omega * np.exp(-omega / 0.08)
    expected = lineshape_function(t, omega, spectral_density, beta=10.0)
    original_outer = np.outer

    def bounded_outer(left, right):
        assert np.size(left) <= 64
        return original_outer(left, right)

    monkeypatch.setattr(absorption_module.np, "outer", bounded_outer)
    actual = lineshape_function(t, omega, spectral_density, beta=10.0)

    assert np.allclose(actual, expected)


def test_absorption_cli_writes_normalized_spectrum(tmp_path):
    omega = np.linspace(0.0, 0.4, 65)
    sd = SpectralDensityResult(
        state=5, temperature=300.0, sample_count=100, dt_fs=4.0,
        mean_energy_ev=2.2, mean_oscillator_strength=0.5,
        std_energy_ev=0.05, omega_ev=omega,
        wavenumber_cm=omega * 8065.543937,
        J_ev=0.002 * omega * np.exp(-omega / 0.08),
    )
    source = tmp_path / "sd.dat"
    output = tmp_path / "absorption.dat"
    write_spectral_density(source, sd)

    result = absorption_main([
        str(source), "--output", str(output), "--t-max-fs", "20",
        "--dt-fs", "0.2", "--energy-min", "1.5",
        "--energy-max", "3.0", "--n-energy", "101",
    ])

    data = np.loadtxt(output, comments="#")
    assert result == 0
    assert data.shape == (101, 4)
    assert np.all(np.isfinite(data))
    assert np.isclose(np.max(data[:, 3]), 1.0)
