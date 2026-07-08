import numpy as np
import pytest

from q4mdkit.analysis.energy_gap_spectral_density import (
    HBAR_EV_FS,
    ExcitationStateSeries,
    build_spectral_density,
    read_spectral_density,
    read_excitation_state,
    write_spectral_density,
)
from q4mdkit.analysis.build_excitation_spectral_density import main as build_main


def test_reads_selected_state_from_horizontal_excitation_file(tmp_path):
    path = tmp_path / "excitations.dat"
    path.write_text(
        "# TD-DFTB Excitation Data\n"
        "# Frame  Time(fs)  E1(eV)  f1  E2(eV)  f2\n"
        "0  0.0  1.0  0.1  2.0  0.2\n"
        "1  4.0  1.1  0.1  2.2  0.3\n"
        "2  8.0  0.9  0.1  1.8  0.4\n"
        "3 12.0  1.0  0.1  2.0  0.5\n"
    )

    series = read_excitation_state(path, state=2)

    assert series.state == 2
    assert np.allclose(series.energy_ev, [2.0, 2.2, 1.8, 2.0])
    assert np.allclose(series.oscillator_strength, [0.2, 0.3, 0.4, 0.5])


def test_rejects_unavailable_state(tmp_path):
    path = tmp_path / "excitations.dat"
    path.write_text(
        "# Frame  Time(fs)  E1(eV)  f1\n"
        "0 0.0 1.0 0.1\n1 4.0 1.1 0.1\n2 8.0 0.9 0.1\n3 12.0 1.0 0.1\n"
    )
    with pytest.raises(ValueError, match="state 2"):
        read_excitation_state(path, state=2)


def test_sinusoidal_gap_produces_nonnegative_peak_at_expected_frequency():
    n = 128
    dt_fs = 1.0
    k = 8
    time_fs = np.arange(n) * dt_fs
    omega_expected = 2.0 * np.pi * k * HBAR_EV_FS / (n * dt_fs)
    energy = 2.0 + 0.05 * np.cos(omega_expected * time_fs / HBAR_EV_FS)
    series = ExcitationStateSeries(
        state=3,
        frame=np.arange(n),
        time_fs=time_fs,
        energy_ev=energy,
        oscillator_strength=np.full(n, 0.4),
    )

    result = build_spectral_density(series, temperature=300.0)

    assert result.J_ev[0] == 0.0
    assert np.all(result.J_ev >= 0.0)
    assert np.isclose(result.omega_ev[np.argmax(result.J_ev)], omega_expected)


def test_spectral_density_metadata_round_trips(tmp_path):
    time_fs = np.arange(8, dtype=float) * 2.0
    series = ExcitationStateSeries(
        state=2,
        frame=np.arange(8),
        time_fs=time_fs,
        energy_ev=2.0 + 0.01 * np.cos(np.arange(8)),
        oscillator_strength=np.full(8, 0.25),
    )
    result = build_spectral_density(series, temperature=250.0)
    path = tmp_path / "sd.dat"

    write_spectral_density(path, result)
    loaded = read_spectral_density(path)

    assert loaded.state == 2
    assert loaded.temperature == 250.0
    assert np.isclose(loaded.mean_oscillator_strength, 0.25)
    assert np.allclose(loaded.J_ev, result.J_ev)


def test_spectral_density_cli_selects_requested_state(tmp_path):
    source = tmp_path / "excitations.dat"
    source.write_text(
        "# Frame Time(fs) E1(eV) f1 E2(eV) f2\n"
        "0 0 1.0 .1 2.0 .2\n1 1 1.1 .1 2.1 .2\n"
        "2 2 .9 .1 1.9 .2\n3 3 1.0 .1 2.0 .2\n"
    )
    output = tmp_path / "state2.dat"

    assert build_main([str(source), "--state", "2", "--output", str(output)]) == 0
    assert read_spectral_density(output).state == 2
