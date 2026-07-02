from dataclasses import fields

from q4mdkit.analysis.tddftb import (
    Excitation,
    TDDFTBConfig,
    append_excitation_frame,
    initialize_excitation_output,
)


def test_excitation_output_writes_energy_strength_pairs_horizontally(tmp_path):
    output_path = tmp_path / "excitations.dat"
    frame_data = {
        "frame": 5,
        "time_fs": 20.0,
        "excitations": [
            Excitation(energy_ev=1.1, oscillator_strength=0.2),
            Excitation(energy_ev=1.3, oscillator_strength=0.4),
        ],
    }

    initialize_excitation_output(output_path, n_excitations=2)
    append_excitation_frame(output_path, frame_data, n_excitations=2)

    assert output_path.read_text().splitlines() == [
        "# TD-DFTB Excitation Data",
        "# Frame  Time(fs)  E1(eV)  f1  E2(eV)  f2",
        "     5  20.00  1.1000  0.200000  1.3000  0.400000",
    ]


def test_excitation_output_limits_pairs_to_n_excitations(tmp_path):
    output_path = tmp_path / "excitations.dat"
    frame_data = {
        "frame": 0,
        "time_fs": 0.0,
        "excitations": [
            Excitation(energy_ev=1.0, oscillator_strength=0.1),
            Excitation(energy_ev=2.0, oscillator_strength=0.2),
            Excitation(energy_ev=3.0, oscillator_strength=0.3),
        ],
    }

    initialize_excitation_output(output_path, n_excitations=2)
    append_excitation_frame(output_path, frame_data, n_excitations=2)

    data_line = output_path.read_text().splitlines()[-1]
    assert data_line.split() == [
        "0",
        "0.00",
        "1.0000",
        "0.100000",
        "2.0000",
        "0.200000",
    ]


def test_tddftb_config_has_no_tracking_fields():
    field_names = {field.name for field in fields(TDDFTBConfig)}
    spectrum_fields = {
        "spectra_file",
        "output_average_spectrum",
        "output_individual_spectra",
        "broadening_ev",
        "energy_range_ev",
        "n_energy_points",
    }

    assert not any(
        "tracking" in name or name.startswith("track_")
        for name in field_names
    )
    assert "state_timeseries_file" not in field_names
    assert field_names.isdisjoint(spectrum_fields)
