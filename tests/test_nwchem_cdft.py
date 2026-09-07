from pathlib import Path
import subprocess

import numpy as np
import pytest

import q4mdkit.analysis.nwchem_cdft as nwchem_cdft
from q4mdkit.analysis.nwchem_cdft import (
    NWChemCDFTConfig,
    NWChemConstraint,
    build_cdft_directive,
    load_config,
    parse_atom_range,
    parse_nwchem_energy,
    parse_nwchem_populations,
    render_nwchem_input,
    run_nwchem_cdft_analysis,
)


def test_parse_atom_range_converts_nwchem_1_based_range_to_python_indices():
    assert parse_atom_range("1:36") == [0, 36]
    assert parse_atom_range("37:72") == [36, 72]
    assert parse_atom_range("1:1") == [0, 1]


def test_build_cdft_directive_supports_charge_and_spin_constraints():
    assert (
        build_cdft_directive(NWChemConstraint(kind="charge", atoms="1:36", value=0.5))
        == "cdft 1 36 charge 0.5 pop lowdin"
    )
    assert (
        build_cdft_directive(
            NWChemConstraint(kind="spin", atoms="1:36", value=-1.0, population="becke")
        )
        == "cdft 1 36 spin -1.0 pop becke"
    )


def test_render_nwchem_input_contains_geometry_cdft_and_point_charges():
    constraints = [
        NWChemConstraint(kind="charge", atoms="1:2", value=1.0, population="lowdin"),
        NWChemConstraint(kind="spin", atoms="3:4", value=-1.0, population="becke"),
    ]
    point_charges = np.array([[0.0, 0.0, 3.0, -0.2], [0.0, 0.0, -3.0, 0.2]])

    text = render_nwchem_input(
        atom_types=["C", "C", "O", "O"],
        coords_ang=np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.0, 2.0], [1.2, 0.0, 2.0]]
        ),
        constraints=constraints,
        basis="6-31G*",
        xc="b3lyp",
        charge=0,
        multiplicity=1,
        point_charges=point_charges,
    )

    assert "geometry noautoz nocenter noautosym" in text
    assert "C     0.0000000000    0.0000000000    0.0000000000" in text
    assert "bq     0.0000000000    0.0000000000    3.0000000000   -0.2000000000" in text
    assert "xc b3lyp" in text
    assert "mult 1" in text
    assert "cdft 1 2 charge 1.0 pop lowdin" in text
    assert "cdft 3 4 spin -1.0 pop becke" in text
    assert "task dft energy" in text


def test_render_nwchem_input_can_print_full_mo_vectors():
    text = render_nwchem_input(
        atom_types=["C", "H"],
        coords_ang=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        constraints=[NWChemConstraint(kind="charge", atoms="1:1", value=0.0)],
        print_mo_vectors=True,
    )

    assert "set movecs:tanalyze 0.0" in text
    assert 'print "final vectors analysis"' in text


def test_load_config_reads_trajectory_and_constraints(tmp_path):
    config_path = tmp_path / "nwchem_cdft.yaml"
    config_path.write_text(
        """
input:
  trajectory: traj/nve.dcd
  topology: traj/top.pdb
  qm_atoms: traj/qmatoms
  basis: 6-31G*
  xc: b3lyp
  charge: 0
  multiplicity: 1
output:
  directory: output_nwchem
  energy_file: energies.dat
  population_file: populations.dat
frames:
  start: 10
  n_frames: 5
  t0_fs: 20.0
  dt_fs: 2.0
constraints:
  - kind: charge
    atoms: "1:36"
    value: 1.0
    population: lowdin
  - kind: spin
    atoms: "1:36"
    value: 1.0
    population: becke
nwchem:
  binary: /opt/nwchem/bin/nwchem
  timeout: 120
""",
        encoding="utf-8",
    )

    config = load_config(config_path)

    assert isinstance(config, NWChemCDFTConfig)
    assert config.traj_path == tmp_path / "traj/nve.dcd"
    assert config.topology_path == tmp_path / "traj/top.pdb"
    assert config.qm_atoms_file == tmp_path / "traj/qmatoms"
    assert config.output_dir == tmp_path / "output_nwchem"
    assert config.energy_file == tmp_path / "output_nwchem/energies.dat"
    assert config.population_file == tmp_path / "output_nwchem/populations.dat"
    assert config.start_frame == 10
    assert config.n_frames == 5
    assert config.t0_fs == 20.0
    assert config.dt_fs == 2.0
    assert config.nwchem_binary == "/opt/nwchem/bin/nwchem"
    assert config.timeout == 120
    assert [c.kind for c in config.constraints] == ["charge", "spin"]
    assert [c.value for c in config.constraints] == [1.0, 1.0]


def test_parse_nwchem_energy_extracts_total_dft_energy():
    output = """
     Total DFT energy =       -76.439812345678
     One electron energy = ...
    """
    assert parse_nwchem_energy(output) == pytest.approx(-76.439812345678)


def test_parse_nwchem_energy_rejects_output_without_energy():
    with pytest.raises(ValueError, match="Total DFT energy"):
        parse_nwchem_energy("Total SCF energy = -1.0")


def test_parse_nwchem_populations_reads_common_labels():
    assert parse_nwchem_populations("constraint value = 1.0") == [1.0]
    assert parse_nwchem_populations("CDFT constrained value = -1.0") == [-1.0]
    assert parse_nwchem_populations("Total DFT energy = -1.0") is None


def test_parse_nwchem_populations_reads_final_multipliers():
    output = """
  CDFT multipliers:
      1        0.0842571009
      2       -0.0558976515
      iter =   14
  CDFT final multipliers
      1        0.1780106061
      2       -0.0309081216


         Total DFT energy =    -1684.099692714465
"""
    assert parse_nwchem_populations(output) == pytest.approx(
        [0.1780106061, -0.0309081216]
    )


def _write_runner_config(tmp_path: Path) -> Path:
    (tmp_path / "traj.dcd").write_text("", encoding="utf-8")
    (tmp_path / "top.pdb").write_text("", encoding="utf-8")
    (tmp_path / "qmatoms").write_text("1 2 3 4\n", encoding="utf-8")
    config_path = tmp_path / "nwchem_cdft.yaml"
    config_path.write_text(
        """
input:
  trajectory: traj.dcd
  topology: top.pdb
  qm_atoms: qmatoms
  basis: 6-31G*
  xc: b3lyp
  charge: 0
  multiplicity: 1
output:
  directory: out
frames:
  start: 1
  n_frames: 2
  t0_fs: 10.0
  dt_fs: 2.0
constraints:
  - kind: charge
    atoms: "1:2"
    value: 1.0
    population: lowdin
  - kind: spin
    atoms: "3:4"
    value: -1.0
    population: becke
""",
        encoding="utf-8",
    )
    return config_path


def _qm_frame(coords):
    return np.asarray(coords, dtype=float)


def test_run_nwchem_cdft_analysis_writes_energy_and_population(monkeypatch, tmp_path):
    config_path = _write_runner_config(tmp_path)

    monkeypatch.setattr(
        nwchem_cdft,
        "_extract_atom_types_from_topology",
        lambda topology_path, qm_indices: ["C", "C", "O", "O"],
    )
    frames = [
        (0, None, _qm_frame([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.0, 2.0], [1.2, 0.0, 2.0]])),
        (1, None, _qm_frame([[0.1, 0.0, 0.0], [1.3, 0.0, 0.0], [0.0, 0.0, 2.1], [1.3, 0.0, 2.1]])),
        (2, None, _qm_frame([[0.2, 0.0, 0.0], [1.4, 0.0, 0.0], [0.0, 0.0, 2.2], [1.4, 0.0, 2.2]])),
    ]
    monkeypatch.setattr(nwchem_cdft, "_iter_qm_coordinates", lambda *a, **k: iter(frames))

    def fake_run(work_dir, config):
        return subprocess.CompletedProcess(
            args=[config.nwchem_binary, "nwchem.inp"],
            returncode=0,
            stdout=(
                "Total DFT energy = -76.5000000000\n"
                "constraint value = 1.0000000000\n"
                "constraint value = -1.0000000000\n"
            ),
            stderr="",
        )

    monkeypatch.setattr(nwchem_cdft, "_run_nwchem", fake_run)

    run_nwchem_cdft_analysis(config_path)

    energy_lines = (tmp_path / "out" / "energies.dat").read_text().splitlines()
    population_lines = (tmp_path / "out" / "populations.dat").read_text().splitlines()

    assert len(energy_lines) == 3  # header + 2 frames
    assert "12.000" in energy_lines[1]
    assert "-76.5000000000" in energy_lines[1]
    assert "14.000" in energy_lines[2]

    assert len(population_lines) == 3
    assert "1.00000000" in population_lines[1]
    assert "-1.00000000" in population_lines[1]
    assert "1.00000000" in population_lines[2]
    assert "-1.00000000" in population_lines[2]


def test_run_nwchem_cdft_analysis_writes_nan_on_failure(monkeypatch, tmp_path):
    config_path = _write_runner_config(tmp_path)

    monkeypatch.setattr(
        nwchem_cdft,
        "_extract_atom_types_from_topology",
        lambda topology_path, qm_indices: ["C", "C", "O", "O"],
    )
    monkeypatch.setattr(
        nwchem_cdft,
        "_iter_qm_coordinates",
        lambda *a, **k: iter([(1, None, _qm_frame([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.0, 2.0], [1.2, 0.0, 2.0]]))]),
    )

    def fake_run(work_dir, config):
        return subprocess.CompletedProcess(
            args=[config.nwchem_binary, "nwchem.inp"],
            returncode=1,
            stdout="",
            stderr="boom",
        )

    monkeypatch.setattr(nwchem_cdft, "_run_nwchem", fake_run)

    run_nwchem_cdft_analysis(config_path)

    energy_lines = (tmp_path / "out" / "energies.dat").read_text().splitlines()
    population_lines = (tmp_path / "out" / "populations.dat").read_text().splitlines()

    assert len(energy_lines) == 2  # header + 1 frame
    assert "nan" in energy_lines[1]
    assert len(population_lines) == 2
    assert population_lines[1].count("nan") == 2


def test_run_nwchem_passes_binary_timeout_and_scratch_dir(monkeypatch, tmp_path):
    calls = {}

    def fake_run(cmd, **kwargs):
        calls["cmd"] = cmd
        calls["kwargs"] = kwargs
        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(nwchem_cdft.subprocess, "run", fake_run)

    work_dir = tmp_path / "work"
    work_dir.mkdir()
    scratch = tmp_path / "scratch"
    config = NWChemCDFTConfig(
        traj_path=tmp_path / "t",
        topology_path=tmp_path / "p",
        qm_atoms_file=tmp_path / "q",
        nwchem_binary="/fake/nwchem",
        scratch_dir=scratch,
        timeout=12,
    )

    nwchem_cdft._run_nwchem(work_dir, config)

    assert calls["cmd"] == ["/fake/nwchem", "nwchem.inp", str(scratch)]
    assert calls["kwargs"]["cwd"] == work_dir
    assert calls["kwargs"]["timeout"] == 12
