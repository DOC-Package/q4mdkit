"""ORCA snapshot analysis, with only the external QM executable replaced."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import yaml

from q4mdkit.analysis import tdorca


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/vg30/analysis"
SPECTRUM = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"


@pytest.mark.parametrize("table", [
    "1 16131.088 619.92095 0.1 0.2 0.3\n2 24196.632 413.28063 0.2 0.3 0.4\n",
    "0-1A -> 1-1A 2.0 16131.088 619.92095 0.1 0.2 0.3\n"
    "0-1A -> 2-1A 3.0 24196.632 413.28063 0.2 0.3 0.4\n",
])
def test_reads_only_electric_dipole_absorption_table(tmp_path, table):
    output = tmp_path / "orca.out"
    output.write_text(
        SPECTRUM + "\n-----\nState/Energy/Wavelength/fosc\n-----\n" + table
        + "\nABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\n"
        + table.replace("0.1", "0.9")
    )
    states = tdorca.parse_orca_output(output)
    assert [state.energy_ev for state in states] == pytest.approx([2.0, 3.0])
    assert [state.oscillator_strength for state in states] == pytest.approx([0.1, 0.2])
    assert states[0].wavelength_nm == pytest.approx(619.92095)


def test_xyzfile_template_preserves_charge_and_multiplicity(tmp_path):
    template = tmp_path / "template.inp"
    template.write_text("! B3LYP def2-SVP\n* xyzfile 1 1 original.xyz\n")
    geometry = tmp_path / "qm.xyz"
    geometry.write_text("1\n\nH 1.0 2.0 3.0\n")
    target = tmp_path / "orca.inp"
    tdorca.prepare_orca_input(template, geometry, target, n_roots=2)
    assert f"* xyzfile 1 1 {geometry}" in target.read_text()


def test_single_mm_point_charge_template(tmp_path):
    path = tmp_path / "PCcharges.dat"
    path.write_text("99.0 98.0 97.0 -0.4\n")
    assert tdorca.load_pccharges(path).tolist() == [-0.4]


@pytest.fixture
def snapshot_run(tmp_path, monkeypatch):
    md = pytest.importorskip("mdtraj")
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("VG3", chain)
    topology.add_atom("C", md.element.carbon, residue)
    topology.add_atom("H", md.element.hydrogen, residue)
    solvent = topology.add_residue("ACE", chain)
    topology.add_atom("C", md.element.carbon, solvent)
    topology.add_atom("N", md.element.nitrogen, solvent)
    positions = np.array([[0, 0.1, 0.2], [0.1, 0.2, 0.3],
                          [0.2, 0.3, 0.4], [0.3, 0.4, 0.5]])
    xyz = np.stack([positions + frame * 0.1 for frame in range(4)])
    trajectory = md.Trajectory(xyz, topology)
    trajectory.save_dcd(str(tmp_path / "trajectory.dcd"))
    trajectory[0].save_pdb(str(tmp_path / "topology.pdb"))
    (tmp_path / "qmatoms").write_text("0 1\n")
    (tmp_path / "PCcharges.dat").write_text("999 999 999 -0.2\n888 888 888 0.2\n")
    binary = tmp_path / "fake-orca"
    binary.write_text(
        f"#!{sys.executable}\n" + r'''
import json, os, pathlib, re, sys
text = pathlib.Path(sys.argv[1]).read_text()
assert '* xyz 1 1' in text
assert re.search(r'nroots\s+2', text, re.I)
assert re.search(r'tda\s+false', text, re.I)
coords = text.split('* xyz 1 1')[1].split('*')[0].strip().splitlines()
x = float(coords[0].split()[1])
saved = pathlib.Path('../excitations.dat')
record = {'x': x, 'pc': None,
          'saved_rows': sum(bool(line.strip()) and not line.startswith('#')
                            for line in saved.read_text().splitlines()) if saved.exists() else 0}
match = re.search(r'%pointcharges\s+"([^"]+)"', text, re.I)
if match:
    lines = pathlib.Path(match[1]).read_text().splitlines()
    record['pc'] = [[float(v) for v in line.split()] for line in lines[1:]]
with open(os.environ['ORCA_CALL_LOG'], 'a') as handle:
    handle.write(json.dumps(record) + '\n')
print('ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS')
print('Transition Energy(eV) Energy(cm-1) Wavelength(nm) fosc(D2)')
energies = [2.0 + 0.1*x, 3.0 + 0.1*x]
if os.environ.get('SHORT_X') and abs(x-float(os.environ['SHORT_X'])) < 1e-6:
    energies = energies[:1]
for state, energy in enumerate(energies, 1):
    print(f'0-1A -> {state}-1A {energy:.6f} {energy*8065.544:.3f} '
          f'{1239.8419/energy:.5f} {state*0.1:.6f} 0 0 0 0')
print('')
print('****ORCA TERMINATED NORMALLY****')
if os.environ.get('FAIL_ALL') == '1' or (
    os.environ.get('FAIL_X') and abs(x-float(os.environ['FAIL_X'])) < 1e-6
):
    sys.exit(7)
'''
    )
    binary.chmod(0o755)
    call_log = tmp_path / "calls.jsonl"
    monkeypatch.setenv("ORCA_CALL_LOG", str(call_log))
    config = {
        "input": {"trajectory": "trajectory.dcd", "topology": "topology.pdb",
                  "qm_atoms": "qmatoms",
                  "orca_template": str(EXAMPLE / "input/orca_tddft.inp")},
        "output": {"directory": "results", "excitations_format": "horizontal",
                   "average_spectrum": False},
        "frames": {"start": 1, "n_frames": 2, "t0_fs": 10.0, "dt_fs": 4.0},
        "orca": {"binary": str(binary), "num_procs": 2, "memory_mb": 500},
        "tddft": {"n_roots": 2},
    }
    return tmp_path, config, call_log


def run_snapshot_fixture(fixture):
    directory, config, _ = fixture
    path = directory / "settings.yaml"
    path.write_text(yaml.safe_dump(config))
    return subprocess.run(
        [sys.executable, str(EXAMPLE / "run_tdorca.py"), str(path)],
        cwd=ROOT, capture_output=True, text=True,
    )


@pytest.mark.parametrize("point_charges", [False, True])
def test_snapshot_entrypoint_matches_tddftb_frame_units_and_output(snapshot_run, point_charges):
    directory, config, call_log = snapshot_run
    if point_charges:
        config["input"]["pccharges_template"] = "PCcharges.dat"
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode == 0, result.stdout + result.stderr
    records = [json.loads(line) for line in call_log.read_text().splitlines()]
    assert [record["x"] for record in records] == pytest.approx([1.0, 2.0])
    assert [record["saved_rows"] for record in records] == [0, 1]
    if point_charges:
        np.testing.assert_allclose(records[0]["pc"],
                                   [[-0.2, 3.0, 4.0, 5.0], [0.2, 4.0, 5.0, 6.0]], atol=1e-6)
        np.testing.assert_allclose(records[1]["pc"],
                                   [[-0.2, 4.0, 5.0, 6.0], [0.2, 5.0, 6.0, 7.0]], atol=1e-6)
    else:
        assert all(record["pc"] is None for record in records)
    data = np.loadtxt(directory / "results/excitations.dat")
    np.testing.assert_allclose(data, [[1, 14, 2.1, 0.1, 3.1, 0.2],
                                    [2, 18, 2.2, 0.1, 3.2, 0.2]])
    assert not (directory / "results/absorption_spectrum.dat").exists()


def test_failed_snapshot_does_not_append_results(snapshot_run, monkeypatch):
    directory, _, _ = snapshot_run
    monkeypatch.setenv("FAIL_X", "2.0")
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode == 0, result.stdout + result.stderr
    data = np.loadtxt(directory / "results/excitations.dat")
    assert data.tolist() == [1, 14, 2.1, 0.1, 3.1, 0.2]
    assert "FAILED" in result.stdout


def test_all_failed_snapshots_return_failure(snapshot_run, monkeypatch):
    monkeypatch.setenv("FAIL_ALL", "1")
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode != 0
    assert "No successful ORCA" in result.stderr


@pytest.mark.parametrize("requested_roots", [1, 2, 3])
def test_legacy_long_format_and_spectrum_are_preserved(snapshot_run, requested_roots):
    directory, config, _ = snapshot_run
    config["output"] = {"directory": "results"}
    # Legacy complete templates can request a different count from the YAML.
    template = directory / "fixed_roots.inp"
    template.write_text((EXAMPLE / "input/orca_tddft.inp").read_text().replace("{NROOTS}", "2"))
    config["input"]["orca_template"] = str(template)
    config["tddft"]["n_roots"] = requested_roots
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode == 0, result.stdout + result.stderr
    data = np.loadtxt(directory / "results/excitations.dat")
    assert data.shape == (4, 6)
    np.testing.assert_allclose(data[:, 2], [1, 2, 1, 2])
    assert (directory / "results/absorption_spectrum.dat").is_file()


def test_incomplete_snapshot_does_not_write_a_ragged_horizontal_row(snapshot_run, monkeypatch):
    directory, _, _ = snapshot_run
    monkeypatch.setenv("SHORT_X", "2.0")
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode == 0, result.stdout + result.stderr
    data = np.loadtxt(directory / "results/excitations.dat")
    assert data.tolist() == [1, 14, 2.1, 0.1, 3.1, 0.2]
    assert "INCOMPLETE EXCITATIONS" in result.stdout


def test_point_charge_count_is_checked_before_orca(snapshot_run):
    directory, config, call_log = snapshot_run
    (directory / "PCcharges.dat").write_text("0 0 0 -0.2\n")
    config["input"]["pccharges_template"] = "PCcharges.dat"
    result = run_snapshot_fixture(snapshot_run)
    assert result.returncode != 0
    assert "Point-charge count 1 does not match MM atom count 2" in result.stderr
    assert not call_log.exists()
