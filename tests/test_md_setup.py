"""Shared MD preparation retains restart behavior without importing ASH."""

from pathlib import Path
import xml.etree.ElementTree as ET

import pytest
import yaml

from q4mdkit import md
from q4mdkit.mm.mm_config import MMConfig
from q4mdkit.qmmm.md_config import MDConfig


@pytest.fixture
def run_inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('system.pdb').write_text(
        'CRYST1   44.255   44.255   44.255  90.00  90.00  90.00 P 1           1\n'
        'ATOM      1  C   DYE A   1       1.000   2.000   3.000  1.00  0.00           C\n'
    )
    Path('system.prmtop').write_text('%FLAG POINTERS\n%FORMAT(10I8)\n       1\n')
    Path('mm.yaml').write_text(yaml.safe_dump({
        'paths': {'pdbfile': 'system.pdb', 'amber_prmtop': 'system.prmtop'},
        'openmm': {'periodic': True},
    }))
    Path('md.yaml').write_text(yaml.safe_dump({
        'md': {'timestep': 0.0005, 'simulation_time': 1.0, 'nvt': {'temperature': 300}},
        'output': {'directory': 'output', 'protect_existing': True},
    }))
    Path('state.xml').write_text(
        '<State><PeriodicBoxVectors>'
        '<A x="4.4255050604395985" y="0" z="0"/>'
        '<B x="0" y="4.4255050604395985" z="0"/>'
        '<C x="0" y="0" z="4.4255050604395985"/>'
        '</PeriodicBoxVectors><Parameters MonteCarloPressure="1" MonteCarloTemperature="300"/>'
        '<Positions><Position x="1" y="2" z="3"/></Positions>'
        '<Velocities><Velocity x="0.1" y="0.2" z="0.3"/></Velocities></State>'
    )
    return MMConfig('mm.yaml'), MDConfig('md.yaml')


@pytest.mark.parametrize('ensemble', ['nvt', 'nve'])
def test_prepares_constant_volume_restart_with_exact_box(run_inputs, ensemble):
    theory, settings = run_inputs
    original = Path('state.xml').read_bytes()
    statefile = md.prepare_md_run(theory, settings, ensemble, statefile='state.xml')
    converted = ET.parse(statefile).getroot()
    source = ET.fromstring(original)
    assert converted.find('Parameters') is None
    for name in ['Positions', 'Velocities', 'PeriodicBoxVectors']:
        assert ET.tostring(converted.find(name)) == ET.tostring(source.find(name))
    assert theory.pbc_vectors[0][0] == pytest.approx(44.255050604395985, abs=1e-12)
    assert Path('state.xml').read_bytes() == original
    assert Path('output').is_dir()


def test_npt_retains_original_barostat_state(run_inputs):
    theory, settings = run_inputs
    assert md.prepare_md_run(theory, settings, 'npt', statefile='state.xml') == 'state.xml'
    assert ET.parse('state.xml').getroot().find('Parameters').get('MonteCarloPressure') == '1'
    assert not Path('initial_state.xml').exists()


def test_validate_only_does_not_create_output_or_rewrite_restart(run_inputs, capsys):
    theory, settings = run_inputs
    with pytest.raises(SystemExit) as stop:
        md.prepare_md_run(theory, settings, 'nvt', statefile='state.xml',
                          argv=['--validate-only'], description='Example NVT')
    assert stop.value.code == 0
    assert 'Example NVT inputs are consistent' in capsys.readouterr().out
    assert not Path('initial_state.xml').exists()
    assert not Path('output').exists()


def test_duration_override_changes_only_requested_phase(run_inputs):
    theory, settings = run_inputs
    md.prepare_md_run(theory, settings, 'nve', statefile='state.xml',
                      argv=['--simulation-time', '0.001'])
    assert settings.nve_time == 0.001
    assert settings.nvt_time == 1.0
    assert settings.npt_time == 1.0


@pytest.mark.parametrize('duration', ['nan', '0', '0.0001', '-1'])
def test_rejects_invalid_duration_before_writing_files(run_inputs, duration):
    theory, settings = run_inputs
    with pytest.raises(SystemExit) as stop:
        md.prepare_md_run(theory, settings, 'nvt', statefile='state.xml',
                          argv=['--simulation-time', duration])
    assert stop.value.code == 2
    assert not Path('initial_state.xml').exists()
    assert not Path('output').exists()


def test_output_protection_is_applied_before_restart_conversion(run_inputs):
    theory, settings = run_inputs
    Path('output').mkdir()
    Path('output/nvt.dcd').write_bytes(b'previous results')
    with pytest.raises(SystemExit) as stop:
        md.prepare_md_run(theory, settings, 'nvt', statefile='state.xml')
    assert stop.value.code == 2
    assert Path('output/nvt.dcd').read_bytes() == b'previous results'
    assert not Path('initial_state.xml').exists()
