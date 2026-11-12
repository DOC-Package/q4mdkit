"""Generic AMBER validation for systems other than VG30."""

from pathlib import Path

import pytest
import yaml

from q4mdkit import validation
from q4mdkit.mm.mm_config import MMConfig
from q4mdkit.qmmm.qmmm_config import QMMMConfig


@pytest.fixture
def system_files(tmp_path):
    pdb = tmp_path / "system.pdb"
    lines = ["CRYST1   30.000   31.000   32.000  90.00  90.00  90.00 P 1           1\n"]
    for index, name in enumerate(["DYE", "DYE", "SOL", "SOL"]):
        lines.append(f"ATOM  {index+1:5d}  C   {name} A   1    "
                     "   1.000   2.000   3.000  1.00  0.00           C\n")
    pdb.write_text(''.join(lines) + 'END\n')
    top = tmp_path / "system.prmtop"
    top.write_text("%VERSION VERSION_STAMP = V0001.000\n%FLAG POINTERS\n%FORMAT(10I8)\n       4\n")
    box = tmp_path / "system.box"
    box.write_text("30 0 0\n0 31 0\n0 0 32\n")
    return pdb, top, box


def test_validates_an_arbitrary_system_with_noncontiguous_qm_atoms(system_files):
    pdb, top, box = system_files
    assert validation.validate_amber_system(pdb, top, box, qm_atoms=[0, 2]) == 4
    assert validation.validate_amber_system(
        pdb, top, box, expected_atoms=4,
        residue_blocks=[{"name": "DYE", "count": 2}, {"name": "SOL", "count": 2}],
        qm_atoms=[0, 1], expected_qm_atoms=2, qm_resname="DYE",
    ) == 4


def test_rejects_topology_particle_count_mismatch(system_files):
    pdb, top, box = system_files
    top.write_text(top.read_text().replace("       4", "       5"))
    with pytest.raises(ValueError, match="topology has 5 atoms but PDB has 4"):
        validation.validate_amber_system(pdb, top, box)


@pytest.mark.parametrize("selection", [[0, 0], [-1, 1], [0, 4]])
def test_rejects_invalid_qm_indices(system_files, selection):
    pdb, top, box = system_files
    with pytest.raises(ValueError, match="QM"):
        validation.validate_amber_system(pdb, top, box, qm_atoms=selection)


def test_rejects_wrong_qm_residue_and_wrong_residue_order(system_files):
    pdb, top, box = system_files
    with pytest.raises(ValueError, match="QM"):
        validation.validate_amber_system(pdb, top, box, qm_atoms=[0, 2], qm_resname="DYE")
    with pytest.raises(ValueError, match="residue"):
        validation.validate_amber_system(pdb, top, box,
                                        residue_blocks=[{"name": "SOL", "count": 4}])


def test_checks_box_and_allows_nonperiodic_inputs(system_files):
    pdb, top, box = system_files
    box.write_text("40 0 0\n0 40 0\n0 0 40\n")
    with pytest.raises(ValueError, match="does not match box file"):
        validation.validate_amber_system(pdb, top, box)
    pdb.write_text('\n'.join(pdb.read_text().splitlines()[1:]) + '\n')
    assert validation.validate_amber_system(pdb, top, periodic=False) == 4
    with pytest.raises(ValueError, match="CRYST1"):
        validation.validate_amber_system(pdb, top)


@pytest.mark.parametrize("kind", ["mm", "qmmm"])
@pytest.mark.parametrize("box_format", ["vectors", "cell"])
def test_configuration_validation_uses_yaml_metadata(system_files, tmp_path, kind, box_format):
    pdb, top, box = system_files
    if box_format == "cell":
        box.write_text("30 31 32 90 90 90\n")
    indices = tmp_path / "qmatoms"
    indices.write_text("0 1\n")
    config_path = tmp_path / "settings.yaml"
    config_path.write_text(yaml.safe_dump({
        "paths": {"pdbfile": str(pdb), "amber_prmtop": str(top),
                  "boxfile": str(box), "qatoms": str(indices)},
        "qm_backend": "orca",
        "openmm": {"periodic": True, "use_boxfile": True},
        "validation": {"expected_atoms": 4, "expected_qm_atoms": 2, "qm_resname": "DYE"},
    }))
    cls = MMConfig if kind == "mm" else QMMMConfig
    config = cls(config_path)
    assert config.validate_inputs() == 4
    config.validation['expected_atoms'] = 5
    with pytest.raises(ValueError, match="expected 5"):
        config.validate_inputs()


def test_reads_legacy_triclinic_cell_parameters(tmp_path):
    box = tmp_path / "cell.box"
    box.write_text("30 31 32 90 90 60\n")
    a, b, c = validation.load_box_vectors(box)
    assert a == [30, 0, 0]
    assert b == pytest.approx([15.5, 26.8467875173, 0])
    assert c == pytest.approx([0, 0, 32], abs=1e-12)


@pytest.mark.parametrize("cell", ["30 31 32 90 90 0", "30 31 32 10 10 170"])
def test_rejects_degenerate_legacy_cells(tmp_path, cell):
    box = tmp_path / "cell.box"
    box.write_text(cell + "\n")
    with pytest.raises(ValueError, match="cell"):
        validation.load_box_vectors(box)
