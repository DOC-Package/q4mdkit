import csv
import importlib.util
import math
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
ENERGY_GAP_DIR = (
    ROOT / "examples" / "pentacene775a-mm" / "analysis" / "energy_gap"
)
PENTACENE_PRMTOP = (
    ROOT / "examples" / "pentacene775a-mm" / "input" / "pentacene.prmtop"
)
PENTACENE_DCD = (
    ROOT / "examples" / "pentacene775a-mm" / "nve-mm" / "output" / "nve.dcd"
)
PREPARE_CATION_SH = ENERGY_GAP_DIR / "prepare_cation.sh"
PREPARE_ORCA_INPUT = ENERGY_GAP_DIR / "prepare_orca_input.py"
RUN_ENERGY_GAP_SH = ENERGY_GAP_DIR / "run_energy_gap.sh"


def load_script(name):
    path = ENERGY_GAP_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"pentacene_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atom(name, charge, atomic_number, atom_type="ca", mass=12.01):
    return SimpleNamespace(
        name=name,
        charge=charge,
        atomic_number=atomic_number,
        type=atom_type,
        mass=mass,
        idx=-1,
    )


def structure(residues, parm_data=None):
    atoms = []
    for residue in residues:
        for item in residue.atoms:
            item.idx = len(atoms)
            atoms.append(item)
    return SimpleNamespace(
        residues=residues,
        atoms=atoms,
        parm_data=parm_data or {},
    )


def test_transfer_charges_maps_unique_atom_names_without_copying_types():
    topology = load_script("prepare_charged_topology")
    target = SimpleNamespace(
        name="PEN",
        atoms=[atom("C1", -0.2, 6, "ca"), atom("H2", 0.2, 1, "ha", 1.008)],
    )
    neutral = structure([target])
    charged_mol2 = structure(
        [
            SimpleNamespace(
                name="PEN",
                atoms=[atom("H2", 0.6, 1, "hx", 1.008), atom("C1", 0.4, 6, "cx")],
            )
        ]
    )

    changes = topology.transfer_charges_by_name(neutral, charged_mol2, resid_1based=1)

    assert [item.name for item in changes] == ["C1", "H2"]
    assert [item.neutral_charge for item in changes] == pytest.approx([-0.2, 0.2])
    assert [item.charged_charge for item in changes] == pytest.approx([0.4, 0.6])
    assert [item.type for item in neutral.atoms] == ["ca", "ha"]
    assert sum(item.charge for item in neutral.atoms) == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("mol2_atoms", "message"),
    [
        ([atom("C1", 0.5, 6), atom("C1", 0.5, 6)], "unique"),
        ([atom("C1", 0.5, 6), atom("X2", 0.5, 1)], "atom names"),
        ([atom("C1", 0.5, 7), atom("H2", 0.5, 1)], "atomic number"),
    ],
)
def test_transfer_charges_rejects_ambiguous_or_incompatible_mapping(mol2_atoms, message):
    topology = load_script("prepare_charged_topology")
    neutral = structure(
        [SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6), atom("H2", 0.2, 1)])]
    )
    charged_mol2 = structure([SimpleNamespace(name="PEN", atoms=mol2_atoms)])

    with pytest.raises(ValueError, match=message):
        topology.transfer_charges_by_name(neutral, charged_mol2, resid_1based=1)


def test_topology_validation_accepts_only_target_charge_changes():
    topology = load_script("prepare_charged_topology")
    neutral = structure(
        [
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)]),
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)]),
        ],
        parm_data={"POINTERS": [2, 1], "ATOM_NAME": ["C1", "C1"], "CHARGE": [-3.64446, -3.64446]},
    )
    charged = structure(
        [
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)]),
            SimpleNamespace(name="PEN", atoms=[atom("C1", 0.8, 6)]),
        ],
        parm_data={"POINTERS": [2, 1], "ATOM_NAME": ["C1", "C1"], "CHARGE": [-3.64446, 14.57784]},
    )

    report = topology.validate_topologies(neutral, charged, resid_1based=2)

    assert report.status == "PASS"
    assert report.neutral_target_charge == pytest.approx(-0.2)
    assert report.charged_target_charge == pytest.approx(0.8)
    assert report.delta_target_charge == pytest.approx(1.0)


def test_topology_validation_rejects_nontarget_charge_change():
    topology = load_script("prepare_charged_topology")
    neutral = structure(
        [
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)]),
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)]),
        ],
        parm_data={"ATOM_NAME": ["C1", "C1"], "CHARGE": [-3.64446, -3.64446]},
    )
    charged = structure(
        [
            SimpleNamespace(name="PEN", atoms=[atom("C1", -0.1, 6)]),
            SimpleNamespace(name="PEN", atoms=[atom("C1", 0.8, 6)]),
        ],
        parm_data={"ATOM_NAME": ["C1", "C1"], "CHARGE": [-1.82223, 14.57784]},
    )

    with pytest.raises(ValueError, match="Non-target charge"):
        topology.validate_topologies(neutral, charged, resid_1based=2)


def test_topology_validation_rejects_noncharge_parameter_change():
    topology = load_script("prepare_charged_topology")
    neutral = structure(
        [SimpleNamespace(name="PEN", atoms=[atom("C1", -0.2, 6)])],
        parm_data={"ATOM_NAME": ["C1"], "MASS": [12.01], "CHARGE": [-3.64446]},
    )
    charged = structure(
        [SimpleNamespace(name="PEN", atoms=[atom("C1", 0.8, 6)])],
        parm_data={"ATOM_NAME": ["C1"], "MASS": [13.01], "CHARGE": [14.57784]},
    )

    with pytest.raises(ValueError, match="MASS"):
        topology.validate_topologies(neutral, charged, resid_1based=1)


def test_prepare_rejects_overwriting_authoritative_neutral_topology(tmp_path):
    topology = load_script("prepare_charged_topology")
    neutral = tmp_path / "system.prmtop"

    with pytest.raises(ValueError, match="overwrite the neutral topology"):
        topology.prepare_charged_topology(
            neutral,
            tmp_path / "charged.mol2",
            resid_1based=1,
            output=neutral,
        )


def test_center_energy_gaps_has_zero_mean_and_hand_checked_values():
    energy_gap = load_script("compute_energy_gap")

    mean, centered = energy_gap.center_energy_gaps([10.0, 12.0, 14.0])

    assert mean == pytest.approx(12.0)
    assert centered == pytest.approx([-2.0, 0.0, 2.0])
    assert np.mean(centered) == pytest.approx(0.0, abs=1.0e-14)


def test_frame_time_override_uses_original_frame_index():
    energy_gap = load_script("compute_energy_gap")

    assert energy_gap.frame_time_ps(7, reader_time_ps=99.0, origin_ps=600.004, step_ps=0.004) == pytest.approx(
        600.032
    )
    assert energy_gap.frame_time_ps(7, reader_time_ps=99.0, origin_ps=None, step_ps=None) == 99.0


def test_write_results_csv_includes_centered_gap_and_wavenumbers(tmp_path):
    energy_gap = load_script("compute_energy_gap")
    output = tmp_path / "energy_gap.csv"
    evaluations = [
        energy_gap.EnergyEvaluation(0, 0.0, -100.0, -90.0),
        energy_gap.EnergyEvaluation(1, 0.004, -101.0, -89.0),
    ]

    energy_gap.write_results_csv(output, evaluations)

    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [float(row["gap_kJmol"]) for row in rows] == pytest.approx([10.0, 12.0])
    assert [float(row["gap_fluctuation_kJmol"]) for row in rows] == pytest.approx([-1.0, 1.0])
    assert [float(row["gap_cm-1"]) for row in rows] == pytest.approx(
        [835.934, 1003.121], abs=0.001
    )


def test_configure_periodic_amber_topology_supports_prmtop_without_ifbox():
    energy_gap = load_script("compute_energy_gap")

    class FakeOpenMMTopology:
        def __init__(self):
            self.box_vectors = None

        def setPeriodicBoxVectors(self, vectors):
            self.box_vectors = vectors

    amber = SimpleNamespace(
        topology=FakeOpenMMTopology(),
        _prmtop=SimpleNamespace(_raw_data={"POINTERS": [0] * 32}),
    )
    box_vectors_nm = np.asarray(
        [[4.397, 0.0, 0.0], [0.516, 5.380, 0.0], [0.251, 1.640, 7.035]]
    )

    energy_gap.configure_periodic_amber_topology(
        amber, box_vectors_nm, nanometer_unit=1.0
    )

    assert amber._prmtop._raw_data["POINTERS"][27] == 1
    assert amber.topology.box_vectors == pytest.approx(box_vectors_nm)


def amber_loader(charges, atom_names=("C1", "H2", "C1", "H2")):
    return SimpleNamespace(
        _raw_data={
            "POINTERS": [4, 2],
            "RESIDUE_LABEL": ["PEN", "PEN"],
            "RESIDUE_POINTER": ["1", "3"],
            "ATOM_NAME": list(atom_names),
            "MASS": ["12.01", "1.008", "12.01", "1.008"],
            "CHARGE": [str(value) for value in charges],
        }
    )


def test_energy_evaluation_topology_guard_accepts_only_target_charge_changes():
    energy_gap = load_script("compute_energy_gap")
    neutral = amber_loader([-1.0, 1.0, -1.0, 1.0])
    charged = amber_loader([-1.0, 1.0, 8.0, 10.2223])

    report = energy_gap.validate_amber_topology_data(
        neutral, charged, resid_1based=2
    )

    assert report.changed_target_atoms == 2
    assert report.delta_target_charge == pytest.approx(1.0)


def test_energy_evaluation_topology_guard_rejects_reordering():
    energy_gap = load_script("compute_energy_gap")
    neutral = amber_loader([-1.0, 1.0, -1.0, 1.0])
    charged = amber_loader(
        [-1.0, 1.0, 8.0, 10.2223], atom_names=("C1", "H2", "H2", "C1")
    )

    with pytest.raises(ValueError, match="ATOM_NAME"):
        energy_gap.validate_amber_topology_data(neutral, charged, resid_1based=2)


def test_energy_evaluation_topology_guard_rejects_nontarget_charge_change():
    energy_gap = load_script("compute_energy_gap")
    neutral = amber_loader([-1.0, 1.0, -1.0, 1.0])
    charged = amber_loader([-0.5, 1.0, 8.0, 10.2223])

    with pytest.raises(ValueError, match="Non-target charge"):
        energy_gap.validate_amber_topology_data(neutral, charged, resid_1based=2)


@pytest.mark.skipif(
    importlib.util.find_spec("openmm") is None
    or not PENTACENE_PRMTOP.exists()
    or not PENTACENE_DCD.exists(),
    reason="requires myash OpenMM and the pentacene NVE trajectory",
)
def test_openmm_neutral_vs_neutral_is_zero_with_frame_box():
    energy_gap = load_script("compute_energy_gap")

    evaluations, _ = energy_gap.evaluate_energy_gap(
        PENTACENE_PRMTOP,
        PENTACENE_PRMTOP,
        PENTACENE_DCD,
        resid_1based=156,
        start=0,
        stop=1,
        time_origin_ps=600.004,
        time_step_ps=0.004,
        platform_name="CPU",
        threads=1,
    )

    assert len(evaluations) == 1
    assert evaluations[0].time_ps == pytest.approx(600.004)
    assert evaluations[0].gap_kjmol == pytest.approx(0.0, abs=1.0e-3)
    assert math.isfinite(evaluations[0].neutral_kjmol)


def test_prepare_cation_dry_run_builds_amber_and_parmed_commands(tmp_path):
    result = subprocess.run(
        ["bash", str(PREPARE_CATION_SH), "--dry-run"],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "ORCA AM1 input preparation:" in result.stdout
    assert "prepare_orca_input.py" in result.stdout
    assert "ORCA open-shell AM1 calculation:" in result.stdout
    orca_path = shutil.which("orca")
    assert orca_path is not None
    assert orca_path in result.stdout
    assert "Antechamber AM1-BCC conversion:" in result.stdout
    assert "conda run -n amber antechamber" in result.stdout
    assert str(ROOT / "examples/pentacene775a-mm/input/pentacene_gaff2.mol2") in result.stdout
    assert "-nc 1" in result.stdout
    assert "-m 2" in result.stdout
    assert "-at gaff2" in result.stdout
    assert "-fo orcinp" in result.stdout
    assert "-fi orcout" in result.stdout
    assert "-c mul" in result.stdout
    assert "conda run -n amber am1bcc" in result.stdout
    assert "-ao name" in result.stdout
    assert "-c bcc" not in result.stdout
    assert "-an n -du n -seq n" in result.stdout
    assert "conda run -n parmed python" in result.stdout
    assert "--resid 156" in result.stdout
    assert "pentacene_cation.mol2" in result.stdout
    assert "pentacene_cation.prmtop" in result.stdout


def test_prepare_orca_input_replaces_antechamber_method_and_spin(tmp_path):
    source = tmp_path / "generated.inp"
    output = tmp_path / "am1.inp"
    source.write_text(
        """%MaxCore 4000
! HF 6-31G(d) TightSCF TightOpt KeepDens
%id \"molecule\"

* xyz 0 1
C 0.0 0.0 0.0
H 0.0 0.0 1.0
*
"""
    )

    result = subprocess.run(
        [
            "python",
            str(PREPARE_ORCA_INPUT),
            "--input",
            str(source),
            "--output",
            str(output),
            "--charge",
            "1",
            "--multiplicity",
            "2",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    converted = output.read_text()
    assert "! AM1 ROHF SlowConv TightSCF Mulliken" in converted
    assert "TightOpt" not in converted
    assert "%scf\n  MaxIter 500\nend" in converted
    assert "* xyz 1 2" in converted
    assert "6-31G" not in converted


def test_run_energy_gap_dry_run_builds_myash_command(tmp_path):
    output = tmp_path / "energy_gap.csv"
    result = subprocess.run(
        [
            "bash",
            str(RUN_ENERGY_GAP_SH),
            "--dry-run",
            "--neutral-prmtop",
            str(PENTACENE_PRMTOP),
            "--charged-prmtop",
            str(ROOT / "examples/pentacene775a-mm/input/pentacene_cation.prmtop"),
            "--trajectory",
            str(PENTACENE_DCD),
            "--output",
            str(output),
            "--resid",
            "156",
            "--start",
            "0",
            "--stop",
            "10",
            "--stride",
            "2",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "conda run -n myash python" in result.stdout
    assert "compute_energy_gap.py" in result.stdout
    assert "--resid 156" in result.stdout
    assert "--start 0" in result.stdout
    assert "--stop 10" in result.stdout
    assert "--stride 2" in result.stdout
    assert str(output) in result.stdout


def test_prepare_cation_refuses_to_overwrite_input_even_with_force():
    input_mol2 = ROOT / "examples/pentacene775a-mm/input/pentacene_gaff2.mol2"
    result = subprocess.run(
        [
            "bash",
            str(PREPARE_CATION_SH),
            "--dry-run",
            "--force",
            "--cation-mol2",
            str(input_mol2),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "Refusing to overwrite the input MOL2" in result.stderr


@pytest.mark.parametrize(
    ("output_option", "protected_input"),
    [
        (
            "--cation-mol2",
            ROOT / "examples/pentacene775a-mm/input/pentacene.prmtop",
        ),
        (
            "--cation-prmtop",
            ROOT / "examples/pentacene775a-mm/input/pentacene_gaff2.mol2",
        ),
    ],
)
def test_prepare_cation_refuses_cross_type_input_overwrites(
    output_option, protected_input
):
    result = subprocess.run(
        [
            "bash",
            str(PREPARE_CATION_SH),
            "--dry-run",
            "--force",
            output_option,
            str(protected_input),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "Refusing to overwrite protected input" in result.stderr


def test_prepare_cation_stages_both_outputs_until_validation_succeeds(tmp_path):
    cation_mol2 = tmp_path / "result.mol2"
    cation_prmtop = tmp_path / "result.prmtop"
    result = subprocess.run(
        [
            "bash",
            str(PREPARE_CATION_SH),
            "--dry-run",
            "--cation-mol2",
            str(cation_mol2),
            "--cation-prmtop",
            str(cation_prmtop),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--output /tmp/q4mdkit-cation.XXXXXX/pentacene_cation.prmtop" in result.stdout
    assert f"mv -f /tmp/q4mdkit-cation.XXXXXX/pentacene_cation.mol2 {cation_mol2}" in result.stdout
    assert f"mv -f /tmp/q4mdkit-cation.XXXXXX/pentacene_cation.prmtop {cation_prmtop}" in result.stdout


def test_prepare_cation_does_not_publish_outputs_when_validation_fails(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_conda = fake_bin / "conda"
    fake_conda.write_text(
        """#!/usr/bin/env bash
set -eu
if [[ " $* " == *" antechamber "* ]]; then
    while [[ $# -gt 0 ]]; do
        if [[ $1 == -o ]]; then
            printf 'test mol2' > "$2"
            exit 0
        fi
        shift
    done
fi
exit 23
"""
    )
    fake_conda.chmod(0o755)
    cation_mol2 = tmp_path / "result.mol2"
    cation_prmtop = tmp_path / "result.prmtop"
    environment = dict(os.environ)
    environment["PATH"] = f"{fake_bin}:{environment['PATH']}"

    result = subprocess.run(
        [
            "bash",
            str(PREPARE_CATION_SH),
            "--cation-mol2",
            str(cation_mol2),
            "--cation-prmtop",
            str(cation_prmtop),
        ],
        text=True,
        capture_output=True,
        check=False,
        env=environment,
    )

    assert result.returncode == 23
    assert not cation_mol2.exists()
    assert not cation_prmtop.exists()
