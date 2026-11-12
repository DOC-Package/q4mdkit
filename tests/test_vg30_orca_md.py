"""Exercise ORCA example wiring without invoking external QM calculations."""

import runpy
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
VG30_DIR = REPO_ROOT / "examples" / "vg30"


@pytest.fixture
def workflow(tmp_path):
    root = tmp_path / "vg30"
    for relative in (
        "input/qmatoms",
        "input/vg30_acetonitrile.prmtop",
        "nvt-orca/nvt.py", "nve-orca/nve.py",
        "nvt-orca/qmmm_settings.yaml", "nve-orca/qmmm_settings.yaml",
        "nvt-orca/md_settings.yaml", "nve-orca/md_settings.yaml",
    ):
        source = VG30_DIR / relative
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.is_file():
            shutil.copy2(source, target)

    # Synthetic restart fixtures: no generated simulation outputs are required.
    npt_output = root / "npt-mm/output"
    npt_output.mkdir(parents=True)
    pdb = ["CRYST1   44.255   44.255   44.255  90.00  90.00  90.00 P 1           1\n"]
    for index in range(5597):
        residue = "VG3" if index < 77 else "ACE"
        pdb.append(
            f"ATOM  {index + 1:5d}  C   {residue} A   1    "
            "   1.000   2.000   3.000  1.00  0.00           C\n"
        )
    (npt_output / "npt_lastframe.pdb").write_text("".join(pdb) + "END\n")
    state = ET.Element("State", type="State", version="1")
    box = ET.SubElement(state, "PeriodicBoxVectors")
    for name, axis in (("A", "x"), ("B", "y"), ("C", "z")):
        values = dict(x="0", y="0", z="0")
        values[axis] = "4.4255050604395985"
        ET.SubElement(box, name, **values)
    ET.SubElement(state, "Parameters", MonteCarloPressure="1", MonteCarloTemperature="300")
    for group, particle in (("Positions", "Position"), ("Velocities", "Velocity")):
        node = ET.SubElement(state, group)
        for _ in range(5597):
            ET.SubElement(node, particle, x="0.1", y="0.2", z="0.3")
    ET.ElementTree(state).write(root / "npt-mm/OpenMM_MD_final_state.xml")

    (root / "nvt-orca/output").mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / "npt-mm/output/npt_lastframe.pdb",
                 root / "nvt-orca/output/nvt_lastframe.pdb")
    shutil.copy2(root / "npt-mm/OpenMM_MD_final_state.xml",
                 root / "nvt-orca/OpenMM_MD_final_state.xml")
    return root


def validate(workflow, stage):
    return subprocess.run(
        [sys.executable, str(workflow / f"{stage}-orca/{stage}.py"), "--validate-only"],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )


@pytest.mark.parametrize("stage", ["nvt", "nve"])
def test_validate_from_repo_root_without_ash(workflow, stage):
    result = validate(workflow, stage)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "inputs are consistent" in result.stdout
    assert not (workflow / f"{stage}-orca/initial_state.xml").exists()


def test_nve_requires_nvt_velocities(workflow):
    state = workflow / "nvt-orca/OpenMM_MD_final_state.xml"
    tree = ET.parse(state)
    tree.getroot().remove(tree.getroot().find("Velocities"))
    tree.write(state)
    result = validate(workflow, "nve")
    assert result.returncode != 0
    assert "Velocities" in result.stderr


def test_rejects_incorrect_qm_selection(workflow):
    (workflow / "input/qmatoms").write_text("0 1 2\n")
    result = validate(workflow, "nvt")
    assert result.returncode != 0
    assert "QM region has 3 atoms; expected 77" in result.stderr


def test_nve_requires_predecessor_state(workflow):
    (workflow / "nvt-orca/OpenMM_MD_final_state.xml").unlink()
    result = validate(workflow, "nve")
    assert result.returncode != 0
    assert "nvt-orca/OpenMM_MD_final_state.xml" in result.stderr


@pytest.mark.parametrize("group,particle", [("Positions", "Position"), ("Velocities", "Velocity")])
@pytest.mark.parametrize("value", [None, "nan", "invalid"])
def test_rejects_incomplete_or_nonfinite_restart_arrays(workflow, group, particle, value):
    state = workflow / "nvt-orca/OpenMM_MD_final_state.xml"
    tree = ET.parse(state)
    node = tree.getroot().find(f"{group}/{particle}")
    if value is None:
        del node.attrib["x"]
    else:
        node.set("x", value)
    tree.write(state)
    result = validate(workflow, "nve")
    assert result.returncode != 0
    assert group in result.stderr


@pytest.mark.parametrize("existing", ["output/nve.dcd", "OpenMM_MD_final_state.xml"])
def test_rejects_existing_stage_results_before_ash(workflow, existing):
    result_file = workflow / "nve-orca" / existing
    result_file.parent.mkdir(parents=True, exist_ok=True)
    result_file.write_text("previous result")
    result = subprocess.run(
        [sys.executable, str(workflow / "nve-orca/nve.py")],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "Move existing results" in result.stderr
    assert result_file.read_text() == "previous result"
    assert not (workflow / "nve-orca/initial_state.xml").exists()


@pytest.mark.parametrize("stage,integrator", [
    ("nvt", "LangevinMiddleIntegrator"), ("nve", "VerletIntegrator"),
])
def test_orca_md_preserves_restart_box_velocities_and_constraints(
    workflow, monkeypatch, stage, integrator,
):
    calls = {}

    class System:
        constraints = [(0, 1, 0.11), (77, 78, 0.11)]

        def getNumConstraints(self):
            return len(self.constraints)

        def getConstraintParameters(self, index):
            return self.constraints[index]

        def addConstraint(self, *constraint):
            self.constraints.append(constraint)

    system = System()

    def openmm_theory(**kwargs):
        calls["openmm"] = kwargs
        return SimpleNamespace(system=system)

    def orca_theory(**kwargs):
        calls["orca"] = kwargs
        return SimpleNamespace()

    def qmmm_theory(**kwargs):
        calls["qmmm"] = kwargs
        # ASH removes QM constraints while constructing QMMMTheory.
        system.constraints = [(77, 78, 0.11)]
        return SimpleNamespace()

    ash = SimpleNamespace(
        Fragment=lambda **kwargs: SimpleNamespace(
            elems=["C", "H"] + ["C"] * 5595, numatoms=5597,
        ),
        OpenMMTheory=openmm_theory, ORCATheory=orca_theory,
        QMMMTheory=qmmm_theory,
        OpenMM_MD=lambda **kwargs: calls.update(md=kwargs),
    )
    monkeypatch.setitem(sys.modules, "ash", ash)
    monkeypatch.syspath_prepend(str(workflow))
    monkeypatch.chdir(REPO_ROOT)
    if stage == "nvt":
        # NVT execution starts fresh; the fixture's NVT result is for NVE tests.
        (workflow / "nvt-orca/OpenMM_MD_final_state.xml").unlink()
    script = workflow / f"{stage}-orca/{stage}.py"
    monkeypatch.setattr(sys, "argv", [str(script), "--simulation-time", "0.001"])
    # The stage intentionally changes its working directory, as ASH uses it.
    runpy.run_path(str(script), run_name="__main__")

    assert "B3LYP" in calls["orca"]["orcasimpleinput"]
    assert calls["orca"]["numcores"] == 8
    assert (workflow / f"{stage}-orca/output/orca_energy.dat").is_file()
    assert calls["qmmm"]["embedding"] == "elstat"
    assert calls["qmmm"]["qm_charge"] == 1
    assert calls["qmmm"]["qm_mult"] == 1
    assert calls["qmmm"]["qmatoms"] == list(range(77))
    assert calls["openmm"]["PBCvectors"][0][0] == pytest.approx(
        44.255050604395985, abs=1e-12,
    )
    assert system.constraints == [(77, 78, 0.11), (0, 1, 0.11)]
    assert calls["md"]["integrator"] == integrator
    assert calls["md"]["timestep"] == 0.0005
    assert calls["md"]["simulation_time"] == 0.001
    assert "barostat" not in calls["md"]
    if stage == "nve":
        assert "temperature" not in calls["md"]
        source = workflow / "nvt-orca/OpenMM_MD_final_state.xml"
    else:
        assert calls["md"]["temperature"] == 300
        source = workflow / "npt-mm/OpenMM_MD_final_state.xml"
    original = ET.parse(source).getroot()
    converted = ET.parse(calls["md"]["statefile"]).getroot()
    for node in ("Positions", "Velocities", "PeriodicBoxVectors"):
        assert ET.tostring(converted.find(node)) == ET.tostring(original.find(node))
    assert converted.find("Parameters") is None
