import runpy
import sys
from pathlib import Path
from types import SimpleNamespace


def test_qmmm_nvt_converts_state_and_restores_qm_hbond_constraints(monkeypatch):
    nvt_dir = Path(__file__).resolve().parents[1] / "examples" / "vg30" / "nvt"

    calls = []

    class FakeSystem:
        def __init__(self):
            self.constraints = [(0, 1, 0.11), (77, 78, 0.14)]

        def getNumConstraints(self):
            return len(self.constraints)

        def getConstraintParameters(self, index):
            return self.constraints[index]

        def addConstraint(self, atom1, atom2, distance):
            self.constraints.append((atom1, atom2, distance))

    system = FakeSystem()
    omm = SimpleNamespace(system=system)

    class FakeQMMMConfig:
        def load_qmatoms(self):
            return list(range(77))

        def create_openmm_theory(self, **kwargs):
            return omm

        def create_qm_theory(self):
            return object()

        def create_qmmm_theory(self, fragment, qmatoms, **kwargs):
            system.constraints.clear()
            return object()

    class FakeMDConfig:
        def run_nvt(self, fragment, theory, **kwargs):
            calls.append(("run_nvt", kwargs["statefile"]))

        def save_final_structure(self, **kwargs):
            calls.append(("save", kwargs["prefix"]))

    class FakeFragment:
        elems = ["C", "H"] + ["C"] * 77

        def __init__(self, pdbfile):
            self.pdbfile = pdbfile

    def remove_montecarlo_xml(source, destination):
        calls.append(("remove_montecarlo_xml", source, destination))

    modules = {
        "ash": SimpleNamespace(Fragment=FakeFragment),
        "q4mdkit.qmmm.qmmm_config": SimpleNamespace(
            get_config=lambda path: FakeQMMMConfig()
        ),
        "q4mdkit.qmmm.md_config": SimpleNamespace(
            get_md_config=lambda path: FakeMDConfig()
        ),
        "q4mdkit.qmmm.xml": SimpleNamespace(
            remove_montecarlo_xml=remove_montecarlo_xml
        ),
    }
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)

    monkeypatch.chdir(nvt_dir)
    monkeypatch.syspath_prepend(str(nvt_dir))
    monkeypatch.setattr(sys, "argv", [str(nvt_dir / "nvt.py")])

    runpy.run_path(str(nvt_dir / "nvt.py"), run_name="__main__")

    assert calls == [
        (
            "remove_montecarlo_xml",
            "../npt-mm/OpenMM_MD_final_state.xml",
            "nve_initial_state.xml",
        ),
        ("run_nvt", "nve_initial_state.xml"),
        ("save", "nvt"),
    ]
    assert system.constraints == [(0, 1, 0.11)]
