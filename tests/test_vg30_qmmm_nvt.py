"""The VG30 DFTB NVT configuration opts into QM X–H constraints."""

import sys
from pathlib import Path
from types import SimpleNamespace

from q4mdkit.qmmm.qmmm_config import QMMMConfig


def test_vg30_nvt_factory_restores_qm_hbond_constraints(monkeypatch):
    nvt_dir = Path(__file__).resolve().parents[1] / "examples" / "vg30" / "nvt"
    monkeypatch.chdir(nvt_dir)
    config = QMMMConfig("qmmm_settings.yaml")

    class System:
        def __init__(self):
            self.constraints = [(0, 1, 0.11), (77, 78, 0.14)]

        def getNumConstraints(self):
            return len(self.constraints)

        def getConstraintParameters(self, index):
            return self.constraints[index]

        def addConstraint(self, *constraint):
            self.constraints.append(constraint)

    system = System()
    omm = SimpleNamespace(system=system)
    fragment = SimpleNamespace(elems=["C", "H"] + ["C"] * 77)
    result = object()

    def qmmm_theory(**kwargs):
        assert kwargs['mm_theory'] is omm
        system.constraints = [(77, 78, 0.14)]
        return result

    monkeypatch.setitem(sys.modules, "ash", SimpleNamespace(QMMMTheory=qmmm_theory))
    assert config.create_qmmm_theory(
        fragment, config.load_qmatoms(), omm=omm, qm_theory=object()
    ) is result
    assert system.constraints == [(77, 78, 0.14), (0, 1, 0.11)]
