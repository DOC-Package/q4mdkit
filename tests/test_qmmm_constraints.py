"""Constraint preservation through the real QM/MM configuration factory."""

import sys
from types import SimpleNamespace

import pytest
import yaml

from q4mdkit.qmmm.qmmm_config import QMMMConfig


@pytest.mark.parametrize("preserve", [True, False, None])
def test_qmmm_factory_preserves_qm_hydrogen_constraints_only_when_enabled(tmp_path, monkeypatch, preserve):
    settings = {"qm_backend": "orca", "openmm": {}}
    if preserve is not None:
        settings['openmm']['preserve_qm_hbond_constraints'] = preserve
    path = tmp_path / "qmmm.yaml"
    path.write_text(yaml.safe_dump(settings))
    config = QMMMConfig(path)

    class System:
        constraints = [(0, 1, 0.11), (0, 2, 0.14), (3, 4, 0.10)]

        def getNumConstraints(self):
            return len(self.constraints)

        def getConstraintParameters(self, index):
            return self.constraints[index]

        def addConstraint(self, *values):
            self.constraints.append(values)

    system = System()
    omm = SimpleNamespace(system=system)
    fragment = SimpleNamespace(elems=['C', 'H', 'C', 'O', 'H'])
    result = object()

    def qmmm_theory(**kwargs):
        assert kwargs['mm_theory'] is omm
        system.constraints = [(3, 4, 0.10)]  # ASH removes all constraints on QM atoms.
        return result

    monkeypatch.setitem(sys.modules, 'ash', SimpleNamespace(QMMMTheory=qmmm_theory))
    assert config.create_qmmm_theory(fragment, [0, 1, 2], omm=omm, qm_theory=object()) is result
    expected = [(3, 4, 0.10), (0, 1, 0.11)] if preserve else [(3, 4, 0.10)]
    assert system.constraints == expected
