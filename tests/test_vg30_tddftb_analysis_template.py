from pathlib import Path

import pytest

pytest.importorskip("hsd")

from q4mdkit.analysis.tddftb import prepare_tddftb_hsd


REPO_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = REPO_ROOT / "examples" / "vg30" / "analysis" / "input" / "dftb_tddft-b3lyp.hsd"


def casida_block(text):
    start = text.index("Casida {")
    return text[start:text.index("\n  }\n", start)]


def test_rendered_input_excludes_carboxyl_lone_pair_transitions(tmp_path):
    hsd_file = prepare_tddftb_hsd(TEMPLATE, tmp_path, n_excitations=3)
    casida = casida_block(hsd_file.read_text())

    assert "NrOfExcitations = 3" in casida
    assert "ExcludeOccupiedOrbitals {" in casida
    assert "Atoms = O" in casida
    assert "Population = 0.5" in casida
