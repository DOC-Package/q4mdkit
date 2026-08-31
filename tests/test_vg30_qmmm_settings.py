from pathlib import Path
import subprocess
import sys

import pytest

from q4mdkit.qmmm.qmmm_config import QMMMConfig


REPO_ROOT = Path(__file__).resolve().parents[1]
VG30_DIR = REPO_ROOT / "examples" / "vg30"


def test_qmmm_uses_cationic_vg30_with_opv_b3lyp_parameters(monkeypatch):
    npt_dir = VG30_DIR / "npt"
    monkeypatch.chdir(npt_dir)

    config = QMMMConfig("qmmm_settings.yaml")

    assert config.qm_backend == "dftb"
    assert config.qm_charge == 1
    assert config.qm_mult == 1
    assert Path(config.sk_dir) == Path("/home/takahashi/qmmm/sk/opv/b3lyp")
    assert config.hybrid == {
        "method": "CAM",
        "params": {"Screening": {"method": "MatrixBased"}},
    }

    elements = ("C", "H", "N", "O", "Cl")
    expected_pairs = {f"{first}-{second}" for first in elements for second in elements}
    assert set(config.slater_koster_files) == expected_pairs
    assert all(Path(path).is_file() for path in config.slater_koster_files.values())

    assert Path(config.amber_prmtop).name == "vg30_acetonitrile.prmtop"
    assert Path(config.pdbfile).name == "npt_lastframe.pdb"
    assert Path(config.qatoms_file).name == "qmatoms"
    assert Path(config.boxfile).name == "npt_lastframe.box"


def test_npt_final_state_is_converted_to_the_qmmm_box_without_precision_loss(tmp_path):
    state_path = tmp_path / "OpenMM_MD_final_state.xml"
    state_path.write_text(
        """<?xml version='1.0' encoding='utf-8'?>
<State>
  <PeriodicBoxVectors>
    <A x="4.4255050604395985" y="0" z="0" />
    <B x="0" y="4.4255050604395985" z="0" />
    <C x="0" y="0" z="4.4255050604395985" />
  </PeriodicBoxVectors>
</State>
"""
    )
    box_path = tmp_path / "npt_lastframe.box"

    result = subprocess.run(
        [
            sys.executable,
            str(VG30_DIR / "npt-mm" / "sync_box.py"),
            str(state_path),
            str(box_path),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    vectors = [
        [float(value) for value in line.split()]
        for line in box_path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    expected = (
        [44.2550506044, 0.0, 0.0],
        [0.0, 44.2550506044, 0.0],
        [0.0, 0.0, 44.2550506044],
    )
    for vector, expected_vector in zip(vectors, expected):
        assert vector == pytest.approx(expected_vector, abs=1e-6)


def test_qmmm_npt_validates_from_repo_root_before_importing_ash():
    result = subprocess.run(
        [
            sys.executable,
            str(VG30_DIR / "npt" / "npt.py"),
            "--validate-only",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "VG30 OPV-B3LYP QM/MM inputs are consistent" in result.stdout
    assert "../npt-mm/OpenMM_MD_final_state.xml" in result.stdout


def test_mm_npt_validates_from_repo_root_before_importing_ash():
    result = subprocess.run(
        [
            sys.executable,
            str(VG30_DIR / "npt-mm" / "npt.py"),
            "--validate-only",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "VG30 MM NPT inputs are consistent" in result.stdout
