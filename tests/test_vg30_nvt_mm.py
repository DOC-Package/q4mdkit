import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

from q4mdkit.mm.mm_config import MMConfig


REPO_ROOT = Path(__file__).resolve().parents[1]
VG30_DIR = REPO_ROOT / "examples" / "vg30"


def test_nvt_mm_loads_vg30_topology_and_safe_box(monkeypatch):
    nvt_dir = VG30_DIR / "nvt-mm"
    monkeypatch.chdir(nvt_dir)

    config = MMConfig("mm_settings.yaml")

    assert Path(config.amber_prmtop).name == "vg30_acetonitrile.prmtop"
    assert Path(config.pdbfile).name == "frag-minimized.pdb"
    assert Path(config.boxfile).name == "vg30_acetonitrile.box"
    expected_vectors = (
        [43.2, 0.0, 0.0],
        [0.0, 43.2, 0.0],
        [0.0, 0.0, 43.2],
    )
    for vector, expected in zip(config.pbc_vectors, expected_vectors):
        assert vector == pytest.approx(expected, abs=1e-6)


def run_preflight(pdb_path):
    return subprocess.run(
        [
            sys.executable,
            str(VG30_DIR / "nvt-mm" / "validate_inputs.py"),
            "--pdb",
            str(pdb_path),
            "--prmtop",
            str(VG30_DIR / "input" / "vg30_acetonitrile.prmtop"),
            "--box",
            str(VG30_DIR / "input" / "vg30_acetonitrile.box"),
        ],
        capture_output=True,
        text=True,
    )


def test_nvt_preflight_accepts_remimized_43_2_angstrom_system(tmp_path):
    pdb_path = tmp_path / "frag-minimized.pdb"
    pdb_path.write_text((VG30_DIR / "opt-mm" / "frag-minimized.pdb").read_text())

    result = run_preflight(pdb_path)

    assert result.returncode == 0, result.stdout + result.stderr


def test_nvt_preflight_rejects_stale_34_8_angstrom_minimization(tmp_path):
    current_pdb = (VG30_DIR / "opt-mm" / "frag-minimized.pdb").read_text()
    stale_pdb = current_pdb.replace(
        "CRYST1   43.200   43.200   43.200",
        "CRYST1   34.800   34.800   34.800",
        1,
    )
    pdb_path = tmp_path / "frag-minimized.pdb"
    pdb_path.write_text(stale_pdb)

    result = run_preflight(pdb_path)

    assert result.returncode != 0
    assert "does not match box file" in result.stderr


def test_minimization_output_is_the_configured_nvt_input(monkeypatch):
    paths_file = VG30_DIR / "workflow_paths.py"
    spec = importlib.util.spec_from_file_location("vg30_workflow_paths", paths_file)
    workflow_paths = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(workflow_paths)

    nvt_dir = VG30_DIR / "nvt-mm"
    monkeypatch.chdir(nvt_dir)
    config = MMConfig("mm_settings.yaml")

    configured_pdb = (nvt_dir / config.pdbfile).resolve()
    assert configured_pdb == workflow_paths.MINIMIZED_PDB


def test_nvt_direct_launch_from_repo_root_reaches_box_preflight():
    result = subprocess.run(
        [
            sys.executable,
            str(VG30_DIR / "nvt-mm" / "nvt.py"),
            "--validate-only",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "VG30 NVT inputs are consistent" in result.stdout
