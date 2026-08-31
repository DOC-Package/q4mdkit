import os
from pathlib import Path
import shutil
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
VG30_INPUT = REPO_ROOT / "examples" / "vg30" / "input-sub"


def find_amber_bin():
    configured_bin = os.environ.get("AMBER_BIN")
    configured_home = os.environ.get("AMBERHOME")
    candidates = []
    if configured_bin:
        candidates.append(Path(configured_bin))
    if configured_home:
        candidates.append(Path(configured_home) / "bin")

    for amber_bin in candidates:
        if (amber_bin / "antechamber").is_file() and (
            amber_bin / "parmchk2"
        ).is_file():
            return amber_bin

    antechamber = shutil.which("antechamber")
    parmchk2 = shutil.which("parmchk2")
    if antechamber and parmchk2:
        return Path(antechamber).parent
    return None


def test_make_gaff_generates_vg30_parameter_files(tmp_path):
    amber_bin = find_amber_bin()
    if amber_bin is None:
        pytest.skip("AmberTools is not available via AMBER_BIN, AMBERHOME, or PATH")

    shutil.copy(VG30_INPUT / "vg30-6-cl.xyz", tmp_path)
    shutil.copy(VG30_INPUT / "make_gaff.sh", tmp_path)

    env = os.environ.copy()
    env["PATH"] = f"{amber_bin}:{env['PATH']}"
    result = subprocess.run(
        ["/bin/sh", "make_gaff.sh"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    mol2_path = tmp_path / "vg30-6-cl_gaff2.mol2"
    frcmod_path = tmp_path / "vg30-6-cl_gaff2.frcmod"
    pdb_path = tmp_path / "vg30-6-cl.pdb"
    assert mol2_path.is_file() and mol2_path.stat().st_size > 0
    assert frcmod_path.is_file() and frcmod_path.stat().st_size > 0
    assert pdb_path.is_file()

    mol2_lines = mol2_path.read_text().splitlines()
    atom_start = mol2_lines.index("@<TRIPOS>ATOM") + 1
    bond_start = mol2_lines.index("@<TRIPOS>BOND")
    atom_records = mol2_lines[atom_start:bond_start]
    assert len(atom_records) == 77
    assert sum(float(line.split()[8]) for line in atom_records) == pytest.approx(
        1.0, abs=0.01
    )

    pdb_atoms = [
        line
        for line in pdb_path.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    assert len(pdb_atoms) == 77
