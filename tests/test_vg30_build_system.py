import os
from pathlib import Path
import shutil
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def find_amber_path():
    current_path = os.environ.get("PATH", "")
    for configured in (
        os.environ.get("AMBER_BIN"),
        str(Path(os.environ["AMBERHOME"]) / "bin")
        if os.environ.get("AMBERHOME")
        else None,
    ):
        if configured:
            amber_bin = Path(configured)
            search_path = f"{amber_bin}:{current_path}"
            if all(shutil.which(tool, path=search_path) for tool in ("packmol", "tleap")):
                return search_path

    if all(shutil.which(tool) for tool in ("packmol", "tleap")):
        return current_path
    return None


def read_prmtop_charge(prmtop_path):
    lines = prmtop_path.read_text().splitlines()
    charge_flag = next(
        index for index, line in enumerate(lines) if line.strip() == "%FLAG CHARGE"
    )
    scaled_charges = []
    for line in lines[charge_flag + 2 :]:
        if line.lstrip().startswith("%FLAG"):
            break
        scaled_charges.extend(float(value) for value in line.split())
    return sum(scaled_charges) / 18.2223


def test_builds_safely_padded_vg30_in_acetonitrile(tmp_path):
    amber_path = find_amber_path()
    if amber_path is None:
        pytest.skip("Packmol and tleap are not available")

    examples_dir = tmp_path / "examples"
    vg30_dir = examples_dir / "vg30"
    input_dir = vg30_dir / "input"
    input_dir.mkdir(parents=True)
    for name in ("build_system.sh", "make_system.in", "pack_system.inp", "make_qmatoms.py"):
        shutil.copy(REPO_ROOT / "examples" / "vg30" / "input" / name, input_dir)
    shutil.copytree(
        REPO_ROOT / "examples" / "vg30" / "input-sub", vg30_dir / "input-sub"
    )
    shutil.copytree(
        REPO_ROOT / "examples" / "solvent" / "acetonitrile",
        examples_dir / "solvent" / "acetonitrile",
    )
    prep_dir = tmp_path / "q4mdkit" / "prep"
    prep_dir.mkdir(parents=True)
    shutil.copy(REPO_ROOT / "q4mdkit" / "prep" / "pdb2box.py", prep_dir)

    env = os.environ.copy()
    env["PATH"] = amber_path
    env["PYTHONPATH"] = f"{REPO_ROOT}:{env.get('PYTHONPATH', '')}"
    result = subprocess.run(
        ["/bin/bash", "build_system.sh"],
        cwd=vg30_dir / "input",
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    output_names = (
        "vg30_acetonitrile.pdb",
        "vg30_acetonitrile.prmtop",
        "vg30_acetonitrile.inpcrd",
        "vg30_acetonitrile_leap.pdb",
        "vg30_acetonitrile.box",
        "qmatoms",
    )
    for name in output_names:
        output = vg30_dir / "input" / name
        assert output.is_file() and output.stat().st_size > 0

    prmtop_path = vg30_dir / "input" / "vg30_acetonitrile.prmtop"
    assert read_prmtop_charge(prmtop_path) == pytest.approx(1.001, abs=1e-6)

    qmatoms = [int(index) for index in (vg30_dir / "input" / "qmatoms").read_text().split()]
    assert qmatoms == list(range(77))

    leap_pdb = (vg30_dir / "input" / "vg30_acetonitrile_leap.pdb").read_text()
    atom_lines = [
        line for line in leap_pdb.splitlines() if line.startswith(("ATOM  ", "HETATM"))
    ]
    assert len(atom_lines) == 77 + 920 * 6
    assert {line[17:20].strip() for line in atom_lines[:77]} == {"VG3"}
    assert {line[17:20].strip() for line in atom_lines[77:]} == {"ACE"}

    box_lines = [
        line
        for line in (vg30_dir / "input" / "vg30_acetonitrile.box")
        .read_text()
        .splitlines()
        if not line.startswith("#")
    ]
    box_vectors = [[float(value) for value in line.split()] for line in box_lines]
    expected_box_vectors = (
        [43.2, 0.0, 0.0],
        [0.0, 43.2, 0.0],
        [0.0, 0.0, 43.2],
    )
    assert len(box_vectors) == len(expected_box_vectors)
    for vector, expected in zip(box_vectors, expected_box_vectors):
        assert vector == pytest.approx(expected, abs=1e-6)
    box_lengths = [box_vectors[index][index] for index in range(3)]

    vg30_coordinates = [
        tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        for line in atom_lines[:77]
    ]
    boundary_clearances = [
        min(
            min(coordinate[axis] for coordinate in vg30_coordinates),
            box_lengths[axis]
            - max(coordinate[axis] for coordinate in vg30_coordinates),
        )
        for axis in range(3)
    ]
    assert min(boundary_clearances) >= 9.9

    centroid = [
        sum(coordinate[axis] for coordinate in vg30_coordinates)
        / len(vg30_coordinates)
        for axis in range(3)
    ]
    assert centroid == pytest.approx([21.6, 21.6, 21.6], abs=1e-3)


def test_make_qmatoms_fails_without_selected_residue(tmp_path):
    pdb_path = tmp_path / "system.pdb"
    pdb_path.write_text(
        "HETATM    1  C1  ACE A   1       0.000   0.000   0.000  1.00  0.00          C\n"
    )
    output_path = tmp_path / "qmatoms"
    output_path.write_text("0 1 2")

    result = subprocess.run(
        [
            "python3",
            str(REPO_ROOT / "examples" / "vg30" / "input" / "make_qmatoms.py"),
            str(pdb_path),
            "--resname",
            "VG3",
            "-o",
            str(output_path),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert not output_path.exists()
