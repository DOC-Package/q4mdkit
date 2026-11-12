from pathlib import Path
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_DIR = REPO_ROOT / "examples" / "vg30" / "tddftb"


def test_vacuum_inputs_exclude_mm_point_charges():
    for relative_path in ("opt/dftb_in.hsd", "excited/dftb_in.hsd"):
        text = (WORKFLOW_DIR / relative_path).read_text()

        assert "ElectricField" not in text
        assert "PointCharges" not in text
        assert "PCcharges.dat" not in text


def test_optimization_input_describes_cationic_opv_b3lyp_geometry_optimization():
    text = (WORKFLOW_DIR / "opt" / "dftb_in.hsd").read_text()

    assert '<<< "../input/vg30.xyz"' in text
    assert "Charge = 1" in text
    assert "Scc = Yes" in text
    assert "Hybrid = CAM" in text
    assert "Driver = GeometryOptimization" in text
    assert 'OutputPrefix = "geo_end"' in text
    assert "/home/takahashi/qmmm/sk/opv/b3lyp/C-Cl.skf" in text


def test_tddftb_input_uses_optimized_geometry_and_ten_singlet_states():
    text = (WORKFLOW_DIR / "excited" / "dftb_in.hsd").read_text()

    assert '<<< "../opt/geo_end.xyz"' in text
    assert "Charge = 1" in text
    assert "Hybrid = None" in text
    assert "/home/takahashi/qmmm/sk/3ob/C-Cl.skf" in text
    assert "/home/takahashi/qmmm/sk/opv/b3lyp/" not in text
    assert "Casida" in text
    assert "NrOfExcitations = 10" in text
    assert "Symmetry = Singlet" in text
    assert "Diagonaliser = Stratmann" in text


def test_starting_geometry_contains_the_77_vg30_atoms():
    lines = (WORKFLOW_DIR / "input" / "vg30.xyz").read_text().splitlines()

    assert int(lines[0]) == 77
    atom_lines = [line for line in lines[2:] if line.strip()]
    assert len(atom_lines) == 77
    assert {line.split()[0] for line in atom_lines} == {"C", "H", "N", "O", "Cl"}
    assert all(len(line.split()) == 4 for line in atom_lines)


def make_fake_dftb(tmp_path):
    fake_dftb = tmp_path / "fake-dftb"
    fake_dftb.write_text(
        """#!/bin/sh
set -eu
printf '%s\\n' "$(basename "$PWD")" >> "$DFTB_CALL_LOG"
case "$(basename "$PWD")" in
  opt)
    cp ../input/vg30.xyz geo_end.xyz
    ;;
  excited)
    test -f ../opt/geo_end.xyz
    printf 'fake excitation output\\n' > EXC.DAT
    ;;
  *)
    exit 90
    ;;
esac
"""
    )
    fake_dftb.chmod(0o755)
    return fake_dftb


def test_optimization_runner_only_runs_optimization(tmp_path):
    workflow = tmp_path / "tddftb"
    shutil.copytree(WORKFLOW_DIR, workflow)
    call_log = tmp_path / "calls.log"
    fake_dftb = make_fake_dftb(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "opt" / "run.sh"), str(fake_dftb)],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "DFTB_CALL_LOG": str(call_log)},
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert call_log.read_text().splitlines() == ["opt"]
    assert (workflow / "opt" / "geo_end.xyz").is_file()
    assert "Geometry optimization completed" in result.stdout


def test_excited_runner_only_runs_excited_state_calculation(tmp_path):
    workflow = tmp_path / "tddftb"
    shutil.copytree(WORKFLOW_DIR, workflow)
    shutil.copyfile(workflow / "input" / "vg30.xyz", workflow / "opt" / "geo_end.xyz")
    call_log = tmp_path / "calls.log"
    fake_dftb = make_fake_dftb(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "excited" / "run.sh"), str(fake_dftb)],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "DFTB_CALL_LOG": str(call_log)},
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert call_log.read_text().splitlines() == ["excited"]
    assert (workflow / "excited" / "EXC.DAT").is_file()
    assert "TD-DFTB calculation completed" in result.stdout


def test_excited_runner_requires_optimized_geometry(tmp_path):
    workflow = tmp_path / "tddftb"
    shutil.copytree(WORKFLOW_DIR, workflow)
    (workflow / "opt" / "geo_end.xyz").unlink(missing_ok=True)
    fake_dftb = make_fake_dftb(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "excited" / "run.sh"), str(fake_dftb)],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "DFTB_CALL_LOG": str(tmp_path / "calls.log")},
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Run opt/run.sh first" in result.stderr
