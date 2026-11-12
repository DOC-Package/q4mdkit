from pathlib import Path
import os
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_DIR = REPO_ROOT / "examples" / "vg30" / "orca"


def make_fake_orca(tmp_path):
    fake_orca = tmp_path / "fake-orca"
    fake_orca.write_text(
        r"""#!/usr/bin/env bash
set -eu
input=$1
stage=$(basename "$PWD")
printf '%s\n' "$stage" >> "$ORCA_CALL_LOG"
grep -q 'B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF' "$input"
grep -q '^  nprocs 6$' "$input"
grep -q '^%maxcore 2500$' "$input"
if grep -Eiq 'pointcharges|pcfile|qmmm' "$input"; then
  exit 91
fi
case "$stage" in
  opt)
    grep -q ' Opt' "$input"
    grep -q '\* xyzfile 1 1 ../input/vg30.xyz' "$input"
    cp ../input/vg30.xyz vg30_opt.xyz
    printf '%s\n' '***        THE OPTIMIZATION HAS CONVERGED     ***'
    printf '%s\n' '****ORCA TERMINATED NORMALLY****'
    ;;
  excited)
    grep -q '^  nroots 10$' "$input"
    grep -q '^  tda false$' "$input"
    grep -q '^  triplets false$' "$input"
    grep -q '\* xyzfile 1 1 ../opt/vg30_opt.xyz' "$input"
    test -s ../opt/vg30_opt.xyz
    printf '%s\n' 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS'
    printf '%s\n' '****ORCA TERMINATED NORMALLY****'
    ;;
  *)
    exit 90
    ;;
esac
"""
    )
    fake_orca.chmod(0o755)
    return fake_orca


def runner_environment(tmp_path):
    environment = dict(os.environ)
    environment.update(
        {
            "ORCA_CALL_LOG": str(tmp_path / "calls.log"),
            "ORCA_NPROCS": "6",
            "ORCA_MAXCORE_MB": "2500",
        }
    )
    return environment


def test_starting_geometry_is_coordinate_only_vg30_cation_structure():
    lines = (WORKFLOW_DIR / "input" / "vg30.xyz").read_text().splitlines()

    assert int(lines[0]) == 77
    atoms = [line.split() for line in lines[2:] if line.strip()]
    assert len(atoms) == 77
    assert all(len(atom) == 4 for atom in atoms)
    assert atoms[0] == ["C", "8.558826", "0.855104", "-0.183722"]


def test_optimization_runner_builds_orca_input_and_optimized_geometry(tmp_path):
    workflow = tmp_path / "orca"
    shutil.copytree(WORKFLOW_DIR, workflow)
    fake_orca = make_fake_orca(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "opt" / "run.sh"), str(fake_orca)],
        cwd=tmp_path,
        env=runner_environment(tmp_path),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "calls.log").read_text().splitlines() == ["opt"]
    assert (workflow / "opt" / "vg30_opt.xyz").is_file()
    assert "ORCA geometry optimization completed" in result.stdout


def test_excited_runner_builds_full_tddft_input_for_ten_singlets(tmp_path):
    workflow = tmp_path / "orca"
    shutil.copytree(WORKFLOW_DIR, workflow)
    shutil.copyfile(workflow / "input" / "vg30.xyz", workflow / "opt" / "vg30_opt.xyz")
    fake_orca = make_fake_orca(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "excited" / "run.sh"), str(fake_orca)],
        cwd=tmp_path,
        env=runner_environment(tmp_path),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "calls.log").read_text().splitlines() == ["excited"]
    assert "ORCA TDDFT calculation completed" in result.stdout
    output = (workflow / "excited" / "vg30_tddft.out").read_text()
    assert "ABSORPTION SPECTRUM" in output


def test_excited_runner_requires_converged_optimized_geometry(tmp_path):
    workflow = tmp_path / "orca"
    shutil.copytree(WORKFLOW_DIR, workflow)
    fake_orca = make_fake_orca(tmp_path)

    result = subprocess.run(
        ["bash", str(workflow / "excited" / "run.sh"), str(fake_orca)],
        cwd=tmp_path,
        env=runner_environment(tmp_path),
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Run opt/run.sh first" in result.stderr
    assert not (tmp_path / "calls.log").exists()
