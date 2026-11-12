"""Check the launcher environment without starting ASH or an MD run."""

import os
from pathlib import Path
import subprocess

import pytest


VG30_DIR = Path(__file__).resolve().parents[1] / "examples" / "vg30"


def executable(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\nset -eu\n" + body)
    path.chmod(0o755)
    return path


@pytest.mark.parametrize("launcher,script", [
    ("nvt-orca/run.sh", "nvt-orca/nvt.py"),
    ("nve-orca/run.sh", "nve-orca/nve.py"),
    ("analysis/run_tdorca.sh", "analysis/run_tdorca.py"),
])
def test_launcher_selects_openmpi_and_preserves_python(tmp_path, launcher, script):
    selected_python = executable(
        tmp_path / "ash/bin/python3",
        'printf "python=%s\\n" "$0"\n'
        'printf "mpi=%s\\n" "$(command -v mpirun)"\n'
        'printf "arg=%s\\n" "$@"\n',
    )
    mpi = executable(tmp_path / "openmpi/bin/mpirun", 'echo "mpirun (Open MPI) 4.1.6"\n')
    executable(mpi.parent / "python3", "exit 91\n")
    intel = executable(tmp_path / "intel/bin/mpirun", 'echo "Intel MPI"\n')
    env = dict(os.environ)
    env.pop("PYTHON", None)
    env.update(PATH=f"{selected_python.parent}:{intel.parent}:/usr/bin:/bin",
               ORCA_MPI_BINDIR=str(mpi.parent))

    result = subprocess.run(
        ["/bin/bash", str(VG30_DIR / launcher),
         "--simulation-time", "0.001"],
        env=env, cwd=tmp_path, capture_output=True, text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert f"python={selected_python}" in result.stdout
    assert f"mpi={mpi}" in result.stdout
    assert f"arg={VG30_DIR / script}" in result.stdout
    assert "arg=--simulation-time\narg=0.001\n" in result.stdout


@pytest.mark.parametrize("launcher", [
    "nvt-orca/run.sh", "nve-orca/run.sh", "analysis/run_tdorca.sh",
])
def test_launcher_rejects_intel_mpi_override(tmp_path, launcher):
    mpi = executable(tmp_path / "intel/bin/mpirun", 'echo "Intel MPI"\n')
    python = executable(tmp_path / "ash/bin/python3", 'echo "MD started"\n')
    env = dict(os.environ, PYTHON=str(python), ORCA_MPI_BINDIR=str(mpi.parent))

    result = subprocess.run(
        ["/bin/bash", str(VG30_DIR / launcher)],
        env=env, cwd=tmp_path, capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "OpenMPI" in result.stderr
    assert "MD started" not in result.stdout
