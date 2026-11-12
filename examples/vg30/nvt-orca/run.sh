#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Resolve the active ASH Python before changing the MPI search path.
python_bin=$(command -v "${PYTHON:-python3}")
mpi_bindir=${ORCA_MPI_BINDIR:-/usr/bin}
if ! mpi_version=$("$mpi_bindir/mpirun" --version 2>&1); then
  echo "Cannot run OpenMPI from $mpi_bindir/mpirun. Set ORCA_MPI_BINDIR." >&2
  exit 1
fi
if [[ "$mpi_version" != *"Open MPI"* && "$mpi_version" != *"OpenRTE"* ]]; then
  echo "ORCA requires OpenMPI; $mpi_bindir/mpirun is not OpenMPI." >&2
  exit 1
fi
export PATH="$mpi_bindir:$PATH"
exec "$python_bin" "$script_dir/nvt.py" "$@"
