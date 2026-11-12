#!/usr/bin/env bash
set -euo pipefail

excited_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
orca_binary=${1:-${ORCA_BINARY:-orca}}
nprocs=${ORCA_NPROCS:-8}
maxcore_mb=${ORCA_MAXCORE_MB:-3000}
geometry="$excited_dir/../opt/vg30_opt.xyz"

if [[ ! -s "$geometry" ]]; then
  echo "Converged optimized geometry not found: $geometry" >&2
  echo "Run opt/run.sh first." >&2
  exit 1
fi
if [[ ! "$nprocs" =~ ^[1-9][0-9]*$ ]]; then
  echo "ORCA_NPROCS must be a positive integer: $nprocs" >&2
  exit 2
fi
if [[ ! "$maxcore_mb" =~ ^[1-9][0-9]*$ ]]; then
  echo "ORCA_MAXCORE_MB must be a positive integer: $maxcore_mb" >&2
  exit 2
fi

if [[ "$orca_binary" == */* ]]; then
  [[ -x "$orca_binary" ]] || {
    echo "ORCA executable is not executable: $orca_binary" >&2
    exit 2
  }
  orca_binary=$(cd -- "$(dirname -- "$orca_binary")" && pwd -P)/$(basename -- "$orca_binary")
else
  orca_binary=$(command -v "$orca_binary") || {
    echo "ORCA executable not found: $orca_binary" >&2
    echo "Pass its path as the first argument or set ORCA_BINARY." >&2
    exit 2
  }
fi

template="$excited_dir/vg30_tddft.inp.template"
input="$excited_dir/vg30_tddft.inp"
output="$excited_dir/vg30_tddft.out"

sed \
  -e "s/@NPROCS@/$nprocs/g" \
  -e "s/@MAXCORE@/$maxcore_mb/g" \
  "$template" > "$input"

echo "Running ORCA full TDDFT for ten singlet states with $nprocs processes..."
if ! (cd "$excited_dir" && "$orca_binary" "$(basename -- "$input")" > "$(basename -- "$output")" 2>&1); then
  echo "ORCA TDDFT calculation failed. See: $output" >&2
  exit 1
fi

if ! grep -q "ORCA TERMINATED NORMALLY" "$output"; then
  echo "ORCA did not terminate normally. See: $output" >&2
  exit 1
fi

echo "ORCA TDDFT calculation completed. Results: $output"
