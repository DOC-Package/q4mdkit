#!/usr/bin/env bash
set -euo pipefail

opt_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
orca_binary=${1:-${ORCA_BINARY:-orca}}
nprocs=${ORCA_NPROCS:-8}
maxcore_mb=${ORCA_MAXCORE_MB:-3000}

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

template="$opt_dir/vg30_opt.inp.template"
input="$opt_dir/vg30_opt.inp"
output="$opt_dir/vg30_opt.out"
geometry="$opt_dir/vg30_opt.xyz"

sed \
  -e "s/@NPROCS@/$nprocs/g" \
  -e "s/@MAXCORE@/$maxcore_mb/g" \
  "$template" > "$input"

if [[ -e "$geometry" ]]; then
  backup="$geometry.previous.$(date +%Y%m%d%H%M%S).$$"
  mv -- "$geometry" "$backup"
  echo "Previous optimized geometry moved to: $backup"
fi

echo "Running ORCA geometry optimization with $nprocs processes..."
if ! (cd "$opt_dir" && "$orca_binary" "$(basename -- "$input")" > "$(basename -- "$output")" 2>&1); then
  echo "ORCA geometry optimization failed. See: $output" >&2
  exit 1
fi

if ! grep -q "ORCA TERMINATED NORMALLY" "$output"; then
  echo "ORCA did not terminate normally. See: $output" >&2
  exit 1
fi
if ! grep -q "THE OPTIMIZATION HAS CONVERGED" "$output"; then
  echo "ORCA terminated, but the geometry optimization did not converge." >&2
  echo "See: $output" >&2
  exit 1
fi
if [[ ! -s "$geometry" ]]; then
  echo "ORCA did not produce opt/vg30_opt.xyz." >&2
  exit 1
fi

echo "ORCA geometry optimization completed. Result: $geometry"
