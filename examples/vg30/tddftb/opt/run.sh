#!/usr/bin/env bash
set -euo pipefail

opt_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
dftb_binary=${1:-${DFTB_BINARY:-dftb+}}

if [[ "$dftb_binary" == */* ]]; then
  [[ -x "$dftb_binary" ]] || {
    echo "DFTB+ executable is not executable: $dftb_binary" >&2
    exit 2
  }
else
  command -v "$dftb_binary" >/dev/null 2>&1 || {
    echo "DFTB+ executable not found: $dftb_binary" >&2
    echo "Pass its path as the first argument or set DFTB_BINARY." >&2
    exit 2
  }
fi

geometry="$opt_dir/geo_end.xyz"
if [[ -e "$geometry" ]]; then
  backup="$geometry.previous.$(date +%Y%m%d%H%M%S).$$"
  mv -- "$geometry" "$backup"
  echo "Previous optimized geometry moved to: $backup"
fi

echo "Running VG30+ geometry optimization..."
(cd "$opt_dir" && "$dftb_binary" > dftb+.out 2>&1)

if [[ ! -s "$geometry" ]]; then
  echo "DFTB+ did not produce opt/geo_end.xyz." >&2
  exit 1
fi

echo "Geometry optimization completed. Result: $geometry"
