#!/usr/bin/env bash
set -euo pipefail

excited_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
dftb_binary=${1:-${DFTB_BINARY:-dftb+}}
geometry="$excited_dir/../opt/geo_end.xyz"

if [[ ! -s "$geometry" ]]; then
  echo "Optimized geometry not found: $geometry" >&2
  echo "Run opt/run.sh first." >&2
  exit 1
fi

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

echo "Running TD-DFTB/Casida on the optimized geometry..."
(cd "$excited_dir" && "$dftb_binary" > dftb+.out 2>&1)

echo "TD-DFTB calculation completed. Results are in: $excited_dir"
