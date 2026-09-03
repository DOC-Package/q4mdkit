#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
example_dir="$(cd -- "$script_dir/../.." && pwd -P)"

neutral_prmtop="$example_dir/input/pentacene.prmtop"
charged_prmtop="$example_dir/input/pentacene_cation.prmtop"
trajectory="$example_dir/nve-mm/output/nve.dcd"
output="$example_dir/nve-mm/output/energy_gap.csv"
resid=156
start=0
stop=""
stride=1
time_origin_ps=600.004
time_step_ps=0.004
platform=CPU
precision=default
threads=1
nonbonded_method=PME
cutoff_angstrom=9.0
ewald_error_tolerance=5.0e-4
chunk_size=100
conda_env=myash
dry_run=false
force=false

usage() {
    cat <<'EOF'
Usage: run_energy_gap.sh [options]

Run compute_energy_gap.py for the pentacene Amber/OpenMM trajectory.

Options:
  --neutral-prmtop PATH       Neutral Amber topology
  --charged-prmtop PATH       Cation Amber topology
  --trajectory PATH           OpenMM trajectory (DCD)
  --output PATH               Output CSV
  --resid N                   1-based target residue (default: 156)
  --start N                   First frame (default: 0)
  --stop N                    Exclusive final frame (default: all frames)
  --stride N                  Frame stride (default: 1)
  --time-origin-ps VALUE      Time of frame 0 (default: 600.004)
  --time-step-ps VALUE        Saved-frame interval (default: 0.004)
  --platform NAME             OpenMM platform (default: CPU)
  --precision NAME             OpenMM precision (default: default)
  --threads N                 CPU threads (default: 1)
  --nonbonded-method NAME     PME, Ewald, CutoffPeriodic, or NoCutoff
  --cutoff-angstrom VALUE     Periodic cutoff (default: 9.0)
  --ewald-error-tolerance V   Ewald tolerance (default: 5.0e-4)
  --chunk-size N              Trajectory chunk size (default: 100)
  --conda-env NAME            Conda environment (default: myash)
  --force                     Overwrite an existing output CSV
  --dry-run                   Print the command without running it
  -h, --help                  Show this help
EOF
}

require_value() {
    if [[ $# -lt 2 || -z ${2:-} ]]; then
        echo "Missing value for $1" >&2
        usage >&2
        exit 2
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --neutral-prmtop) require_value "$@"; neutral_prmtop=$2; shift 2 ;;
        --charged-prmtop) require_value "$@"; charged_prmtop=$2; shift 2 ;;
        --trajectory) require_value "$@"; trajectory=$2; shift 2 ;;
        --output) require_value "$@"; output=$2; shift 2 ;;
        --resid) require_value "$@"; resid=$2; shift 2 ;;
        --start) require_value "$@"; start=$2; shift 2 ;;
        --stop) require_value "$@"; stop=$2; shift 2 ;;
        --stride) require_value "$@"; stride=$2; shift 2 ;;
        --time-origin-ps) require_value "$@"; time_origin_ps=$2; shift 2 ;;
        --time-step-ps) require_value "$@"; time_step_ps=$2; shift 2 ;;
        --platform) require_value "$@"; platform=$2; shift 2 ;;
        --precision) require_value "$@"; precision=$2; shift 2 ;;
        --threads) require_value "$@"; threads=$2; shift 2 ;;
        --nonbonded-method) require_value "$@"; nonbonded_method=$2; shift 2 ;;
        --cutoff-angstrom) require_value "$@"; cutoff_angstrom=$2; shift 2 ;;
        --ewald-error-tolerance) require_value "$@"; ewald_error_tolerance=$2; shift 2 ;;
        --chunk-size) require_value "$@"; chunk_size=$2; shift 2 ;;
        --conda-env) require_value "$@"; conda_env=$2; shift 2 ;;
        --force) force=true; shift ;;
        --dry-run) dry_run=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

is_nonnegative_integer() { [[ $1 =~ ^[0-9]+$ ]]; }
if ! is_nonnegative_integer "$resid" || (( resid < 1 )); then
    echo "--resid must be a positive integer: $resid" >&2; exit 2
fi
if ! is_nonnegative_integer "$start"; then
    echo "--start must be a non-negative integer: $start" >&2; exit 2
fi
if [[ -n $stop ]]; then
    if ! is_nonnegative_integer "$stop" || (( stop <= start )); then
        echo "--stop must be an integer greater than --start: $stop" >&2; exit 2
    fi
fi
for numeric_option in stride threads chunk_size; do
    numeric_value=${!numeric_option}
    if ! is_nonnegative_integer "$numeric_value" || (( numeric_value < 1 )); then
        echo "--${numeric_option//_/-} must be a positive integer: $numeric_value" >&2
        exit 2
    fi
done

neutral_prmtop="$(realpath -m -- "$neutral_prmtop")"
charged_prmtop="$(realpath -m -- "$charged_prmtop")"
trajectory="$(realpath -m -- "$trajectory")"
output="$(realpath -m -- "$output")"
for input_path in "$neutral_prmtop" "$charged_prmtop" "$trajectory"; do
    if [[ ! -f $input_path ]]; then
        echo "Input file not found: $input_path" >&2
        exit 1
    fi
done
if [[ $output == "$neutral_prmtop" || $output == "$charged_prmtop" || $output == "$trajectory" ]]; then
    echo "Refusing to overwrite an input file: $output" >&2
    exit 2
fi
if [[ -e $output && $force != true ]]; then
    echo "Output already exists (use --force to overwrite): $output" >&2
    exit 1
fi

compute_command=(
    conda run -n "$conda_env" python "$script_dir/compute_energy_gap.py"
    --neutral-prmtop "$neutral_prmtop"
    --charged-prmtop "$charged_prmtop"
    --trajectory "$trajectory"
    --output "$output"
    --resid "$resid"
    --start "$start"
    --stride "$stride"
    --time-origin-ps "$time_origin_ps"
    --time-step-ps "$time_step_ps"
    --platform "$platform"
    --precision "$precision"
    --threads "$threads"
    --nonbonded-method "$nonbonded_method"
    --cutoff-angstrom "$cutoff_angstrom"
    --ewald-error-tolerance "$ewald_error_tolerance"
    --chunk-size "$chunk_size"
)
if [[ -n $stop ]]; then
    compute_command+=(--stop "$stop")
fi

echo "Energy-gap command:"
printf '%q ' "${compute_command[@]}"
printf '\n'

if [[ $dry_run == true ]]; then
    exit 0
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "conda was not found in PATH; use --dry-run or activate Conda first" >&2
    exit 1
fi
mkdir -p -- "$(dirname -- "$output")"
"${compute_command[@]}"
