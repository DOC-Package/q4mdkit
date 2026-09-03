#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
example_dir=$(cd -- "$script_dir/../.." && pwd -P)

input_mol2="$example_dir/input/picene_gaff2.mol2"
neutral_prmtop="$example_dir/input/picene.prmtop"
cation_mol2="$example_dir/input/picene_cation.mol2"
cation_prmtop="$example_dir/input/picene_cation.prmtop"
resid=202
net_charge=1
multiplicity=2
orca_executable=${ORCA_EXE:-orca}
dry_run=false
force=false

usage() {
    cat <<'EOF'
Usage: prepare_cation.sh [options]

Generate radical-cation AM1-BCC charges with ORCA and AmberTools, then create
and validate an Amber topology whose selected residue differs only in charge.

Options:
  --input-mol2 PATH       Single-molecule input MOL2
  --neutral-prmtop PATH   Authoritative neutral Amber topology
  --cation-mol2 PATH      Output cation MOL2
  --cation-prmtop PATH    Output cation Amber topology
  --resid N               Target residue number, 1-based (default: 202)
  --net-charge N          Cation net charge (default: 1)
  --multiplicity N        Spin multiplicity 2S+1 (default: 2)
  --orca PATH             ORCA executable (default: ORCA_EXE or orca)
  --force                 Allow existing cation outputs to be replaced
  --dry-run               Print commands without running them
  -h, --help              Show this help
EOF
}

require_value() {
    if [[ $# -lt 2 ]]; then
        echo "Missing value for $1" >&2
        exit 2
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-mol2)
            require_value "$@"
            input_mol2=$2
            shift 2
            ;;
        --neutral-prmtop)
            require_value "$@"
            neutral_prmtop=$2
            shift 2
            ;;
        --cation-mol2)
            require_value "$@"
            cation_mol2=$2
            shift 2
            ;;
        --cation-prmtop)
            require_value "$@"
            cation_prmtop=$2
            shift 2
            ;;
        --resid)
            require_value "$@"
            resid=$2
            shift 2
            ;;
        --net-charge)
            require_value "$@"
            net_charge=$2
            shift 2
            ;;
        --multiplicity)
            require_value "$@"
            multiplicity=$2
            shift 2
            ;;
        --orca)
            require_value "$@"
            orca_executable=$2
            shift 2
            ;;
        --force)
            force=true
            shift
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

for value_name in resid net_charge multiplicity; do
    value=${!value_name}
    if [[ ! $value =~ ^-?[0-9]+$ ]]; then
        echo "$value_name must be an integer: $value" >&2
        exit 2
    fi
done
if (( resid < 1 )); then
    echo "resid must be at least 1" >&2
    exit 2
fi
if (( multiplicity < 1 )); then
    echo "multiplicity must be at least 1" >&2
    exit 2
fi

input_mol2=$(realpath -m -- "$input_mol2")
neutral_prmtop=$(realpath -m -- "$neutral_prmtop")
cation_mol2=$(realpath -m -- "$cation_mol2")
cation_prmtop=$(realpath -m -- "$cation_prmtop")

if [[ ! -f $input_mol2 ]]; then
    echo "Input MOL2 not found: $input_mol2" >&2
    exit 1
fi
if [[ ! -f $neutral_prmtop ]]; then
    echo "Neutral topology not found: $neutral_prmtop" >&2
    exit 1
fi
if [[ $cation_mol2 == "$input_mol2" ]]; then
    echo "Refusing to overwrite the input MOL2: $input_mol2" >&2
    exit 1
fi
if [[ $cation_prmtop == "$neutral_prmtop" ]]; then
    echo "Refusing to overwrite the neutral topology: $neutral_prmtop" >&2
    exit 1
fi
if [[ $cation_mol2 == "$neutral_prmtop" || $cation_prmtop == "$input_mol2" ]]; then
    echo "Refusing to overwrite protected input with a different output type" >&2
    exit 1
fi
if [[ $cation_mol2 == "$cation_prmtop" ]]; then
    echo "Cation MOL2 and prmtop outputs must be different paths" >&2
    exit 1
fi
if [[ $force == false && $dry_run == false ]]; then
    for output in "$cation_mol2" "$cation_prmtop"; do
        if [[ -e $output ]]; then
            echo "Output already exists: $output (use --force to replace it)" >&2
            exit 1
        fi
    done
fi

print_command() {
    printf '  '
    printf '%q ' "$@"
    printf '\n'
}

if resolved_orca=$(command -v -- "$orca_executable" 2>/dev/null); then
    orca_executable=$resolved_orca
elif [[ $dry_run == false ]]; then
    echo "ORCA executable was not found: $orca_executable" >&2
    exit 1
fi

if [[ $dry_run == true ]]; then
    work_dir="${TMPDIR:-/tmp}/q4mdkit-cation.XXXXXX"
else
    if ! command -v conda >/dev/null 2>&1; then
        echo "conda was not found in PATH" >&2
        exit 1
    fi
    mkdir -p -- "$(dirname -- "$cation_mol2")" "$(dirname -- "$cation_prmtop")"
    work_dir=$(mktemp -d "${TMPDIR:-/tmp}/q4mdkit-cation.XXXXXX")
    cleanup() {
        rm -rf -- "$work_dir"
    }
    trap cleanup EXIT
fi

temporary_mol2="$work_dir/picene_cation.mol2"
temporary_prmtop="$work_dir/picene_cation.prmtop"
orca_scaffold="$work_dir/picene_orca_scaffold.inp"
orca_input="$work_dir/picene_cation.inp"
orca_output="$work_dir/picene_cation.out"
am1_mol2="$work_dir/picene_am1.mol2"
am1_ac="$work_dir/picene_am1.ac"
bcc_ac="$work_dir/picene_bcc.ac"
antechamber_scaffold_command=(
    conda run -n amber antechamber
    -i "$input_mol2" -fi mol2
    -o "$orca_scaffold" -fo orcinp
    -at gaff2
    -nc "$net_charge" -m "$multiplicity" -rn PIC
    -an n -du n -seq n
)
prepare_orca_command=(
    conda run -n parmed python "$script_dir/prepare_orca_input.py"
    --input "$orca_scaffold" --output "$orca_input"
    --charge "$net_charge" --multiplicity "$multiplicity"
)
orca_command=("$orca_executable" "$orca_input")
antechamber_mulliken_command=(
    conda run -n amber antechamber
    -i "$orca_output" -fi orcout
    -o "$am1_mol2" -fo mol2
    -at gaff2 -c mul
    -nc "$net_charge" -m "$multiplicity" -rn PIC
    -an n -du n -seq n
)
antechamber_ac_command=(
    conda run -n amber antechamber
    -i "$am1_mol2" -fi mol2
    -o "$am1_ac" -fo ac
    -at bcc
    -nc "$net_charge" -m "$multiplicity" -rn PIC
    -an n -du n -seq n
)
am1bcc_command=(
    conda run -n amber am1bcc
    -i "$am1_ac" -o "$bcc_ac" -f ac -j 4
)
antechamber_output_command=(
    conda run -n amber antechamber
    -i "$bcc_ac" -fi ac
    -o "$temporary_mol2" -fo mol2
    -at gaff2
    -nc "$net_charge" -m "$multiplicity" -rn PIC
    -an n -du n -seq n
    -a "$input_mol2" -fa mol2 -ao name
)
parmed_command=(
    conda run -n parmed python "$script_dir/prepare_charged_topology.py"
    --prmtop "$neutral_prmtop"
    --charged-mol2 "$temporary_mol2"
    --resid "$resid"
    --expected-delta-charge "$net_charge"
    --charge-tolerance 0.01
    --output "$temporary_prmtop"
)
publish_mol2_command=(mv -f "$temporary_mol2" "$cation_mol2")
publish_prmtop_command=(mv -f "$temporary_prmtop" "$cation_prmtop")

echo "Antechamber ORCA input scaffold:"
print_command "${antechamber_scaffold_command[@]}"
echo "ORCA AM1 input preparation:"
print_command "${prepare_orca_command[@]}"
echo "ORCA open-shell AM1 calculation:"
print_command "${orca_command[@]}"
echo "Antechamber Mulliken-charge conversion:"
print_command "${antechamber_mulliken_command[@]}"
echo "Antechamber BCC input conversion:"
print_command "${antechamber_ac_command[@]}"
echo "Amber BCC correction:"
print_command "${am1bcc_command[@]}"
echo "Antechamber AM1-BCC conversion:"
print_command "${antechamber_output_command[@]}"
echo "Charged topology preparation and validation:"
print_command "${parmed_command[@]}"
echo "Publish validated outputs:"
print_command "${publish_mol2_command[@]}"
print_command "${publish_prmtop_command[@]}"

if [[ $dry_run == true ]]; then
    exit 0
fi

(
    cd -- "$work_dir"
    "${antechamber_scaffold_command[@]}"
)
"${prepare_orca_command[@]}"
(cd -- "$work_dir" && "${orca_command[@]}" > "$orca_output")
(cd -- "$work_dir" && "${antechamber_mulliken_command[@]}")
(cd -- "$work_dir" && "${antechamber_ac_command[@]}")
(cd -- "$work_dir" && "${am1bcc_command[@]}")
(cd -- "$work_dir" && "${antechamber_output_command[@]}")
"${parmed_command[@]}"
"${publish_mol2_command[@]}"
"${publish_prmtop_command[@]}"

echo
echo "Created cation MOL2:   $cation_mol2"
echo "Created cation prmtop: $cation_prmtop"
