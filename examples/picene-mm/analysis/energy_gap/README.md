# Amber/OpenMM energy-gap workflow

This workflow evaluates the vertical MM energy gap

\[
\Delta E(t)=U_{\mathrm{charged}}(\mathbf R_t)-U_{\mathrm{neutral}}(\mathbf R_t)
\]

on the existing neutral trajectory. It does not minimize structures, propagate
dynamics, or create a charged-state trajectory. The charged topology must differ
from the neutral topology only in the partial charges of one selected residue.

All commands below are run from the repository root.

## Automated Amber stage

The complete cation-charge and charged-topology preparation can be inspected
without running external calculations:

```bash
bash examples/picene-mm/analysis/energy_gap/prepare_cation.sh --dry-run
```

Run it with the defaults for a +1 doublet on one-based residue 202:

```bash
bash examples/picene-mm/analysis/energy_gap/prepare_cation.sh
```

This creates `picene_cation.mol2` and `picene_cation.prmtop` in the input
directory by default. The script uses ORCA's open-shell ROHF-AM1 calculation
for the cation, converts its Mulliken charges with Antechamber, applies Amber's
`am1bcc` bond-charge corrections, and restores the original atom names. ORCA
must be available on `PATH`, or can be selected with `--orca PATH` (and
`ORCA_EXE`). Antechamber and ORCA run in a temporary directory so that
intermediate files do not pollute the repository. Existing outputs are refused
unless `--force` is supplied. Both outputs remain in the temporary directory
until ParmEd validation passes, so a failed calculation does not publish a
partial or unvalidated result. Use `--help` to see input/output and
charge-state overrides.

The script assumes conda environments named `amber` (Antechamber and
`am1bcc`) and `parmed` (the topology helper). Change those environment names in
the script if they differ locally. Path normalization uses GNU `realpath -m`,
as provided by the Linux environment used for this example.

## Run the energy gap with the wrapper

`run_energy_gap.sh` supplies the paths and timing defaults for this example and
then runs `compute_energy_gap.py` in the `myash` environment. Inspect the exact
command without starting a calculation:

```bash
bash examples/picene-mm/analysis/energy_gap/run_energy_gap.sh --dry-run
```

Run a short ten-frame check:

```bash
bash examples/picene-mm/analysis/energy_gap/run_energy_gap.sh \
    --stop 10 \
    --output /tmp/picene_energy_gap_test.csv
```

For the full saved trajectory, omit `--stop` (the default output is
`analysis/energy_gap/energy_gap.csv`). Use `--stride N` for a pilot calculation and
`--force` only when intentionally replacing an existing output. All options,
including overrides for the topology, trajectory, residue, timing, OpenMM
platform, and conda environment, are listed by `--help`.

The following sections show the equivalent individual steps.

## 1. Prepare charged-state charges

The automated script first asks Antechamber to create an ORCA input scaffold,
then `prepare_orca_input.py` changes it to a single-point ROHF-AM1 calculation
with the requested charge and multiplicity. ORCA's output is converted to
Mulliken-charge MOL2, passed through `am1bcc`, and converted back to MOL2 with
the original atom names. Amber SQM is not used here: its AM1 implementation
accepts singlets only, so the direct `antechamber -c bcc -m 2` command fails for
this radical cation.

The resulting MOL2 must have unique atom names matching the selected Amber
residue. Only its charges are transferred; atom types and other parameters are
validated independently against the authoritative topology.

## 2. Create the charged topology

The current `qmatoms` selection contains zero-based atoms 7236--7271. This is
Amber/ParmEd residue 202 in the script's one-based `--resid` convention.

```bash
conda run -n parmed python \
    examples/picene-mm/analysis/energy_gap/prepare_charged_topology.py \
    --prmtop examples/picene-mm/input/picene.prmtop \
    --charged-mol2 target_charged.mol2 \
    --resid 202 \
    --expected-delta-charge 1 \
    --charge-tolerance 0.01 \
    --output examples/picene-mm/analysis/energy_gap/picene_cation.prmtop
```

The script maps atoms by unique name, checks atom counts and atomic numbers,
prints every charge replacement, writes the new topology, reloads it, and checks
that every Amber topology FLAG except `CHARGE` is identical. It also verifies
that all non-target charges remain unchanged. `--expected-delta-charge` is
optional and is not hard-coded to +1. For this picene topology the neutral
target residue has a small Amber-rounded charge of `+0.006`, so a
`0.01` charge tolerance is used to accept the resulting `+0.994` charge
change while still rejecting a mis-assigned target.

## 3. Validate an existing topology pair

```bash
conda run -n parmed python \
    examples/picene-mm/analysis/energy_gap/validate_topologies.py \
    --neutral-prmtop examples/picene-mm/input/picene.prmtop \
    --charged-prmtop examples/picene-mm/analysis/energy_gap/picene_cation.prmtop \
    --resid 202 \
    --expected-delta-charge 1 \
    --charge-tolerance 0.01
```

Do not proceed unless the status is `PASS`.

## 4. Run the neutral-vs-neutral sanity check

The original MM calculation used PME, a 9 Å cutoff, Ewald error tolerance
`5e-4`, HBonds constraints, non-rigid water, and the CPU platform. These are the
energy script defaults. The DCD stores triclinic box vectors per frame; the
script applies the same frame box to both OpenMM contexts.

The DCD reader does not recover the physical timestamps. From `nvt.csv`, frame
0 is at 0.05 ps and the saved-frame interval is 0.05 ps, so both time options
are specified explicitly:

```bash
conda run -n myash python \
    examples/picene-mm/analysis/energy_gap/compute_energy_gap.py \
    --neutral-prmtop examples/picene-mm/input/picene.prmtop \
    --charged-prmtop examples/picene-mm/input/picene.prmtop \
    --trajectory examples/picene-mm/nvt-mm/output/nvt.dcd \
    --resid 202 \
    --start 0 \
    --stop 10 \
    --time-origin-ps 0.05 \
    --time-step-ps 0.05 \
    --platform CPU \
    --threads 1 \
    --output neutral_vs_neutral.csv
```

All gaps must be zero within a stated numerical tolerance; `1e-3 kJ/mol` is
used by the integration test for this system. CPU 1-thread evaluation is
recommended for validation. Context-dependent floating-point reductions can
give tiny residuals: observed one-frame neutral-vs-neutral results ranged from
zero to `2.44e-4 kJ/mol`.

## 5. Test the charged topology on a short interval

```bash
conda run -n myash python \
    examples/picene-mm/analysis/energy_gap/compute_energy_gap.py \
    --neutral-prmtop examples/picene-mm/input/picene.prmtop \
    --charged-prmtop examples/picene-mm/analysis/energy_gap/picene_cation.prmtop \
    --trajectory examples/picene-mm/nvt-mm/output/nvt.dcd \
    --resid 202 \
    --start 0 \
    --stop 10 \
    --stride 1 \
    --time-origin-ps 0.05 \
    --time-step-ps 0.05 \
    --platform CPU \
    --threads 1 \
    --output energy_gap_test.csv
```

Confirm that all energies and gaps are finite and that the mean of
`gap_fluctuation_kJmol` is zero within floating-point precision.

## 6. Run the selected trajectory

Remove `--stop`, or select a range and stride appropriate for the analysis:

```bash
conda run -n myash python \
    examples/picene-mm/analysis/energy_gap/compute_energy_gap.py \
    --neutral-prmtop examples/picene-mm/input/picene.prmtop \
    --charged-prmtop examples/picene-mm/analysis/energy_gap/picene_cation.prmtop \
    --trajectory examples/picene-mm/nvt-mm/output/nvt.dcd \
    --resid 202 \
    --time-origin-ps 0.05 \
    --time-step-ps 0.05 \
    --platform CPU \
    --threads 1 \
    --output energy_gap.csv
```

For a faster multi-thread or GPU production run, compare several frames against
the CPU 1-thread result first. CUDA and OpenCL precision can be selected with
`--precision single`, `mixed`, or `double`; unsupported platform properties are
not applied and the script prints the properties actually used.

The CSV contains frame/time, both total potential energies, the energy gap, its
mean-centered fluctuation, and both gap columns converted to cm\(^{-1}\).
Before constructing either OpenMM System, the energy script independently
checks all non-`CHARGE` Amber FLAGs and verifies that charge differences are
restricted to `--resid`. This prevents a same-size but reordered or otherwise
incompatible charged topology from being evaluated silently.
