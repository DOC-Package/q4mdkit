# Implementation Request: Energy-Gap Fluctuation from an Existing Amber/OpenMM Trajectory

## Objective

Implement a small, reliable workflow to calculate the **vertical energy-gap fluctuation** between a neutral state and a charged state for one selected molecule in an already existing MD trajectory.

The existing MD trajectory was generated for the **neutral/ground-state system**. Do **not** run new MD simulations, geometry optimizations, or minimizations.

For every trajectory snapshot \(\mathbf R_t\), calculate

\[
\Delta E(t)
=
U_{\mathrm{charged}}(\mathbf R_t)
-
U_{\mathrm{neutral}}(\mathbf R_t),
\]

and then calculate the mean-centered fluctuation

\[
\delta\Delta E(t)
=
\Delta E(t)-\langle\Delta E\rangle.
\]

The charged Hamiltonian should differ from the neutral Hamiltonian **only in the partial atomic charges of one selected molecule**.

---

# Physical model

Suppose molecule \(M\) is the molecule to be charged.

The neutral system has charges

\[
q_i^{(0)}.
\]

For atoms belonging to \(M\), replace them with charged-state partial charges

\[
q_i^{(0)} \rightarrow q_i^{(+)}.
\]

For all atoms outside \(M\),

\[
q_j^{(+)}=q_j^{(0)}.
\]

Do not change any other force-field parameters.

In particular, retain the original:

- GAFF2 atom types
- bond parameters
- angle parameters
- dihedral/torsion parameters
- Lennard-Jones parameters
- exclusions
- 1-4 scaling parameters
- solvent parameters
- ion parameters
- atom ordering
- residue ordering

Thus the neutral and charged topologies should differ only in the partial charges of the selected molecule.

Because the molecular geometry fluctuates along the trajectory, this definition includes both

\[
\Delta E_{\mathrm{intra,elec}}(t)
\]

and

\[
\Delta E_{\mathrm{molecule-env,elec}}(t).
\]

This is intentional. I want the intramolecular MM electrostatic contribution to the energy-gap fluctuation to be retained.

Do not simplify the calculation to only

\[
\sum_i \Delta q_i \phi_i^{\mathrm{env}}.
\]

Use the full MM potential-energy difference.

---

# Software environments

There are two existing Conda environments.

## Amber / AmberTools environment

```bash
conda activate amber
```

AmberTools, Antechamber, and ParmEd should be accessed from this environment.

Prefer reproducible commands such as

```bash
conda run -n amber ...
```

when appropriate.

## OpenMM environment

```bash
conda activate myash
```

OpenMM is installed in this environment.

Use this environment for all OpenMM energy evaluations.

Prefer

```bash
conda run -n myash ...
```

when appropriate.

Before implementing, inspect the actual installed versions:

```bash
conda run -n amber python -c "import parmed; print(parmed.__version__)"
conda run -n myash python -c "import openmm; print(openmm.__version__)"
```

Also check whether MDTraj or MDAnalysis is already installed in `myash`.

Do not modify the Conda environments automatically unless absolutely necessary. If a required package is missing, report it first.

---

# Force field

The target molecule uses **GAFF2**.

The existing neutral Amber topology should be treated as authoritative for all GAFF2 bonded and nonbonded parameters.

The charged topology should therefore NOT be rebuilt from scratch with `tleap`.

Instead:

1. start from the existing neutral `prmtop`;
2. identify the selected target residue/molecule;
3. replace only its atomic partial charges;
4. write a new `charged.prmtop`.

This avoids accidentally changing atom types, bonded parameters, exclusions, or atom ordering.

---

# Charged-state partial charges

Charged-state partial charges may be supplied as a MOL2 file, for example

```text
target_charged.mol2
```

containing the charged-state atomic charges.

For a simple \(+1\) radical cation, such a file may be generated separately using Antechamber/AM1-BCC in the `amber` environment, e.g. conceptually:

```bash
conda run -n amber antechamber \
    -i target.pdb \
    -fi pdb \
    -o target_charged.mol2 \
    -fo mol2 \
    -at gaff2 \
    -c bcc \
    -nc 1 \
    -m 2 \
    -rn MOL
```

However, the topology-generation program should **not depend on Antechamber being run internally**.

Its main input should simply be an already prepared charged-state MOL2 file.

Only the charges from this MOL2 should be transferred to the existing Amber topology. Do not transfer GAFF2 atom types or other parameters from it.

---

# Part 1: Generate the charged Amber topology

Implement something like

```text
prepare_charged_topology.py
```

This script runs in the `amber` environment.

Example interface:

```bash
conda run -n amber python prepare_charged_topology.py \
    --prmtop system.prmtop \
    --charged-mol2 target_charged.mol2 \
    --resid 25 \
    --output system_charged.prmtop
```

It should use ParmEd.

## Required behavior

Load

```text
system.prmtop
```

and select exactly one residue/molecule specified by the user.

Load

```text
target_charged.mol2
```

and transfer only

```python
atom.charge
```

from the MOL2 atoms to the selected residue.

Do NOT transfer:

```python
atom.type
```

or any other parameter.

Then write

```text
system_charged.prmtop
```

using ParmEd.

---

# Atom mapping

Do not silently assume atom mapping is correct.

At minimum validate:

1. same number of atoms;
2. matching atom names;
3. matching elements/atomic numbers if available;
4. unique mapping;
5. no target atoms are missing;
6. no additional atoms appear in the charged MOL2.

Prefer mapping by atom name when the names are unique.

If atom names are not unique, fail clearly rather than guessing.

Print a table such as

```text
index   atom   q_neutral       q_charged       delta_q
1234    C1     -0.1234          -0.0521          0.0713
1235    C2      0.1045           0.1520          0.0475
...
```

before writing the topology.

---

# Mandatory topology validation

After generating `system_charged.prmtop`, verify:

### Target molecular charge

Calculate

\[
Q_M^{(0)}=\sum_{i\in M}q_i^{(0)}
\]

and

\[
Q_M^{(+)}=\sum_{i\in M}q_i^{(+)}.
\]

Print both values and

\[
\Delta Q_M=Q_M^{(+)}-Q_M^{(0)}.
\]

For a neutral-to-cation transformation, this should be approximately

\[
\Delta Q_M=+1.
\]

Do not hard-code +1 because I may later use another charge state.

### Non-target charges

Verify programmatically that

\[
q_j^{(+)}=q_j^{(0)}
\]

for every atom outside the selected molecule.

### Topology consistency

Verify that neutral and charged systems have exactly the same:

- number of atoms;
- atom ordering;
- residue ordering;
- atom types;
- masses;
- bond definitions;
- angle definitions;
- torsion definitions;
- Lennard-Jones parameters.

The only intended difference is the target molecule's atomic charges.

If practical, implement a topology comparison utility or validation function and report something like:

```text
Topology validation
-------------------
Atoms:                    identical
Residues:                 identical
Atom types:               identical
Bond parameters:          identical
Angle parameters:         identical
Dihedral parameters:      identical
LJ parameters:            identical
Non-target charges:       identical
Target charges:           changed
Target total charge:      0.0000 -> 1.0000
Status:                   PASS
```

Fail loudly if unexpected differences are found.

---

# Part 2: Calculate the energy-gap trajectory

Implement a second script:

```text
compute_energy_gap.py
```

This script must run in the `myash` environment.

Example:

```bash
conda run -n myash python compute_energy_gap.py \
    --neutral-prmtop system.prmtop \
    --charged-prmtop system_charged.prmtop \
    --trajectory trajectory.nc \
    --output energy_gap.csv
```

Allow common trajectory formats if the installed trajectory library supports them, especially:

- Amber NetCDF
- DCD
- XTC

Inspect the packages already installed in `myash`.

Prefer MDTraj if already installed. MDAnalysis is also acceptable.

Do not add an unnecessary trajectory dependency if the existing environment already contains a suitable reader.

---

# OpenMM system construction

Construct two OpenMM systems:

```text
neutral_system
charged_system
```

using the two Amber topology files.

The OpenMM settings must be identical between them.

If an existing OpenMM MD script exists in the repository, inspect it and reproduce its system-construction settings, especially:

- nonbonded method;
- cutoff;
- PME settings;
- periodic boundary conditions;
- constraints;
- rigid water;
- Ewald error tolerance.

Do not guess these settings if the original simulation setup is available.

If no simulation setup is available, expose the relevant choices as command-line options and document the assumptions.

For a periodic system, ensure the periodic box vectors are handled correctly for every frame. This is particularly important for NPT trajectories with variable box vectors.

---

# Energy evaluation

Create separate OpenMM `Context`s for neutral and charged systems.

For each trajectory frame \(\mathbf R_t\):

1. set identical coordinates in both contexts;
2. set identical periodic box vectors;
3. evaluate

\[
U_0(t)
=
U_{\mathrm{neutral}}(\mathbf R_t);
\]

4. evaluate

\[
U_+(t)
=
U_{\mathrm{charged}}(\mathbf R_t);
\]

5. calculate

\[
\Delta E(t)=U_+(t)-U_0(t).
\]

Do NOT:

- minimize the snapshot;
- propagate dynamics;
- modify coordinates;
- equilibrate the charged system;
- generate a charged-state trajectory.

This is a vertical single-point calculation on the existing neutral trajectory.

---

# Numerical precision

Because

\[
U_+(t)-U_0(t)
\]

is obtained by subtracting two potentially large total MM energies, numerical precision matters.

Make the OpenMM platform configurable, for example

```bash
--platform CPU
```

or

```bash
--platform CUDA
```

where available.

For validation, compare several frames using a high-precision setting, preferably CPU or CUDA double precision if supported.

Do not silently rely on low precision if it produces noticeable noise in the energy difference.

Document the OpenMM platform and precision used in the output metadata/log.

---

# Output

Write a CSV file containing at least:

```text
frame
time_ps
E_neutral_kJmol
E_charged_kJmol
gap_kJmol
gap_fluctuation_kJmol
gap_cm-1
gap_fluctuation_cm-1
```

where

\[
\mathrm{gap}(t)
=
E_{\mathrm{charged}}(t)
-
E_{\mathrm{neutral}}(t)
\]

and

\[
\mathrm{gap\_fluctuation}(t)
=
\mathrm{gap}(t)
-
\langle\mathrm{gap}\rangle.
\]

Also print summary statistics:

```text
Number of frames
Mean energy gap
Standard deviation
Minimum
Maximum
```

Use OpenMM units internally wherever possible rather than hard-coded unit-conversion constants.

Saving an additional NumPy `.npz` file is welcome but not required.

---

# Frame selection

Support basic trajectory selection options:

```text
--start
--stop
--stride
```

so that I can first test the calculation on a few frames.

Example:

```bash
conda run -n myash python compute_energy_gap.py \
    --neutral-prmtop system.prmtop \
    --charged-prmtop system_charged.prmtop \
    --trajectory trajectory.nc \
    --start 0 \
    --stop 20 \
    --stride 1 \
    --output energy_gap_test.csv
```

---

# Important validation tests

Implement automated or easily runnable tests for the following.

## Test 1: neutral vs neutral

Run the energy-gap program using the same topology for both states:

```text
neutral.prmtop
neutral.prmtop
```

The result must satisfy

\[
\Delta E(t)\approx0
\]

for every frame to numerical precision.

This is the most important basic sanity check.

## Test 2: topology charge difference

Check that only target-molecule partial charges differ between the two topology files.

## Test 3: total molecular charge

Confirm that the target molecular charge changes by the expected amount.

## Test 4: atom ordering

Confirm that atom ordering and atom count in

```text
neutral.prmtop
charged.prmtop
trajectory
```

are identical/compatible.

## Test 5: short trajectory

Run the full neutral-vs-charged calculation on approximately 5-20 frames and verify that all energies and gaps are finite.

## Test 6: centered fluctuation

Numerically verify

\[
\left\langle\delta\Delta E\right\rangle \approx 0.
\]

---

# Suggested project structure

Keep the implementation simple, for example:

```text
energy_gap/
├── prepare_charged_topology.py
├── compute_energy_gap.py
├── validate_topologies.py
├── README.md
└── tests/
    ├── test_topology.py
    └── test_energy_gap.py
```

Avoid unnecessary abstraction.

The scripts should also be usable directly from the command line.

---

# README

Write a concise README explaining the complete workflow.

For example:

## Step 1: prepare charged-state charges

Run Antechamber in

```text
amber
```

if needed.

## Step 2: modify the topology

```bash
conda run -n amber python prepare_charged_topology.py ...
```

## Step 3: validate topology

```bash
conda run -n amber python validate_topologies.py ...
```

## Step 4: test several trajectory frames

```bash
conda run -n myash python compute_energy_gap.py \
    ... \
    --stop 10
```

## Step 5: run the full trajectory

```bash
conda run -n myash python compute_energy_gap.py ...
```

Clearly state that this workflow computes

\[
\Delta E(t)
=
U_{\mathrm{charged}}(\mathbf R_t)
-
U_{\mathrm{neutral}}(\mathbf R_t)
\]

on the original neutral trajectory.

---

# Scope

For now, implement **only the energy-gap trajectory and its mean-centered fluctuation**.

Do NOT implement:

- spectral density;
- Fourier transforms;
- quantum correction factors;
- charged-state MD;
- Marcus free-energy calculations;
- reorganization energies;
- QM/MM calculations.

These can be added later.

---

# Implementation approach

Before modifying or creating code:

1. inspect the current repository;
2. inspect relevant existing topology/trajectory handling scripts;
3. inspect the `amber` and `myash` environments;
4. identify the actual trajectory format;
5. identify how the original OpenMM system was constructed, if that information is available.

If Superpowers is available, follow its normal development workflow. Keep the implementation focused and test the core physics assumptions explicitly.

The most important requirements are:

\[
\boxed{
\text{charged topology}
=
\text{neutral topology with only target-molecule charges changed}
}
\]

and

\[
\boxed{
\Delta E(t)
=
U_{\mathrm{charged}}(\mathbf R_t)
-
U_{\mathrm{neutral}}(\mathbf R_t)
}
\]

using exactly the same coordinates and box for the two energy evaluations.