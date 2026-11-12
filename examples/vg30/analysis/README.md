# VG30 trajectory snapshots: ORCA TDDFT

`run_tdorca.py` follows `run_tddftb.py`: read the trajectory, extract the
selected QM atoms at each snapshot, run the QM calculation, and append
excitation energies and oscillator strengths. The ORCA version uses
B3LYP-D3BJ/def2-SVP with RIJCOSX and full singlet TDDFT (`tda false`), VG30
charge +1, and multiplicity 1. It computes 3 states by default, matching the
TDDFTB settings. The model is editable in `input/orca_tddft.inp`.

## Run

Activate an environment containing q4mdkit, NumPy, PyYAML, and MDTraj. ASH and
DFTB+ are not required for this analysis. Set `orca.binary` in the YAML to the
full path of the ORCA driver. Then, from the repository root:

```bash
./examples/vg30/analysis/run_tdorca.sh
```

The launcher selects OpenMPI from `/usr/bin` and retains your active Python.
Override these with `ORCA_MPI_BINDIR=/path/to/openmpi/bin` and
`PYTHON=/path/to/python` as needed. It uses the same MPI setup as the ORCA MD
launchers. To run with an already configured MPI environment directly:

```bash
python examples/vg30/analysis/run_tdorca.py
```

## Trajectories and frames

Paths in `tdorca_settings.yaml` are relative to that YAML file, matching
`tddftb_settings.yaml`. The default trajectory is `../nve/output/nve.dcd`.
Point `input.trajectory` and `input.topology` at the trajectory and matching
PDB/topology you want to analyze. For example, an ORCA NVT trajectory uses
`../nvt-orca/output/nvt.dcd` and `../nvt-orca/output/nvt_firstframe.pdb`.

Use a separate settings file and `output.directory` for each trajectory:

```bash
./examples/vg30/analysis/run_tdorca.sh examples/vg30/analysis/tdorca_nvt.yaml
```

`frames.start` is the first zero-based frame; `frames.n_frames: null` processes
all remaining frames. A positive count restricts the calculation to that many
snapshots. Set `frames.dt_fs` to the saved-frame interval; `Time(fs)` is
`t0_fs + frame_index * dt_fs`. Its initial value of 4 fs follows the TDDFTB
example. For a 0.5 fs MD timestep saving every 50 steps, use **25 fs**.

QM indices come from the same zero-based `input/qmatoms`. Coordinates are
converted from MDTraj nm to angstrom. The analysis preserves the snapshots'
coordinates as in the TDDFTB workflow, so use trajectories in which the
QM molecule is whole across periodic boundaries.

## Optional solvent point charges

The reference TDDFTB configuration disables point charges. This is also the
ORCA default. To include the MM solvent, enable:

```yaml
input:
  pccharges_template: "input/PCcharges.dat"
```

The file has columns `x y z charge`, one row per non-QM atom in topology order.
Charges stay fixed and solvent coordinates are replaced for every snapshot.
The runner converts them to ORCA's `count` followed by `charge x y z` format.
A point-charge count mismatch raises an error before any QM calculation.

## Output

`tdorca/excitations.dat` has the same horizontal layout as TDDFTB:

```text
# Frame  Time(fs)  E1(eV)  f1  E2(eV)  f2  E3(eV)  f3
```

Each successful frame is appended immediately. States are energy ordered;
no state tracking or spectrum broadening is applied by the example defaults.
Failed or incomplete calculations are reported and skipped; the process
exits with an error if no snapshots succeed.

`tdorca/work` holds the current `qm_coords.xyz`, `orca.inp`, `orca.out`, and
optional `pointcharges.pc`. These working files are reused between snapshots;
ORCA can use the previous snapshot's orbitals through AutoStart. Each new
analysis run recreates the excitation table. Keep separate output directories
for independent analyses. Other users of the common tdorca module retain
its default long-format table and average-spectrum output.
