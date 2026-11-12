# VG30 in acetonitrile: ORCA QM/MM MD

The ORCA examples use the 77-atom VG30 cation (indices 0–76, charge +1,
singlet) as the QM region and the 920 acetonitrile molecules as the AMBER MM
region, with electrostatic embedding and periodic boundaries.

The stages are `npt-mm` → `nvt-orca` → `nve-orca`. The existing `nvt` and
`npt` directories use DFTB. Isolated-molecule optimization and TDDFT examples
are documented separately in [orca/README.md](orca/README.md).

## Setup

Use an environment with q4mdkit, ASH, OpenMM, PyYAML, mdtraj, and ORCA.
Configure ORCA discovery in ASH or set `orca.orcadir` to the directory
containing the ORCA driver in **both** stages' `qmmm_settings.yaml` files.
The parallel ORCA installation must have its matching MPI runtime available.

Complete the MM preparation and NPT stages first. NVT requires:

- `input/vg30_acetonitrile.prmtop` and `input/qmatoms`
- `npt-mm/output/npt_lastframe.pdb`
- `npt-mm/OpenMM_MD_final_state.xml`, including positions and velocities

The default QM method is B3LYP-D3BJ/def2-SVP with RIJCOSX and def2/J,
8 processes, and `%maxcore 3000` (MB per process). Each stage runs for 1 ps
with a 0.5 fs timestep and writes trajectory frames every 50 steps. NVT uses
300 K LangevinMiddle dynamics; NVE uses Verlet without a thermostat or barostat.
Edit the YAML settings for your production run and keep the QM method,
constraints, hydrogen masses, and timestep consistent between stages.
X–H constraints on VG30 are retained in both stages, matching the existing
VG30 constrained QM/MM setup.

All six MM and QM/MM MD scripts use the same explicit sequence: load YAML,
prepare inputs, create the Fragment and theories, run MD, and save the final
structure. Shared preparation is provided by `q4mdkit.md.prepare_md_run`.
Input checks use `MMConfig.validate_inputs()` / `QMMMConfig.validate_inputs()`;
VG30 atom counts and residue names are set in each YAML `validation` section.
QM X–H constraints are restored by the QM/MM factory when
`openmm.preserve_qm_hbond_constraints: true` (DFTB NVT and ORCA NVT/NVE).
See [the reusable MD utilities](../../docs/md.md) for other systems.

## Run

From the repository root, validate the NVT inputs without importing ASH or
starting an ORCA calculation:

```bash
./examples/vg30/nvt-orca/run.sh --validate-only
```

Run NVT, then validate and run NVE after NVT finishes:

```bash
./examples/vg30/nvt-orca/run.sh
./examples/vg30/nve-orca/run.sh --validate-only
./examples/vg30/nve-orca/run.sh
```

The scripts work from any directory and accept `PYTHON=/path/to/python`.
They select OpenMPI from `/usr/bin` while retaining the Python selected in
your active environment. To use another matching OpenMPI installation, set
`ORCA_MPI_BINDIR=/path/to/openmpi/bin`. The launchers reject a non-OpenMPI
runtime before starting Python.
For a short smoke run, override the duration in ps:

```bash
./examples/vg30/nvt-orca/run.sh --simulation-time 0.001
```

To run in a detached screen session:

```bash
screen -dmS vg30_orca_nvt -L -Logfile vg30_orca_nvt.log \
  bash ./examples/vg30/nvt-orca/run.sh
```

## Restart and output

Each stage loads its predecessor's final PDB and State XML. The XML supplies
positions, velocities, and full precision box vectors; NVE continues the ORCA
NVT velocities. `initial_state.xml` is a local copy with MonteCarlo barostat
parameters removed. The source state is preserved.

The restart reader is also available as a library utility:

```python
from q4mdkit.qmmm.xml import read_restart_box

box_vectors = read_restart_box("OpenMM_MD_final_state.xml", expected_atoms=5597)
```

It returns vectors in angstrom and validates restart positions and velocities.
For other systems, supply their atom count or omit `expected_atoms` to infer
the count from the positions.

Each stage directory contains `OpenMM_MD_final_state.xml` and its ASH checkpoint.
The `output` directory contains `nvt.dcd`/`nve.dcd`, the corresponding CSV energy
log, first/last-frame PDB files, a separate wrapped DCD for visualization, and
`orca_energy.dat` with QM energies. ORCA working files stay in each stage
directory. Each command starts from its predecessor, rather than continuing
the same stage. Before rerunning a completed stage (including after a smoke
run), archive its `output` directory and `OpenMM_MD_final_state.xml`.
The runner refuses to start if its DCD or final State XML already exists,
preventing ASH from appending a separate run to the previous trajectory.

## MPI startup errors

The installed ORCA binary on this machine links to OpenMPI, while the shell
may select Intel MPI from `/opt/intel/oneapi/mpi/2021.14/bin/mpirun`.
This mismatch can make every worker run as rank 0 of a separate one-process
MPI world. The workers then access the same ORCA working files independently,
causing startup failures such as "the number of points read from the grid
does not match the expectation" and SHARK matrix errors.

For this installation, `/usr/bin/mpirun` is OpenMPI 4.1.6. Both `run.sh`
launchers select it automatically, resolving the active Python first so the
MPI path change does not switch to the system Python. After activating your
ASH Python environment, run:

```bash
./examples/vg30/nvt-orca/run.sh
```

For NVE, use `nve-orca/run.sh`. On other machines, set `ORCA_MPI_BINDIR`
to the MPI runtime matching the ORCA installation. See the
[ORCA parallel execution instructions](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html).

Before retrying this failed startup, archive its working files and output
from within `nvt-orca` (after the failed job has stopped):

```bash
backup_dir="failed-mpi-$(date +%Y%m%d-%H%M%S)"
mkdir "$backup_dir"
mv orca.* output "$backup_dir/"
```

This retains the failed log and clears the empty DCD that otherwise prevents
a fresh run. It also avoids reusing incomplete working files.
