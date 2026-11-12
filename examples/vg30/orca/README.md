# VG30 ORCA optimization and TDDFT

This workflow treats VG30 as an isolated closed-shell cation. It does not use
MM point charges.

The calculations are separated into two stages:

1. `opt/run.sh` performs a vacuum B3LYP-D3BJ/def2-SVP geometry optimization.
2. `excited/run.sh` reads `opt/vg30_opt.xyz` and calculates ten singlet states
   using full TDDFT at the same level (not TDA).

RIJCOSX and the def2/J auxiliary basis are used to accelerate the hybrid DFT
calculations. By default, each calculation uses 8 processes and 3000 MB per
process. Override these values with `ORCA_NPROCS` and `ORCA_MAXCORE_MB`.

Run the optimization first:

```bash
ORCA_NPROCS=8 ORCA_MAXCORE_MB=3000 \
  ./examples/vg30/orca/opt/run.sh /path/to/orca
```

After `opt/vg30_opt.xyz` has been produced, run TDDFT separately:

```bash
ORCA_NPROCS=8 ORCA_MAXCORE_MB=3000 \
  ./examples/vg30/orca/excited/run.sh /path/to/orca
```

The ORCA path can instead be set once with `ORCA_BINARY`. Use the full path to
the ORCA driver when running parallel calculations.

For a detached `screen` optimization:

```bash
screen -dmS vg30-orca-opt bash -lc '
  cd /path/to/q4mdkit
  ORCA_NPROCS=8 ORCA_MAXCORE_MB=3000 \
    ./examples/vg30/orca/opt/run.sh /path/to/orca
'
```

Replace `opt` with `excited` only after the optimization has converged.
