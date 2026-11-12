# VG30 DFTB optimization and TD-DFTB

This workflow treats VG30 as an isolated cation. It deliberately does not
read the existing `PCcharges.dat`; no MM point charges are included.

The two stages are:

1. `opt/dftb_in.hsd` optimizes the 77-atom cation with the reference CAM
   Hamiltonian and writes `opt/geo_end.xyz`.
2. `excited/dftb_in.hsd` reads that optimized geometry and computes ten singlet
   TD-DFTB/Casida states with the non-hybrid 3ob SK set. The global CAM hybrid
   and OPV-B3LYP SK set are intentionally omitted from this stage because
   DFTB+ linear response does not support global hybrids or hybrid SK files.

Run the geometry optimization from any directory with:

```bash
./examples/vg30/tddftb/opt/run.sh /path/to/dftb+
```

After it produces `opt/geo_end.xyz`, run the excited-state calculation
separately with:

```bash
./examples/vg30/tddftb/excited/run.sh /path/to/dftb+
```

For either script, the argument can be omitted if `dftb+` is already on
`PATH`. The optimization stage must be able to read the OPV-B3LYP Slater–Koster files under
`/home/takahashi/qmmm/sk/opv/b3lyp`; the TD-DFTB stage uses the installed 3ob
set under `/home/takahashi/qmmm/sk/3ob`.

With the SK files currently installed on this machine, use the development
DFTB+ build (for example
`/home/takahashi/anaconda3/envs/myash/bin/dftb+`). The packaged DFTB+ 24.1
binary rejects the legacy OPV-B3LYP hybrid metadata before starting the
calculation.
