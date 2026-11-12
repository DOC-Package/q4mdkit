# MM and QM/MM preparation

`MMConfig.validate_inputs()` and `QMMMConfig.validate_inputs()` check the
configured PDB and AMBER topology without importing ASH or running a QM
calculation. They return the atom count. QM/MM validation also checks QM
indices and the existence of supplied DFTB Slater–Koster files.

System-specific checks are optional YAML settings:

```yaml
validation:
  expected_atoms: 104
  residue_blocks:
    - name: DYE
      count: 4
    - name: SOL
      count: 100
  expected_qm_atoms: 4  # QM/MM only
  qm_resname: DYE      # QM/MM only
```

`residue_blocks` describes ordered blocks of **atoms**, not numbers of
molecules. Omit these checks for unrestricted compositions. QM selections
can be noncontiguous. Periodic inputs require CRYST1; when `use_boxfile` is
enabled, the box-vector lengths must match CRYST1 within 0.001 angstrom.
Box files accept three Cartesian vectors or one `a b c alpha beta gamma` row.
Paths follow the existing config convention: relative to the working directory.
For direct use, import `validate_amber_system` from `q4mdkit.validation`.

## Simple run scripts

Keep theory construction and the MD call explicit; use `prepare_md_run` for
the common validation, output-directory and restart handling:

```python
from ash import Fragment
from q4mdkit.md import prepare_md_run
from q4mdkit.qmmm.qmmm_config import get_config
from q4mdkit.qmmm.md_config import get_md_config

cfg = get_config("qmmm_settings.yaml")
md = get_md_config("md_settings.yaml")
state = prepare_md_run(cfg, md, "nvt", statefile="../npt/final_state.xml")
frag = Fragment(pdbfile=cfg.pdbfile)
qmatoms = cfg.load_qmatoms()
omm = cfg.create_openmm_theory(fragment=frag)
qm = cfg.create_qm_theory()
theory = cfg.create_qmmm_theory(frag, qmatoms, omm=omm, qm_theory=qm)
md.run_nvt(frag, theory, statefile=state)
md.save_final_structure(prefix="nvt")
```

The ensemble can be `nvt`, `npt`, or `nve`. Omit `statefile` to start from
the configured structure. For restart runs, positions and velocities must
match the input atom count; the full precision State XML box is passed to
OpenMM. NVT/NVE receive a local `initial_state.xml` with MonteCarlo barostat
parameters removed; the source file remains unchanged.

In command-line scripts, pass `argv=None` and import ASH **after** preparation,
as in [VG30 NVT](../examples/vg30/nvt-orca/nvt.py). This enables
`--validate-only` without ASH, and `--simulation-time PS`. The default
`argv=()` ignores process arguments for programmatic use. Set
`output.protect_existing: true` in the MD YAML to reject an existing ensemble
DCD or `OpenMM_MD_final_state.xml` before starting a separate run.

## QM bond constraints

ASH removes constraints involving QM atoms during `QMMMTheory` construction.
Set `openmm.preserve_qm_hbond_constraints: true` in the QM/MM YAML to capture
and restore existing covalent X–H constraints automatically. Here “HBonds”
has OpenMM's meaning: bonds involving hydrogen. The default is `false`,
preserving the behavior of existing configurations. This option retains
constraints already present in the OpenMM system; it does not create new bonds.

The lower-level helpers are
`q4mdkit.qmmm.constraints.capture_qm_hbond_constraints` and
`q4mdkit.qmmm.constraints.restore_constraints`.
