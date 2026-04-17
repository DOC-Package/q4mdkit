# ORCA QM/MM MD Example

This example demonstrates QM/MM molecular dynamics using ORCA as the QM backend.

## Prerequisites

- ASH environment with OpenMM
- ORCA 5.x installed
- AMBER topology files

## Setup

1. Activate the ash environment:
   ```bash
   conda activate ash
   ```

2. Edit `qmmm_settings.yaml`:
   - Set `qm_backend: "orca"`
   - Set `orca.orcadir` to your ORCA bin directory
   - Set `orca.orcasimpleinput` to desired method (e.g., `! B3LYP def2-SVP D3BJ TightSCF`)
   - Configure paths to your topology files

3. Edit `md_settings.yaml`:
   - Set timestep (0.5 fs recommended for ORCA)
   - Set simulation time and temperature

4. Prepare input files in `input/` directory:
   - `system.prmtop`: AMBER topology
   - `system.inpcrd`: AMBER coordinates
   - `qmatoms`: QM region atom indices (space-separated)
   - `actatoms`: Active region atom indices

## Running

```bash
python run_qmmm_md.py
```

## ORCA Method Options

Common ORCA method lines for QM/MM:

```yaml
# Standard DFT
orcasimpleinput: "! B3LYP def2-SVP D3BJ TightSCF"

# RI-accelerated (faster for large QM regions)
orcasimpleinput: "! RI-B3LYP def2-SVP def2/J RIJCOSX D3BJ TightSCF"

# Range-separated hybrid (better for CT states)
orcasimpleinput: "! CAM-B3LYP def2-SVP D3BJ TightSCF"

# wB97X-D3 (alternative range-separated)
orcasimpleinput: "! wB97X-D3 def2-SVP TightSCF"
```

## Output

- `output/nvt.dcd`: Trajectory file
- `output/nvt.csv`: Energy data
- `output/orca_energy.dat`: QM energies (if `log_enabled: true`)
