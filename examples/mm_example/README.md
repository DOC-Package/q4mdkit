# MM Example

Pure force field (OpenMM) MD simulation example.

## Directory Structure

```
mm_example/
├── input/           # Input files (prmtop, inpcrd, pdb)
├── nvt/             # NVT equilibration
├── npt/             # NPT equilibration  
└── nve/             # NVE production
```

## Usage

1. Prepare input files in `input/` directory:
   - `system.prmtop` - AMBER topology file
   - `system.inpcrd` - AMBER coordinate file
   - `system.pdb` - PDB structure file

2. Run NVT equilibration:
   ```bash
   cd nvt
   python nvt.py
   ```

3. Run NPT equilibration:
   ```bash
   cd npt
   python npt.py
   ```

4. Run NVE production:
   ```bash
   cd nve
   python nve.py
   ```

## Configuration

- `mm_settings.yaml` - OpenMM force field settings
- `md_settings.yaml` - MD simulation parameters (timestep, temperature, etc.)
