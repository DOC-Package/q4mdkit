# QM/MM Simulation with pyCHARMM

## Overview

`run_qmmm.py` implements QM/MM (Quantum Mechanics / Molecular Mechanics) molecular dynamics using pyCHARMM.

## Features

- **QM Region**: Central 2 molecules (residues 27, 28) treated with semi-empirical QM
- **MM Region**: Remaining 52 molecules treated with CHARMM36 CGenFF force field
- **QM Method**: PM7 (MOPAC) - fast semi-empirical Hamiltonian
- **Boundary**: Dividing scheme for QM/MM interface
- **Simulation**: Energy minimization + NVT equilibration with trajectory output

## System Details

- Total: 54 pentacene molecules (1944 atoms)
- QM atoms: 72 (2 molecules × 36 atoms/molecule)
- MM atoms: 1872 (52 molecules × 36 atoms/molecule)
- Periodic boundary conditions with PME electrostatics
- Triclinic unit cell

## Requirements

### Software
1. **pyCHARMM** with QM/MM support
   - CHARMM must be compiled with QM interface enabled
   - Check with: `python test_qmmm_support.py`

2. **MOPAC** (or other QM program)
   - Install MOPAC: http://openmopac.net/
   - Set environment variable: `export MOPAC_LICENSE=/path/to/license`
   - Ensure `mopac` is in PATH

### Alternative QM Programs
If MOPAC is not available, modify `QM_PROGRAM` and `QM_METHOD` in `run_qmmm.py`:

- **GAMESS**: `QM_PROGRAM = "GAMESS"`, `QM_METHOD = "HF/6-31G*"`
- **ORCA**: `QM_PROGRAM = "ORCA"`, `QM_METHOD = "HF"`
- **Gaussian**: `QM_PROGRAM = "GAUSSIAN"`, `QM_METHOD = "HF/6-31G*"`

## Usage

### Basic Run
```bash
# Activate pyCHARMM environment
source ~/anaconda3/bin/activate pycharmm

# Run QM/MM simulation
python run_qmmm.py
```

### Key Parameters (edit in script)

```python
# QM region selection
QM_RESIDUES = [27, 28]  # Residue IDs for QM treatment

# QM method
QM_METHOD = "PM7"       # MOPAC method
QM_PROGRAM = "MOPAC"    # QM program name
QM_CHARGE = 0           # Total charge of QM region
QM_MULT = 1             # Spin multiplicity (1=singlet)

# Simulation parameters
MIN_SD = 50             # Steepest descent steps
MIN_ABNR = 100          # ABNR minimization steps
NVT_STEPS = 500         # NVT dynamics steps (1 ps)
DT_PS = 0.002           # Timestep (2 fs)
TEMP_K = 300.0          # Temperature (K)
```

## Output Files

```
out/
├── qmmm_min.coor.crd           # Minimized coordinates
├── qmmm_min.coor.pdb           # Minimized structure (PDB)
├── qmmm_nvt.dcd                # NVT trajectory (binary)
├── qmmm_after_nvt.coor.crd     # Final coordinates
└── qmmm_after_nvt.coor.pdb     # Final structure (PDB)
```

## Performance Notes

QM/MM calculations are **significantly slower** than pure MM:

- **PM7 (semi-empirical)**: ~10-100× slower than MM
- **DFT methods**: ~100-1000× slower than MM
- **Ab initio methods**: ~1000-10000× slower than MM

Recommendations:
- Use semi-empirical methods (PM7, AM1) for dynamics
- Reduce simulation length for testing
- Increase `nsavc` (trajectory save frequency) to reduce I/O overhead
- Consider using fewer QM atoms if possible

## QM Region Selection

The script selects the **central 2 molecules** (residues 27-28 out of 54 total).

To modify QM region:

```python
# Select different residues
QM_RESIDUES = [20, 21, 22]  # 3 molecules

# Select by position (requires manual PSF inspection)
# Option 1: Select by atom index range
# Option 2: Select by geometric criteria (distance from center)
```

## Validation

1. **Check QM/MM setup**:
   ```bash
   python test_qmmm_support.py
   ```

2. **Compare energies**:
   - Pure MM energy: Run `run_md.py` briefly
   - QM/MM energy: Should be different (QM region treated differently)

3. **Visualize QM region**:
   ```bash
   # Load minimized structure in VMD
   vmd -pdb out/qmmm_min.coor.pdb
   
   # In VMD Tk Console:
   mol modselect 0 0 "resid 27 28"
   mol modcolor 0 0 ColorID 1  # Red for QM region
   ```

## Troubleshooting

### Issue: "QUANTUM command not recognized"
- **Solution**: CHARMM not compiled with QM support. Rebuild CHARMM with `-DQUANTUM` flag or use pre-built QM-enabled version.

### Issue: "MOPAC executable not found"
- **Solution**: Install MOPAC and add to PATH, or specify full path in CHARMM input.

### Issue: QM/MM energy very high
- **Solution**: 
  1. Check QM region charge/multiplicity
  2. Perform longer minimization
  3. Check for broken bonds at QM/MM boundary

### Issue: Very slow performance
- **Solution**:
  1. Use PM7 instead of DFT
  2. Reduce QM region size
  3. Reduce trajectory save frequency (`nsavc`)
  4. Use smaller timestep if unstable

## Advanced Usage

### Custom QM Selection
```python
# Example: Select molecules within 5 Å of center
def select_central_molecules(crd_path, cutoff=5.0):
    # Read coordinates
    # Calculate center of mass
    # Find molecules within cutoff
    # Return residue IDs
    pass

QM_RESIDUES = select_central_molecules("pentacene.crd", cutoff=5.0)
```

### Different QM Methods
```python
# MOPAC options
QM_METHOD = "PM7"      # Default, balanced
QM_METHOD = "PM6"      # Older, faster
QM_METHOD = "PM3"      # Even older
QM_METHOD = "AM1"      # Alternative parameterization

# For GAMESS (if available)
QM_PROGRAM = "GAMESS"
QM_METHOD = "RHF"      # Restricted Hartree-Fock
# or
QM_METHOD = "B3LYP"    # DFT (much slower)
```

### QM/MM Boundary Options
```python
# In setup_qmmm function:
boundary="DIV"   # Dividing scheme (default, no link atoms)
boundary="LINK"  # Link atom scheme (adds H at boundaries)
```

## References

1. CHARMM QM/MM Documentation: https://www.charmm.org/charmm/documentation/by-version/c47a1/params/doc/quantum/
2. MOPAC Manual: http://openmopac.net/Manual/
3. QM/MM methodology: Senn & Thiel, Angew. Chem. Int. Ed. 2009, 48, 1198-1229

## Example Workflow

```bash
# 1. Test QM/MM support
python test_qmmm_support.py

# 2. Run QM/MM simulation (short test)
# Edit run_qmmm.py: Set NVT_STEPS = 100
python run_qmmm.py

# 3. Visualize results
vmd out/qmmm_min.coor.pdb
# or
vmd pentacene.psf out/qmmm_nvt.dcd

# 4. Compare with pure MM
python run_md.py  # Original MM simulation
# Compare energies in output logs
```

## Notes

- QM/MM is primarily for **energy calculations and properties**, not long production runs
- For production MD, consider pure MM or use QM/MM for specific reaction coordinates
- Always validate QM region selection and energies before production runs
- Consider using implicit solvent or larger MM buffer around QM region for better accuracy
