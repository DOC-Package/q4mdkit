from ash import *
import os
import sys
from qm4d4crystal.qmmm.qmmm_config import get_config

# Load QM/MM configuration
qmmm_config = get_config("qmmm_settings.yaml")

# Input from Step 1 (MM-optimized structure)
input_pdb = "../input/spiro-ometad.pdb"

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Load fragment from MM-optimized structure
frag = Fragment(pdbfile=input_pdb)
print(f"\nLoaded MM-optimized structure: {frag.numatoms} atoms")

# Load QM atom indices
qmatoms = qmmm_config.load_qmatoms()
print(f"QM region: {len(qmatoms)} atoms")

# Create theories
omm = qmmm_config.create_openmm_theory()
qm_dftb = qmmm_config.create_dftb_theory()
qmmm = qmmm_config.create_qmmm_theory(frag, qmatoms, omm, qm_dftb)

# Print system info
qmmm_config.print_system_info(frag, qmatoms)

# Run QM/MM geometry optimization
# Using LooseTZ for faster initial convergence, can tighten later if needed
Optimizer(
        fragment=frag, 
        theory=qmmm, 
        ActiveRegion=True, 
        actatoms=qmatoms, 
        maxiter=200,
        mult=6
)

# Save optimized structure
output_pdb = f"{output_dir}/spiro-ometad_opt_raw.pdb"
frag.write_pdbfile(output_pdb)
print(f"\nSaved QM/MM-optimized structure: {output_pdb}")

# Also save XYZ for inspection
output_xyz = f"{output_dir}/spiro-ometad_opt.xyz"
frag.write_xyzfile(output_xyz)
print(f"Saved XYZ file: {output_xyz}")

# Save PDB with original format (only QM region coordinates replaced)
output_pdb_orig_format = f"{output_dir}/spiro-ometad_opt.pdb"

def update_pdb_qm_coords(original_pdb, optimized_frag, qm_indices, output_pdb):
    """
    Update only the QM region coordinates in the original PDB file.
    Preserves all formatting, residue names, chain IDs, etc.
    
    Parameters
    ----------
    original_pdb : str
        Path to the original PDB file
    optimized_frag : Fragment
        ASH Fragment with optimized coordinates
    qm_indices : list
        List of QM atom indices (0-based)
    output_pdb : str
        Path to output PDB file
    """
    # Get optimized coordinates (convert to Angstrom if needed)
    opt_coords = optimized_frag.coords  # Should be in Angstrom
    
    # Create set for fast lookup (convert to 1-based for PDB atom numbering)
    qm_set = set(qm_indices)
    
    with open(original_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                # PDB format: columns 7-11 are atom serial number (1-based)
                atom_serial = int(line[6:11].strip())
                atom_idx = atom_serial - 1  # Convert to 0-based
                
                if atom_idx in qm_set:
                    # Replace coordinates (columns 31-54: x, y, z each 8.3f)
                    x, y, z = opt_coords[atom_idx]
                    new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                    fout.write(new_line)
                else:
                    fout.write(line)
            else:
                fout.write(line)

update_pdb_qm_coords(input_pdb, frag, qmatoms, output_pdb_orig_format)
print(f"Saved PDB with original format (QM coords updated): {output_pdb_orig_format}")

