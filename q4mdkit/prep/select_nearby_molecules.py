#!/usr/bin/env python3
"""
Select molecules near the QM region from a GRO or PDB file.
Output atom indices and molecule indices of surrounding molecules.
"""

import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union

# Import common functions from select_qmatoms
try:
    from .select_qmatoms import (
        read_structure_file,
        calculate_molecule_centers,
        write_atom_list
    )
except ImportError:
    from select_qmatoms import (
        read_structure_file,
        calculate_molecule_centers,
        write_atom_list
    )


def find_qm_molecules(molecules: List[Dict], qm_atom_indices: List[int]) -> List[Dict]:
    """
    Find molecules that contain QM atoms.
    
    Args:
        molecules: list of molecule dictionaries
        qm_atom_indices: list of atom indices in the QM region
    
    Returns:
        List of molecules that contain QM atoms
    """
    qm_set = set(qm_atom_indices)
    qm_molecules = []
    
    for mol in molecules:
        mol_atoms = set(mol['atom_indices'])
        if mol_atoms & qm_set:  # If there's any overlap
            qm_molecules.append(mol)
    
    return qm_molecules


def calculate_molecule_distance(mol1: Dict, mol2: Dict, box: np.ndarray = None) -> float:
    """
    Calculate minimum distance between two molecules.
    Uses center-to-center distance by default.
    
    Args:
        mol1: first molecule dictionary
        mol2: second molecule dictionary
        box: optional box dimensions for PBC (nm)
    
    Returns:
        Distance in nm
    """
    diff = mol1['center'] - mol2['center']
    
    # Apply minimum image convention if box is provided
    if box is not None and np.all(box > 0):
        diff = diff - box * np.round(diff / box)
    
    return np.linalg.norm(diff)


def find_nearby_molecules(molecules: List[Dict], qm_molecules: List[Dict], 
                          cutoff: float, box: np.ndarray = None) -> List[Dict]:
    """
    Find molecules within cutoff distance from QM molecules.
    
    Args:
        molecules: all molecules
        qm_molecules: molecules in QM region
        cutoff: distance cutoff in nm
        box: box dimensions for PBC (nm)
    
    Returns:
        List of nearby molecules (excluding QM molecules themselves)
    """
    qm_mol_ids = {mol['id'] for mol in qm_molecules}
    nearby = []
    
    for mol in molecules:
        # Skip if this is a QM molecule
        if mol['id'] in qm_mol_ids:
            continue
        
        # Check distance to any QM molecule
        for qm_mol in qm_molecules:
            dist = calculate_molecule_distance(mol, qm_mol, box)
            if dist <= cutoff:
                mol['min_distance_to_qm'] = dist
                nearby.append(mol)
                break
    
    # Sort by distance
    nearby.sort(key=lambda x: x.get('min_distance_to_qm', float('inf')))
    
    return nearby


def find_nearby_molecules_by_shells(molecules: List[Dict], qm_molecules: List[Dict],
                                    n_shells: int, box: np.ndarray = None) -> Dict[int, List[Dict]]:
    """
    Find molecules by coordination shells around QM molecules.
    
    Args:
        molecules: all molecules
        qm_molecules: molecules in QM region
        n_shells: number of coordination shells to find
        box: box dimensions for PBC (nm)
    
    Returns:
        Dictionary mapping shell number (1-indexed) to list of molecules in that shell
    """
    qm_mol_ids = {mol['id'] for mol in qm_molecules}
    
    # Calculate distances from each molecule to nearest QM molecule
    for mol in molecules:
        if mol['id'] in qm_mol_ids:
            mol['min_distance_to_qm'] = 0.0
            continue
        
        min_dist = float('inf')
        for qm_mol in qm_molecules:
            dist = calculate_molecule_distance(mol, qm_mol, box)
            min_dist = min(min_dist, dist)
        mol['min_distance_to_qm'] = min_dist
    
    # Sort non-QM molecules by distance
    non_qm = [m for m in molecules if m['id'] not in qm_mol_ids]
    non_qm.sort(key=lambda x: x['min_distance_to_qm'])
    
    # Group into shells
    # First, estimate typical nearest-neighbor distance
    if len(non_qm) > 0:
        first_dist = non_qm[0]['min_distance_to_qm']
    else:
        return {}
    
    shells = {}
    current_shell = 1
    shell_base_distance = first_dist
    shell_tolerance = first_dist * 0.3  # 30% tolerance for same shell
    
    for mol in non_qm:
        if current_shell > n_shells:
            break
        
        dist = mol['min_distance_to_qm']
        
        # Check if this molecule belongs to current shell or a new one
        if dist > shell_base_distance + shell_tolerance:
            current_shell += 1
            shell_base_distance = dist
            shell_tolerance = first_dist * 0.3
        
        if current_shell <= n_shells:
            if current_shell not in shells:
                shells[current_shell] = []
            shells[current_shell].append(mol)
    
    return shells


def read_atom_indices(filepath: str) -> List[int]:
    """
    Read atom indices from a file.
    Supports space-separated or newline-separated formats.
    
    Args:
        filepath: path to file containing atom indices
    
    Returns:
        List of atom indices
    """
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Parse indices (supports both space and newline separated)
    indices = []
    for token in content.split():
        try:
            indices.append(int(token))
        except ValueError:
            continue
    
    return indices


def write_molecule_indices(nearby_molecules: List[Dict], filename: str):
    """
    Write atom indices per molecule to a file.
    Each line contains: mol_id atom_index1 atom_index2 ...
    
    Args:
        nearby_molecules: list of molecule dictionaries with 'id' and 'atom_indices'
        filename: output file path
    """
    with open(filename, 'w') as f:
        for mol in sorted(nearby_molecules, key=lambda x: x['id']):
            atom_indices_str = ' '.join(map(str, sorted(mol['atom_indices'])))
            f.write(f"{mol['id']} {atom_indices_str}\n")


def select_nearby_molecules(structure_file: str, 
                            qm_atoms_file: str = None,
                            qm_atom_indices: List[int] = None,
                            cutoff: float = 1.0,
                            n_shells: int = None,
                            output: str = "nearby_atoms",
                            mol_output: str = "nearby_molecules",
                            use_pbc: bool = True,
                            verbose: bool = False) -> Dict:
    """
    Select molecules near the QM region.
    
    Args:
        structure_file: Input GRO or PDB file path
        qm_atoms_file: File containing QM atom indices (space or newline separated)
        qm_atom_indices: Direct list of QM atom indices (alternative to file)
        cutoff: Distance cutoff in nm (default: 1.0)
        n_shells: Number of coordination shells to select (alternative to cutoff)
        output: Output file for nearby atom indices
        mol_output: Output file for per-molecule atom indices (each line: mol_id atoms...)
        use_pbc: Apply periodic boundary conditions
        verbose: Show detailed information
    
    Returns:
        dict with 'nearby_atoms', 'nearby_molecules', 'qm_molecules'
    """
    print("=" * 70)
    print("Nearby Molecule Selection: Finding molecules around QM region")
    print("=" * 70)
    
    # Read QM atom indices
    if qm_atom_indices is not None:
        qm_atoms = list(qm_atom_indices)
        print(f"\nUsing provided QM atom indices: {len(qm_atoms)} atoms")
    elif qm_atoms_file is not None:
        qm_atoms = read_atom_indices(qm_atoms_file)
        print(f"\nReading QM atoms from: {qm_atoms_file}")
        print(f"  Found {len(qm_atoms)} atoms")
    else:
        raise ValueError("Either qm_atoms_file or qm_atom_indices must be provided")
    
    # Read structure file
    print(f"\nReading structure: {structure_file}")
    coords, residues, box, natoms, title = read_structure_file(structure_file)
    
    print(f"  Title: {title}")
    print(f"  Total atoms: {natoms}")
    print(f"  Number of molecules: {len(np.unique(residues))}")
    print(f"  Box size: {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} nm")
    
    if not use_pbc:
        box = None
        print("  PBC: disabled")
    else:
        print("  PBC: enabled")
    
    # Calculate molecule centers
    print(f"\nCalculating molecular centers...")
    molecules = calculate_molecule_centers(coords, residues)
    atoms_per_molecule = len(molecules[0]['atom_indices'])
    print(f"  Atoms per molecule: {atoms_per_molecule}")
    
    # Find QM molecules
    print(f"\nIdentifying QM molecules...")
    qm_molecules = find_qm_molecules(molecules, qm_atoms)
    print(f"  Found {len(qm_molecules)} QM molecules")
    print(f"  QM molecule IDs: {[m['id'] for m in qm_molecules]}")
    
    # Find nearby molecules
    if n_shells is not None:
        print(f"\nFinding molecules in {n_shells} coordination shell(s)...")
        shells = find_nearby_molecules_by_shells(molecules, qm_molecules, n_shells, box)
        
        nearby_molecules = []
        for shell_num in sorted(shells.keys()):
            shell_mols = shells[shell_num]
            print(f"\n  Shell {shell_num}: {len(shell_mols)} molecules")
            if verbose:
                for mol in shell_mols:
                    print(f"    Mol {mol['id']}: {mol['min_distance_to_qm']:.3f} nm from QM")
            nearby_molecules.extend(shell_mols)
    else:
        print(f"\nFinding molecules within {cutoff:.2f} nm of QM region...")
        nearby_molecules = find_nearby_molecules(molecules, qm_molecules, cutoff, box)
        print(f"  Found {len(nearby_molecules)} nearby molecules")
    
    # Display results
    if len(nearby_molecules) > 0:
        print(f"\nNearby molecules:")
        print(f"  {'Mol ID':<8} {'Distance (nm)':<15} {'Atom range':<15}")
        print("  " + "-" * 40)
        
        for mol in nearby_molecules[:10]:  # Show first 10
            atom_min = min(mol['atom_indices'])
            atom_max = max(mol['atom_indices'])
            dist_str = f"{mol.get('min_distance_to_qm', 0):.3f}"
            print(f"  {mol['id']:<8} {dist_str:<15} {atom_min:4d}-{atom_max:4d}")
        
        if len(nearby_molecules) > 10:
            print(f"  ... and {len(nearby_molecules) - 10} more molecules")
    
    # Collect all nearby atom indices
    all_nearby_atoms = []
    all_nearby_mol_ids = []
    for mol in nearby_molecules:
        all_nearby_atoms.extend(mol['atom_indices'])
        all_nearby_mol_ids.append(mol['id'])
    all_nearby_atoms.sort()
    all_nearby_mol_ids.sort()
    
    print(f"\nSummary:")
    print(f"  QM molecules: {len(qm_molecules)}")
    print(f"  Nearby molecules: {len(nearby_molecules)}")
    print(f"  Total nearby atoms: {len(all_nearby_atoms)}")
    
    if len(all_nearby_atoms) > 0:
        print(f"  Nearby atom range: {min(all_nearby_atoms)} - {max(all_nearby_atoms)}")
    
    # Write output files
    if len(all_nearby_atoms) > 0:
        write_atom_list(all_nearby_atoms, output)
        write_molecule_indices(nearby_molecules, mol_output)
        
        print(f"\nOutput files:")
        print(f"  {output} - atom indices of nearby molecules ({len(all_nearby_atoms)} atoms)")
        print(f"  {mol_output} - per-molecule atom indices ({len(nearby_molecules)} molecules)")
    else:
        print(f"\nNo nearby molecules found. Check cutoff distance or QM atom selection.")
    
    # Verbose output
    if verbose:
        print(f"\nDetailed information:")
        print(f"\n  QM molecules:")
        for mol in qm_molecules:
            center_str = f"({mol['center'][0]:.3f}, {mol['center'][1]:.3f}, {mol['center'][2]:.3f})"
            print(f"    Mol {mol['id']}: center = {center_str}")
        
        print(f"\n  All nearby molecules:")
        for mol in nearby_molecules:
            center_str = f"({mol['center'][0]:.3f}, {mol['center'][1]:.3f}, {mol['center'][2]:.3f})"
            print(f"    Mol {mol['id']}: center = {center_str}, dist = {mol.get('min_distance_to_qm', 0):.3f} nm")
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)
    
    return {
        'nearby_atoms': all_nearby_atoms,
        'nearby_molecules': all_nearby_mol_ids,
        'qm_molecules': [m['id'] for m in qm_molecules],
        'nearby_molecule_data': nearby_molecules
    }


def get_combined_region(structure_file: str,
                        qm_atoms_file: str = None,
                        qm_atom_indices: List[int] = None,
                        cutoff: float = 1.0,
                        output: str = "qm_and_nearby_atoms",
                        use_pbc: bool = True) -> Dict:
    """
    Get combined atom indices of QM region and nearby molecules.
    
    Args:
        structure_file: Input GRO or PDB file path
        qm_atoms_file: File containing QM atom indices
        qm_atom_indices: Direct list of QM atom indices
        cutoff: Distance cutoff in nm
        output: Output file for combined atom indices
        use_pbc: Apply periodic boundary conditions
    
    Returns:
        dict with 'combined_atoms', 'qm_atoms', 'nearby_atoms'
    """
    # Read QM atoms
    if qm_atom_indices is not None:
        qm_atoms = list(qm_atom_indices)
    elif qm_atoms_file is not None:
        qm_atoms = read_atom_indices(qm_atoms_file)
    else:
        raise ValueError("Either qm_atoms_file or qm_atom_indices must be provided")
    
    # Get nearby molecules
    result = select_nearby_molecules(
        structure_file=structure_file,
        qm_atom_indices=qm_atoms,
        cutoff=cutoff,
        use_pbc=use_pbc,
        verbose=False
    )
    
    # Combine QM and nearby atoms
    combined = sorted(set(qm_atoms) | set(result['nearby_atoms']))
    
    # Write output
    write_atom_list(combined, output)
    print(f"\nCombined region written to: {output}")
    print(f"  QM atoms: {len(qm_atoms)}")
    print(f"  Nearby atoms: {len(result['nearby_atoms'])}")
    print(f"  Total: {len(combined)}")
    
    return {
        'combined_atoms': combined,
        'qm_atoms': qm_atoms,
        'nearby_atoms': result['nearby_atoms']
    }


def generate_vmd_script(structure_file: str,
                        qm_atoms: List[int],
                        nearby_atoms: List[int],
                        output_script: str,
                        show_all_system: bool = False) -> str:
    """
    Generate VMD script to visualize QM and nearby molecules.
    
    Args:
        structure_file: Input PDB or GRO file path
        qm_atoms: List of QM atom indices (0-based)
        nearby_atoms: List of nearby atom indices (0-based)
        output_script: Output VMD script file path
        show_all_system: If True, show entire system as background
    
    Returns:
        Path to the generated script
    """
    script_dir = Path(output_script).parent
    viewpoint_file = script_dir / "viewpoint_nearby.tcl"
    
    # Convert 0-based to 1-based for VMD (VMD uses 0-based index)
    # Actually VMD uses 0-based index, so no conversion needed
    qm_atoms_str = ' '.join(map(str, sorted(qm_atoms)))
    nearby_atoms_str = ' '.join(map(str, sorted(nearby_atoms)))
    
    # Combined atoms for centering
    all_selected = sorted(set(qm_atoms) | set(nearby_atoms))
    all_selected_str = ' '.join(map(str, all_selected))
    
    with open(output_script, 'w') as f:
        f.write(f"""# VMD visualization script for QM and nearby molecules
# Generated by select_nearby_molecules.py
# Usage: vmd -e {output_script}

# Load structure
mol new {structure_file} type pdb waitfor all

# Delete default representation
mol delrep 0 top

""")
        
        if show_all_system:
            f.write("""# Background: entire system as thin lines
mol representation Lines 0.5
mol color ColorID 8
mol selection {all}
mol material Transparent
mol addrep top

""")
        
        f.write(f"""# Nearby molecules: Licorice representation (gray)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 2
mol selection {{index {nearby_atoms_str}}}
mol material Opaque
mol addrep top

# QM molecule: CPK representation (colored by element)
mol representation CPK 1.0 0.3 12.0 12.0
mol color Name
mol selection {{index {qm_atoms_str}}}
mol material Opaque
mol addrep top

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on selected molecules
set sel [atomselect top "index {all_selected_str}"]
set center [measure center $sel]
molinfo top set center [list $center]

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {{rotate_matrix}} {{[molinfo top get rotate_matrix]}}"
    puts $fp "molinfo top set {{center_matrix}} {{[molinfo top get center_matrix]}}"
    puts $fp "molinfo top set {{scale_matrix}} {{[molinfo top get scale_matrix]}}"
    puts $fp "molinfo top set {{global_matrix}} {{[molinfo top get global_matrix]}}"
    close $fp
    puts "Viewpoint saved to: $filename"
}}

proc load_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    if {{[file exists $filename]}} {{
        source $filename
        puts "Viewpoint loaded from: $filename"
    }} else {{
        puts "Viewpoint file not found: $filename"
    }}
}}

# Auto-load viewpoint if exists
if {{[file exists "{viewpoint_file}"]}} {{
    load_viewpoint
    puts "Auto-loaded viewpoint from {viewpoint_file}"
}}

puts ""
puts "QM molecule: CPK (colored by element) - {len(qm_atoms)} atoms"
puts "Nearby molecules: Licorice (gray) - {len(nearby_atoms)} atoms"
puts ""
puts "Toggle representations:"
""")
        
        rep_idx = 0
        if show_all_system:
            f.write(f'puts "  mol showrep top {rep_idx} off/on  - Background system"\n')
            rep_idx += 1
        f.write(f'puts "  mol showrep top {rep_idx} off/on  - Nearby molecules"\n')
        rep_idx += 1
        f.write(f'puts "  mol showrep top {rep_idx} off/on  - QM molecule"\n')
        
        f.write("""puts ""
puts "Viewpoint commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
""")
    
    print(f"✓ VMD script created: {output_script}")
    print(f"  Usage: vmd -e {output_script}")
    
    return output_script


