#!/usr/bin/env python3
"""
Select two adjacent molecules for QM region from pentacene crystal.
Shows bond type information (Type 1-4) for each pair.

Pentacene crystal bond types:
  Type 1: A-B herringbone    ~4.89 Å (red)
  Type 2: A-A/B-B a-axis     ~6.28 Å (blue)
  Type 3: A-B diagonal       ~7.29 Å (green)
  Type 4: A-A/B-B b-axis     ~7.71 Å (orange)
"""

import numpy as np
import argparse
from pathlib import Path


# Pentacene crystal bond type parameters (distances in Å)
# Based on actual crystal structure analysis:
#   Type 1: ~4.73-4.89 Å (herringbone, shortest neighbor)
#   Type 2: ~5.20-6.28 Å (second nearest)
#   Type 3: ~6.28 Å (a-axis direction) - actually overlaps with Type 2
#   Type 4: ~7.71 Å (b-axis direction)
#
# Simplified classification by distance only (molecular type A/B is hard to 
# determine from cif2pdb output since molecules are reordered):
BOND_TYPES = {
    1: {'name': 'herringbone (nearest)', 'dist': 4.8, 'range': (4.5, 5.0), 'color': 'red'},
    2: {'name': 'second nearest', 'dist': 5.2, 'range': (5.0, 5.5), 'color': 'cyan'},
    3: {'name': 'a-axis', 'dist': 6.28, 'range': (6.0, 6.6), 'color': 'blue'},
    4: {'name': 'b-axis', 'dist': 7.71, 'range': (7.5, 8.0), 'color': 'orange'},
}


def read_gro_file(grofile):
    """Read a GRO file and return coordinates, residue numbers, and box information."""
    with open(grofile, 'r') as f:
        lines = f.readlines()
    
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    coords = []
    residues = []
    
    for i in range(2, 2 + natoms):
        line = lines[i]
        resnum = int(line[0:5].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        coords.append([x, y, z])
        residues.append(resnum)
    
    coords = np.array(coords)
    residues = np.array(residues)
    
    box_line = lines[2 + natoms].strip().split()
    box = np.array([float(box_line[0]), float(box_line[1]), float(box_line[2])])
    
    return coords, residues, box, natoms, title


def read_pdb_file(pdbfile):
    """Read a PDB file and return coordinates, residue numbers, and box information."""
    with open(pdbfile, 'r') as f:
        lines = f.readlines()
    
    coords = []
    residues = []
    title = "PDB file"
    box = np.array([0.0, 0.0, 0.0])
    
    for line in lines:
        if line.startswith("TITLE"):
            title = line[10:].strip()
        elif line.startswith("CRYST1"):
            a = float(line[6:15].strip())
            b = float(line[15:24].strip())
            c = float(line[24:33].strip())
            # Keep in Angstrom for consistency with coordinates
            box = np.array([a, b, c])
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            resnum = int(line[22:26].strip())
            # Keep coordinates in Angstrom (PDB native unit)
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x, y, z])
            residues.append(resnum)
    
    coords = np.array(coords)
    residues = np.array(residues)
    natoms = len(coords)
    
    return coords, residues, box, natoms, title


def read_structure_file(filepath):
    """Read a GRO or PDB file based on extension."""
    ext = Path(filepath).suffix.lower()
    if ext == '.gro':
        return read_gro_file(filepath)
    elif ext == '.pdb':
        return read_pdb_file(filepath)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .gro or .pdb")


def calculate_molecule_info(coords, residues):
    """
    Calculate molecule centers and determine molecule type (A or B).
    
    In pentacene crystal: even residue numbers = Type A, odd = Type B
    (based on how ASE reads the CIF with symmetry operations)
    """
    molecules = []
    unique_residues = np.unique(residues)
    
    for res_id in unique_residues:
        mask = residues == res_id
        mol_coords = coords[mask]
        atom_indices = np.where(mask)[0]
        center = mol_coords.mean(axis=0)
        
        # Determine molecule type: A (even index in unit cell) or B (odd)
        # residue 1,2 -> cell 0 (mol 1=A, mol 2=B)
        # residue 3,4 -> cell 1 (mol 3=A, mol 4=B), etc.
        mol_type = 'A' if (res_id - 1) % 2 == 0 else 'B'
        
        molecules.append({
            'id': int(res_id),
            'center': center,
            'atom_indices': atom_indices.tolist(),
            'type': mol_type
        })
    
    return molecules


def classify_bond(dist, mol1_type=None, mol2_type=None):
    """
    Classify the bond type between two molecules based on distance only.
    
    Args:
        dist: distance in Angstrom
        mol1_type, mol2_type: ignored (kept for compatibility)
    
    Returns:
        bond_type: int (1-4) or None if not a recognized bond
    """
    for btype, info in BOND_TYPES.items():
        dmin, dmax = info['range']
        if dmin < dist < dmax:
            return btype
    
    return None


def find_pairs_with_bond_types(molecules, coords):
    """
    Find all molecule pairs with their bond types.
    
    Returns:
        pairs: list of dicts with pair info and bond type
    """
    system_center = coords.mean(axis=0)
    
    # Calculate distance from center for each molecule
    for mol in molecules:
        mol['distance_to_center'] = np.linalg.norm(mol['center'] - system_center)
    
    pairs = []
    n_mols = len(molecules)
    
    for i in range(n_mols):
        for j in range(i + 1, n_mols):
            mol1 = molecules[i]
            mol2 = molecules[j]
            
            # Distance between molecule centers
            dist = np.linalg.norm(mol1['center'] - mol2['center'])
            
            # Classify bond type
            bond_type = classify_bond(dist, mol1['type'], mol2['type'])
            
            if bond_type is not None:
                # Distance of pair center from system center
                pair_center = (mol1['center'] + mol2['center']) / 2
                pair_dist_to_center = np.linalg.norm(pair_center - system_center)
                
                pairs.append({
                    'mol1': mol1,
                    'mol2': mol2,
                    'distance': dist,
                    'bond_type': bond_type,
                    'pair_center': pair_center,
                    'pair_dist_to_center': pair_dist_to_center,
                    'mol_ids': sorted([mol1['id'], mol2['id']])
                })
    
    return pairs, system_center


def write_atom_list(atom_indices, filename):
    """Write atom indices to a file."""
    with open(filename, 'w') as f:
        f.write(' '.join(map(str, atom_indices)) + '\n')


def select_qmatoms_pentacene(structure_file, bond_type=None, 
                              output="qmatoms", verbose=False, choice=None):
    """
    Select molecule pairs for QM region in pentacene crystal.
    
    Args:
        structure_file: Input GRO or PDB file path
        bond_type: Filter by bond type (1-4). If None, show all types.
        output: Output file for QM atom indices
        verbose: Show detailed information
        choice: If multiple candidates exist, select this index (1-based)
    
    Returns:
        dict with 'qm_atoms', 'mol_ids', 'bond_type'
    """
    print("=" * 75)
    print("Pentacene QM Region Selection: Select molecule pairs by bond type")
    print("=" * 75)
    
    # Print bond type legend
    print("\nBond types in pentacene crystal:")
    for btype, info in BOND_TYPES.items():
        print(f"  Type {btype}: {info['name']:<20} ~{info['dist']:.2f} Å ({info['color']})")
    
    # Read structure file
    print(f"\nReading: {structure_file}")
    coords, residues, box, natoms, title = read_structure_file(structure_file)
    
    print(f"  Title: {title}")
    print(f"  Total atoms: {natoms}")
    n_molecules = len(np.unique(residues))
    print(f"  Number of molecules: {n_molecules}")
    print(f"  Box size: {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} Å")
    
    # Calculate molecule centers and types
    print(f"\nAnalyzing molecular structure...")
    molecules = calculate_molecule_info(coords, residues)
    atoms_per_molecule = len(molecules[0]['atom_indices'])
    print(f"  Atoms per molecule: {atoms_per_molecule}")
    
    # Count A and B molecules
    n_type_a = sum(1 for m in molecules if m['type'] == 'A')
    n_type_b = sum(1 for m in molecules if m['type'] == 'B')
    print(f"  Type A molecules: {n_type_a}")
    print(f"  Type B molecules: {n_type_b}")
    
    # Find all pairs with bond types
    print(f"\nFinding molecule pairs...")
    pairs, system_center = find_pairs_with_bond_types(molecules, coords)
    print(f"  System center: ({system_center[0]:.2f}, {system_center[1]:.2f}, {system_center[2]:.2f}) Å")
    
    # Count pairs by bond type
    type_counts = {1: 0, 2: 0, 3: 0, 4: 0}
    for p in pairs:
        type_counts[p['bond_type']] += 1
    
    print(f"\n  Pairs by bond type:")
    for btype in sorted(type_counts.keys()):
        info = BOND_TYPES[btype]
        print(f"    Type {btype} ({info['name']}): {type_counts[btype]} pairs")
    
    # Filter by bond type if specified
    if bond_type is not None:
        pairs = [p for p in pairs if p['bond_type'] == bond_type]
        if not pairs:
            print(f"\nNo pairs found with bond type {bond_type}!")
            return None
        print(f"\nFiltered to Type {bond_type} pairs: {len(pairs)} pairs")
    
    # Sort by distance to center
    pairs.sort(key=lambda x: x['pair_dist_to_center'])
    
    # Display candidates
    print(f"\n{'='*75}")
    print(f"Candidate pairs (sorted by distance to center):")
    print(f"{'='*75}")
    print(f"  {'#':<4} {'Mol IDs':<12} {'Type':<6} {'Bond Type':<25} {'Dist (Å)':<10} {'From Center (Å)':<15}")
    print("  " + "-" * 70)
    
    max_display = 20 if verbose else 10
    for i, p in enumerate(pairs[:max_display], 1):
        mol_ids_str = f"{p['mol1']['id']}-{p['mol2']['id']}"
        types_str = f"{p['mol1']['type']}-{p['mol2']['type']}"
        bond_info = BOND_TYPES[p['bond_type']]
        bond_str = f"Type {p['bond_type']}: {bond_info['name']}"
        print(f"  {i:<4} {mol_ids_str:<12} {types_str:<6} {bond_str:<25} {p['distance']:<10.2f} {p['pair_dist_to_center']:<15.2f}")
    
    if len(pairs) > max_display:
        print(f"  ... and {len(pairs) - max_display} more pairs")
    
    print()
    
    # Select pair
    if len(pairs) == 1:
        selected_idx = 0
        print(f"Only one candidate, selecting pair {pairs[0]['mol_ids']}")
    elif choice is not None:
        selected_idx = choice - 1
        print(f"Using choice {choice} (specified by parameter)")
    else:
        try:
            user_input = input(f"Select pair (1-{min(len(pairs), max_display)}) [default: 1]: ").strip()
            selected_idx = int(user_input) - 1 if user_input else 0
        except (ValueError, EOFError):
            selected_idx = 0
            print("Using default choice: 1")
    
    if not (0 <= selected_idx < len(pairs)):
        print(f"Invalid choice, using first pair")
        selected_idx = 0
    
    selected_pair = pairs[selected_idx]
    
    # Display selected pair info
    print(f"\n{'='*75}")
    print(f"Selected pair:")
    print(f"{'='*75}")
    bond_info = BOND_TYPES[selected_pair['bond_type']]
    print(f"  Molecules: {selected_pair['mol1']['id']} ({selected_pair['mol1']['type']}) and "
          f"{selected_pair['mol2']['id']} ({selected_pair['mol2']['type']})")
    print(f"  Bond type: Type {selected_pair['bond_type']} - {bond_info['name']}")
    print(f"  Distance: {selected_pair['distance']:.3f} Å")
    print(f"  Distance from center: {selected_pair['pair_dist_to_center']:.2f} Å")
    
    # Collect atom indices
    qm_atoms1 = sorted(selected_pair['mol1']['atom_indices'])
    qm_atoms2 = sorted(selected_pair['mol2']['atom_indices'])
    all_qm_atoms = sorted(qm_atoms1 + qm_atoms2)
    
    print(f"\nQM region:")
    print(f"  Total atoms: {len(all_qm_atoms)}")
    print(f"  Molecule 1 (res {selected_pair['mol1']['id']}): {len(qm_atoms1)} atoms, index {min(qm_atoms1)} - {max(qm_atoms1)}")
    print(f"  Molecule 2 (res {selected_pair['mol2']['id']}): {len(qm_atoms2)} atoms, index {min(qm_atoms2)} - {max(qm_atoms2)}")
    print(f"  Atom index range: {min(all_qm_atoms)} - {max(all_qm_atoms)}")
    
    # Write to files
    write_atom_list(all_qm_atoms, output)
    output1 = output + "1"
    output2 = output + "2"
    write_atom_list(qm_atoms1, output1)
    write_atom_list(qm_atoms2, output2)
    print(f"\nOutput files:")
    print(f"  {output}  - atom indices for QM region ({len(all_qm_atoms)} atoms)")
    print(f"  {output1} - atom indices for molecule 1 ({len(qm_atoms1)} atoms)")
    print(f"  {output2} - atom indices for molecule 2 ({len(qm_atoms2)} atoms)")
    
    print("\n" + "=" * 75)
    print("Done!")
    print("=" * 75)
    
    return {
        'qm_atoms': all_qm_atoms,
        'qm_atoms1': qm_atoms1,
        'qm_atoms2': qm_atoms2,
        'mol_ids': selected_pair['mol_ids'],
        'bond_type': selected_pair['bond_type'],
        'bond_name': bond_info['name'],
        'distance': selected_pair['distance']  # in Å
    }


def main():
    parser = argparse.ArgumentParser(
        description='Select molecule pairs by bond type from pentacene crystal for QM region'
    )
    parser.add_argument('structure_file', help='Input GRO or PDB file')
    parser.add_argument('-t', '--bond-type', type=int, choices=[1, 2, 3, 4], default=None,
                        help='Filter by bond type (1: A-B herringbone, 2: A-A a-axis, '
                             '3: A-B diagonal, 4: A-A b-axis). If not specified, show all.')
    parser.add_argument('-o', '--output', default='qmatoms',
                        help='Output file name (default: qmatoms)')
    parser.add_argument('-c', '--choice', type=int, default=None,
                        help='Select pair by number (1-based), skip interactive prompt')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed information')
    
    args = parser.parse_args()
    
    select_qmatoms_pentacene(
        structure_file=args.structure_file,
        bond_type=args.bond_type,
        output=args.output,
        verbose=args.verbose,
        choice=args.choice
    )


if __name__ == "__main__":
    main()
