#!/usr/bin/env python3
"""
Select two adjacent molecules closest to the center of the system from a GRO or PDB file.
Output atom indices to be used for the QM region.
"""

import numpy as np
import argparse
from pathlib import Path


def read_gro_file(grofile):
    """
    Read a GRO file and return coordinates, residue numbers, and box information.
    
    Returns:
        coords: numpy array (N, 3) - atomic coordinates (nm)
        residues: numpy array (N,) - residue numbers
        box: numpy array (3,) - box dimensions (nm)
        natoms: int - total number of atoms
        title: str - title from file
    """
    with open(grofile, 'r') as f:
        lines = f.readlines()
    
    # Title and number of atoms
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    # Read coordinates and residue numbers
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
    
    # Box dimensions
    box_line = lines[2 + natoms].strip().split()
    box = np.array([float(box_line[0]), float(box_line[1]), float(box_line[2])])
    
    return coords, residues, box, natoms, title


def read_pdb_file(pdbfile):
    """
    Read a PDB file and return coordinates, residue numbers, and box information.
    
    Returns:
        coords: numpy array (N, 3) - atomic coordinates (nm, converted from Angstrom)
        residues: numpy array (N,) - residue numbers
        box: numpy array (3,) - box dimensions (nm)
        natoms: int - total number of atoms
        title: str - title from file
    """
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
            # CRYST1 record: a, b, c in Angstrom
            a = float(line[6:15].strip())
            b = float(line[15:24].strip())
            c = float(line[24:33].strip())
            # Convert to nm
            box = np.array([a / 10.0, b / 10.0, c / 10.0])
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB format: columns are fixed width
            # Residue number: columns 23-26 (1-indexed)
            # X: columns 31-38, Y: columns 39-46, Z: columns 47-54 (in Angstrom)
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip()) / 10.0  # Convert to nm
            y = float(line[38:46].strip()) / 10.0
            z = float(line[46:54].strip()) / 10.0
            coords.append([x, y, z])
            residues.append(resnum)
    
    coords = np.array(coords)
    residues = np.array(residues)
    natoms = len(coords)
    
    return coords, residues, box, natoms, title


def read_structure_file(filepath):
    """
    Read a GRO or PDB file based on extension.
    
    Returns:
        coords: numpy array (N, 3) - atomic coordinates (nm)
        residues: numpy array (N,) - residue numbers
        box: numpy array (3,) - box dimensions (nm)
        natoms: int - total number of atoms
        title: str - title from file
    """
    ext = Path(filepath).suffix.lower()
    if ext == '.gro':
        return read_gro_file(filepath)
    elif ext == '.pdb':
        return read_pdb_file(filepath)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .gro or .pdb")


def calculate_molecule_centers(coords, residues):
    """
    Calculate the center of mass for each molecule.
    
    Returns:
        molecules: list of dict with keys:
            - id: residue number
            - center: center of mass coordinates
            - atom_indices: list of atom indices
    """
    molecules = []
    unique_residues = np.unique(residues)
    
    for res_id in unique_residues:
        # Mask for atoms belonging to this residue
        mask = residues == res_id
        mol_coords = coords[mask]
        atom_indices = np.where(mask)[0]
        
        # Calculate center of mass
        center = mol_coords.mean(axis=0)
        
        molecules.append({
            'id': int(res_id),
            'center': center,
            'atom_indices': atom_indices.tolist()
        })
    
    return molecules


def find_closest_molecules(molecules, coords, n_molecules=2):
    """
    Find n molecules closest to the geometric center of all atoms.
    
    Args:
        molecules: list of molecules
        coords: all atomic coordinates (N, 3)
        n_molecules: number of molecules to select
    
    Returns:
        selected_molecules: list of selected molecules (sorted by distance)
        system_center: the geometric center of all atoms
    """
    # Use geometric center of all coordinates (not box center)
    # This is more robust for triclinic cells where atoms may extend beyond box
    system_center = coords.mean(axis=0)
    
    # Calculate distance from each molecule center to the system center
    for mol in molecules:
        distance = np.linalg.norm(mol['center'] - system_center)
        mol['distance_to_center'] = distance
    
    # Sort by distance
    molecules_sorted = sorted(molecules, key=lambda x: x['distance_to_center'])
    
    # Select the n closest molecules
    selected = molecules_sorted[:n_molecules]

    return selected, system_center, molecules_sorted


def find_molecule_combinations(molecules_sorted, n_molecules=2, distance_threshold=0.3):
    """
    Find all possible combinations of n molecules near the center within a distance threshold.
    
    Args:
        molecules_sorted: list of molecules sorted by distance to center
        n_molecules: number of molecules to select
        distance_threshold: maximum additional distance from the nth molecule (nm)
    
    Returns:
        list of candidate molecule combinations, each with total distance
    """
    from itertools import combinations
    
    if len(molecules_sorted) < n_molecules:
        return []
    
    # Find the distance of the nth closest molecule
    base_distance = molecules_sorted[n_molecules - 1]['distance_to_center']
    
    # Find all molecules within threshold of the nth molecule's distance
    candidates = []
    for mol in molecules_sorted:
        if mol['distance_to_center'] <= base_distance + distance_threshold:
            candidates.append(mol)
    
    # Generate all combinations
    all_combinations = []
    for combo in combinations(candidates, n_molecules):
        total_distance = sum(m['distance_to_center'] for m in combo)
        avg_distance = total_distance / n_molecules
        all_combinations.append({
            'molecules': list(combo),
            'total_distance': total_distance,
            'avg_distance': avg_distance,
            'mol_ids': sorted([m['id'] for m in combo])
        })
    
    # Sort by total distance
    all_combinations.sort(key=lambda x: x['total_distance'])
    
    return all_combinations


def write_atom_list(atom_indices, filename):
    """Write atom indices to a file."""
    with open(filename, 'w') as f:
        f.write(' '.join(map(str, atom_indices)) + '\n')


def select_qmatoms(structure_file, n_molecules=2, n_active_molecules=None, 
                   output="qmatoms", active_output="active_atoms", verbose=False,
                   distance_threshold=0.3, choice=None):
    """
    Select molecules near the center of the system for QM region.
    
    Args:
        structure_file: Input GRO or PDB file path
        n_molecules: Number of molecules for QM region (default: 2)
        n_active_molecules: Number of molecules for active region (default: same as n_molecules)
        output: Output file for QM atom indices (default: "qmatoms")
        active_output: Output file for active atom indices (default: "active_atoms")
        verbose: Show detailed information
        distance_threshold: Maximum additional distance for candidate search (nm, default: 0.3)
        choice: If multiple candidates exist, select this index (1-based). If None, prompt user.
    
    Returns:
        dict with 'qm_atoms' and 'active_atoms' lists
    """
    print("=" * 70)
    print("QM Region Selection: Selecting molecules near the center of the system")
    print("=" * 70)
    
    # Read structure file (GRO or PDB)
    print(f"\nReading: {structure_file}")
    coords, residues, box, natoms, title = read_structure_file(structure_file)
    
    box_center = box / 2.0
    print(f"  Title: {title}")
    print(f"  Total atoms: {natoms}")
    print(f"  Number of molecules: {len(np.unique(residues))}")
    print(f"  Box size: {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} nm")
    print(f"  Box center: ({box_center[0]:.3f}, {box_center[1]:.3f}, {box_center[2]:.3f}) nm")
    
    # Calculate center of mass for each molecule
    print(f"\nCalculating center of mass for each molecule...")
    molecules = calculate_molecule_centers(coords, residues)
    atoms_per_molecule = len(molecules[0]['atom_indices'])
    print(f"  Atoms per molecule: {atoms_per_molecule}")
    
    # Calculate geometric center of all atoms
    system_center = coords.mean(axis=0)
    print(f"  System center (geometric): ({system_center[0]:.3f}, {system_center[1]:.3f}, {system_center[2]:.3f}) nm")

    # Select molecules closest to the center
    print(f"\nSelecting {n_molecules} molecules closest to the system center...")
    selected_molecules, _, molecules_sorted = find_closest_molecules(molecules, coords, n_molecules)
    
    # Find all candidate combinations
    combinations = find_molecule_combinations(molecules_sorted, n_molecules, distance_threshold)
    
    if len(combinations) > 1:
        print(f"\n{'='*70}")
        print(f"Multiple candidate combinations found ({len(combinations)} options):")
        print(f"{'='*70}")
        print(f"  {'#':<4} {'Mol IDs':<20} {'Avg Dist (nm)':<15} {'Total Dist (nm)':<15}")
        print("  " + "-" * 55)
        for i, combo in enumerate(combinations, 1):
            mol_ids_str = ', '.join(map(str, combo['mol_ids']))
            print(f"  {i:<4} {mol_ids_str:<20} {combo['avg_distance']:<15.4f} {combo['total_distance']:<15.4f}")
        print()
        
        # Select which combination to use
        if choice is not None:
            selected_idx = choice - 1
            print(f"Using choice {choice} (specified by parameter)")
        else:
            try:
                user_input = input(f"Select combination (1-{len(combinations)}) [default: 1]: ").strip()
                selected_idx = int(user_input) - 1 if user_input else 0
            except (ValueError, EOFError):
                selected_idx = 0
                print("Using default choice: 1")
        
        if 0 <= selected_idx < len(combinations):
            selected_molecules = combinations[selected_idx]['molecules']
            print(f"\nSelected: molecules {combinations[selected_idx]['mol_ids']}")
        else:
            print(f"Invalid choice, using first combination")
            selected_molecules = combinations[0]['molecules']
    
    # Display results
    print(f"\nSelected molecules:")
    print(f"  {'Mol ID':<8} {'Center (nm)':<30} {'Distance (nm)':<12} {'Atom range':<15}")
    print("  " + "-" * 65)
    
    all_qm_atoms = []
    for mol in selected_molecules:
        center_str = f"({mol['center'][0]:.3f}, {mol['center'][1]:.3f}, {mol['center'][2]:.3f})"
        atom_min = min(mol['atom_indices'])
        atom_max = max(mol['atom_indices'])
        print(f"  {mol['id']:<8} {center_str:<30} {mol['distance_to_center']:<12.3f} "
              f"{atom_min:4d}-{atom_max:4d}")
        all_qm_atoms.extend(mol['atom_indices'])
    
    # Sort atom indices
    all_qm_atoms.sort()
    
    print(f"\nQM region:")
    print(f"  Total atoms: {len(all_qm_atoms)}")
    print(f"  Atom index range: {min(all_qm_atoms)} - {max(all_qm_atoms)}")
    if verbose:
        print(f"  Atom list: {all_qm_atoms[:10]} ... {all_qm_atoms[-10:]}")
    
    # Active region selection (QM + surrounding molecules)
    if n_active_molecules is None:
        n_active_molecules = n_molecules
        
    if n_active_molecules > n_molecules:
        print(f"\nActive region (atoms to move during optimization):")
        print(f"  QM {n_molecules} molecules + surrounding {n_active_molecules - n_molecules} molecules")
        active_molecules, _, _ = find_closest_molecules(molecules, coords, n_active_molecules)
        all_active_atoms = []
        for mol in active_molecules:
            all_active_atoms.extend(mol['atom_indices'])
        all_active_atoms.sort()
        print(f"  Total atoms: {len(all_active_atoms)}")
        print(f"  Atom index range: {min(all_active_atoms)} - {max(all_active_atoms)}")
    else:
        print(f"\nActive region:")
        print(f"  Same as QM region ({len(all_qm_atoms)} atoms)")
        all_active_atoms = all_qm_atoms
    
    # Write to files
    write_atom_list(all_qm_atoms, output)
    write_atom_list(all_active_atoms, active_output)
    
    print(f"\nOutput files:")
    print(f"  {output} - atom indices for QM region ({len(all_qm_atoms)} atoms)")
    print(f"  {active_output} - atom indices for active region in optimization ({len(all_active_atoms)} atoms)")
    
    if len(all_active_atoms) > len(all_qm_atoms):
        print(f"\nRecommendation: Including MM atoms surrounding the QM region in the active region")
        print(f"                reduces strain at the QM-MM boundary for more natural optimization")
    
    # Statistics
    if verbose:
        print(f"\nStatistics:")
        print(f"  Distance range from center for all molecules:")
        distances = [m['distance_to_center'] for m in molecules]
        print(f"    Min: {min(distances):.3f} nm")
        print(f"    Max: {max(distances):.3f} nm")
        print(f"    Mean: {np.mean(distances):.3f} nm")
        
        print(f"\n  Top 10 molecules closest to center:")
        for i, mol in enumerate(sorted(molecules, key=lambda x: x['distance_to_center'])[:10]):
            print(f"    {i+1}. Molecule {mol['id']:2d}: {mol['distance_to_center']:.3f} nm")
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)
    
    return {'qm_atoms': all_qm_atoms, 'active_atoms': all_active_atoms}


def main():
    parser = argparse.ArgumentParser(
        description='Select molecules near the center from a GRO or PDB file to define QM region'
    )
    parser.add_argument('structure_file', help='Input GRO or PDB file')
    parser.add_argument('-n', '--nmolecules', type=int, default=2,
                        help='Number of molecules to select (default: 2)')
    parser.add_argument('--active-molecules', type=int, default=None,
                        help='Number of molecules for active_atoms (QM + surrounding). '
                             'If not specified, same as QM')
    parser.add_argument('-o', '--output', default='qmatoms',
                        help='Output file name (default: qmatoms)')
    parser.add_argument('--active-output', default='active_atoms',
                        help='active_atoms output file name (default: active_atoms)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed information')
    
    args = parser.parse_args()
    
    select_qmatoms(
        structure_file=args.structure_file,
        n_molecules=args.nmolecules,
        n_active_molecules=args.active_molecules,
        output=args.output,
        active_output=args.active_output,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
