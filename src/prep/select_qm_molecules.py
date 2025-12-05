#!/usr/bin/env python3
"""
Select two adjacent molecules closest to the center of the system from a GRO file.
Output atom indices to be used for the QM region.
"""

import numpy as np
import argparse


def read_gro_file(grofile):
    """
    Read a GRO file and return coordinates, residue numbers, and box information.
    
    Returns:
        coords: numpy array (N, 3) - atomic coordinates (nm)
        residues: numpy array (N,) - residue numbers
        box: numpy array (3,) - box dimensions (nm)
        natoms: int - total number of atoms
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


def find_closest_molecules(molecules, box, n_molecules=2):
    """
    Find n molecules closest to the center of the box.
    
    Args:
        molecules: list of molecules
        box: box dimensions
        n_molecules: number of molecules to select
    
    Returns:
        selected_molecules: list of selected molecules (sorted by distance)
    """
    box_center = box / 2.0
    
    # Calculate distance from each molecule center to the box center
    for mol in molecules:
        distance = np.linalg.norm(mol['center'] - box_center)
        mol['distance_to_center'] = distance
    
    # Sort by distance
    molecules_sorted = sorted(molecules, key=lambda x: x['distance_to_center'])
    
    # Select the n closest molecules
    selected = molecules_sorted[:n_molecules]
    
    return selected


def write_atom_list(atom_indices, filename):
    """Write atom indices to a file."""
    with open(filename, 'w') as f:
        f.write(' '.join(map(str, atom_indices)) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Select molecules near the center from a GRO file to define QM region'
    )
    parser.add_argument('grofile', help='Input GRO file')
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
    
    print("=" * 70)
    print("QM Region Selection: Selecting molecules near the center of the system")
    print("=" * 70)
    
    # Read GRO file
    print(f"\nReading: {args.grofile}")
    coords, residues, box, natoms, title = read_gro_file(args.grofile)
    
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
    
    # Select molecules closest to the center
    print(f"\nSelecting {args.nmolecules} molecules closest to the center...")
    selected_molecules = find_closest_molecules(molecules, box, args.nmolecules)
    
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
    if args.verbose:
        print(f"  Atom list: {all_qm_atoms[:10]} ... {all_qm_atoms[-10:]}")
    
    # Active region selection (QM + surrounding molecules)
    n_active_molecules = args.active_molecules if args.active_molecules is not None else args.nmolecules
    if n_active_molecules > args.nmolecules:
        print(f"\nActive region (atoms to move during optimization):")
        print(f"  QM {args.nmolecules} molecules + surrounding {n_active_molecules - args.nmolecules} molecules")
        active_molecules = find_closest_molecules(molecules, box, n_active_molecules)
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
    write_atom_list(all_qm_atoms, args.output)
    write_atom_list(all_active_atoms, args.active_output)
    
    print(f"\nOutput files:")
    print(f"  {args.output} - atom indices for QM region ({len(all_qm_atoms)} atoms)")
    print(f"  {args.active_output} - atom indices for active region in optimization ({len(all_active_atoms)} atoms)")
    
    if len(all_active_atoms) > len(all_qm_atoms):
        print(f"\nRecommendation: Including MM atoms surrounding the QM region in the active region")
        print(f"                reduces strain at the QM-MM boundary for more natural optimization")
    
    # Statistics
    if args.verbose:
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


if __name__ == "__main__":
    main()
