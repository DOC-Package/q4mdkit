#!/usr/bin/env python
"""
Analyze 2D in-plane bonding patterns in pentacene crystal
Pentacene has 2 molecules per unit cell (site A and site B)
This script analyzes all 4 unique dimer types in the ab-plane
"""
import numpy as np
from pathlib import Path
from collections import defaultdict

def read_pdb(pdb_file):
    """Read PDB file and extract coordinates by residue"""
    atoms = []
    box = None
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_data = {
                    'resid': int(line[22:26].strip()),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'element': line[76:78].strip()
                }
                atoms.append(atom_data)
            elif line.startswith('CRYST1'):
                box = [
                    float(line[6:15].strip()),
                    float(line[15:24].strip()),
                    float(line[24:33].strip())
                ]
    
    residues = {}
    for atom in atoms:
        resid = atom['resid']
        if resid not in residues:
            residues[resid] = []
        residues[resid].append(atom)
    
    return residues, box

def get_molecular_orientation(atoms):
    """Calculate the principal axis (long axis) of pentacene molecule"""
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    centered = coords - np.mean(coords, axis=0)
    
    # SVD to find principal axes
    U, S, Vt = np.linalg.svd(centered)
    # First principal axis (longest)
    return Vt[0]

def analyze_crystal_bonds_complete(pdb_file, distance_cutoff=0.55):
    """
    Analyze all bonding patterns in pentacene crystal ab-plane
    """
    residues, box = read_pdb(pdb_file)
    box = np.array(box)
    
    print(f"Box dimensions: {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} Å")
    print(f"Box dimensions: {box[0]/10:.3f} x {box[1]/10:.3f} x {box[2]/10:.3f} nm")
    
    n_molecules = len(residues)
    print(f"Total molecules in supercell: {n_molecules}")
    
    # Calculate COM and orientation for each molecule
    mol_indices = sorted(residues.keys())
    mol_data = []
    
    for resid in mol_indices:
        coords = np.array([[a['x'], a['y'], a['z']] for a in residues[resid]])
        com = np.mean(coords, axis=0)
        orientation = get_molecular_orientation(residues[resid])
        mol_data.append({
            'resid': resid,
            'com': com,
            'orientation': orientation,
            'coords': coords
        })
    
    # Identify the two sites in the unit cell by their z-coordinate pattern
    # In pentacene, molecules alternate between two orientations
    box_center = box / 2
    distances_to_center = [np.linalg.norm(m['com'] - box_center) for m in mol_data]
    central_idx = np.argmin(distances_to_center)
    central_mol = mol_data[central_idx]
    
    print(f"\nCentral molecule: index={central_idx}, resid={central_mol['resid']}")
    print(f"Central COM: {central_mol['com']}")
    
    # Classify molecules into site A and site B based on orientation
    # The two sites have different herringbone angles
    central_orient = central_mol['orientation']
    
    # Find all molecules in the same ab-plane (similar z)
    z_tolerance = 2.0  # Angstrom
    same_layer = []
    cutoff_angstrom = distance_cutoff * 10
    
    for i, mol in enumerate(mol_data):
        if i == central_idx:
            continue
        
        delta = mol['com'] - central_mol['com']
        
        # Apply PBC
        for j in range(3):
            if delta[j] > box[j] / 2:
                delta[j] -= box[j]
            elif delta[j] < -box[j] / 2:
                delta[j] += box[j]
        
        distance = np.linalg.norm(delta)
        
        # Check if in same ab-plane and within cutoff
        if abs(delta[2]) < z_tolerance and distance < cutoff_angstrom:
            # Determine if same site or different site by orientation
            dot_product = abs(np.dot(mol['orientation'], central_orient))
            same_site = dot_product > 0.9  # Similar orientation = same site
            
            same_layer.append({
                'index': i,
                'resid': mol['resid'],
                'distance': distance,
                'delta': delta,
                'same_site': same_site,
                'orientation': mol['orientation'],
                'dot_product': dot_product
            })
    
    # Sort by distance
    same_layer.sort(key=lambda x: x['distance'])
    
    print(f"\n{'='*80}")
    print("IN-PLANE NEIGHBORS (ab-plane)")
    print(f"{'='*80}")
    print(f"\n{'Idx':>4} | {'ResID':>5} | {'Dist(nm)':>8} | {'Δx(Å)':>7} | {'Δy(Å)':>7} | {'Δz(Å)':>7} | Site")
    print("-" * 80)
    
    for n in same_layer:
        site_str = "SAME" if n['same_site'] else "DIFF"
        print(f"{n['index']:4d} | {n['resid']:5d} | {n['distance']/10:8.4f} | {n['delta'][0]:7.2f} | {n['delta'][1]:7.2f} | {n['delta'][2]:7.2f} | {site_str}")
    
    # Group into unique bond types
    print(f"\n{'='*80}")
    print("UNIQUE BOND TYPES IN AB-PLANE")
    print(f"{'='*80}")
    
    bond_types = []
    tolerance_dist = 0.02  # nm
    tolerance_angle = 10  # degrees
    
    for n in same_layer:
        delta = n['delta']
        angle = np.arctan2(delta[1], delta[0]) * 180 / np.pi
        dist = n['distance'] / 10  # nm
        
        # For classification: consider both distance, angle, and site type
        norm_angle = angle % 180
        
        found = False
        for bt in bond_types:
            angle_diff = min(abs(norm_angle - bt['norm_angle']), 180 - abs(norm_angle - bt['norm_angle']))
            same_type = (n['same_site'] == bt['same_site'])
            
            if abs(dist - bt['distance']) < tolerance_dist and angle_diff < tolerance_angle and same_type:
                bt['members'].append(n)
                found = True
                break
        
        if not found:
            bond_types.append({
                'distance': dist,
                'norm_angle': norm_angle,
                'same_site': n['same_site'],
                'members': [n]
            })
    
    # Sort bond types by distance
    bond_types.sort(key=lambda x: x['distance'])
    
    print(f"\nFound {len(bond_types)} unique bond types:")
    
    for i, bt in enumerate(bond_types, 1):
        site_type = "Same-site (A-A or B-B)" if bt['same_site'] else "Cross-site (A-B)"
        print(f"\n{'='*80}")
        print(f"BOND TYPE {i}: {site_type}")
        print(f"{'='*80}")
        print(f"  Distance: {bt['distance']:.4f} nm ({bt['distance']*10:.3f} Å)")
        print(f"  Angle: {bt['norm_angle']:.1f}° (normalized to 0-180°)")
        print(f"  Multiplicity: {len(bt['members'])}")
        
        for m in bt['members']:
            angle = np.arctan2(m['delta'][1], m['delta'][0]) * 180 / np.pi
            print(f"    Mol {m['index']:3d} (res {m['resid']:3d}): "
                  f"Δ=({m['delta'][0]:6.2f}, {m['delta'][1]:6.2f}, {m['delta'][2]:6.2f}) Å, "
                  f"angle={angle:.1f}°")
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY: 2D IN-PLANE BONDING IN PENTACENE")
    print(f"{'='*80}")
    
    same_site_types = [bt for bt in bond_types if bt['same_site']]
    cross_site_types = [bt for bt in bond_types if not bt['same_site']]
    
    print(f"\nSame-site bonds (A-A or B-B): {len(same_site_types)} types")
    for i, bt in enumerate(same_site_types, 1):
        print(f"  Type S{i}: d={bt['distance']:.4f} nm, angle={bt['norm_angle']:.1f}°, count={len(bt['members'])}")
    
    print(f"\nCross-site bonds (A-B): {len(cross_site_types)} types")
    for i, bt in enumerate(cross_site_types, 1):
        print(f"  Type C{i}: d={bt['distance']:.4f} nm, angle={bt['norm_angle']:.1f}°, count={len(bt['members'])}")
    
    print(f"\nTotal unique bond types: {len(bond_types)}")
    
    return {
        'central_mol': central_mol,
        'neighbors': same_layer,
        'bond_types': bond_types,
        'same_site_types': same_site_types,
        'cross_site_types': cross_site_types
    }


if __name__ == "__main__":
    input_dir = Path(__file__).parent
    pdb_file = str(input_dir / "pentacene.pdb")
    
    print("="*80)
    print("PENTACENE CRYSTAL 2D IN-PLANE BONDING ANALYSIS")
    print("(Analyzing all bond types in the ab-plane)")
    print("="*80)
    
    result = analyze_crystal_bonds_complete(pdb_file, distance_cutoff=0.70)
