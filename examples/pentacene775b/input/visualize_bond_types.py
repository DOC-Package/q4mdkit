#!/usr/bin/env python3
"""
Visualize the central pair of each bond type in pentacene crystal using VMD.
Extracts the most central pair for each of the 4 bond types and creates
a VMD visualization script.
"""

import numpy as np
from pathlib import Path


# Bond type definitions (same as select_qmatoms_pentacene.py)
BOND_TYPES = {
    1: {'name': 'herringbone', 'dist': 4.8, 'range': (4.5, 5.0), 'color': 'red'},
    2: {'name': 'second nearest', 'dist': 5.2, 'range': (5.0, 5.5), 'color': 'cyan'},
    3: {'name': 'a-axis', 'dist': 6.28, 'range': (6.0, 6.6), 'color': 'blue'},
    4: {'name': 'b-axis', 'dist': 7.71, 'range': (7.5, 8.0), 'color': 'orange'},
}

# Pentacene unit cell parameters (PENCEN01, P-1)
# a = 6.275 Å, b = 7.714 Å, c = 14.442 Å
# alpha = 76.75°, beta = 88.01°, gamma = 84.52°
def get_unit_cell_vectors():
    """Calculate unit cell vectors from crystallographic parameters."""
    a, b, c = 6.275, 7.714, 14.442
    alpha, beta, gamma = np.radians([76.75, 88.01, 84.52])
    
    # a along x
    a_vec = np.array([a, 0, 0])
    
    # b in xy plane
    b_vec = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    
    # c vector
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    c_vec = np.array([cx, cy, cz])
    
    return a_vec, b_vec, c_vec


def read_pdb_file(pdbfile):
    """Read a PDB file and return coordinates and residue numbers."""
    coords = []
    residues = []
    
    with open(pdbfile, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resnum = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
                residues.append(resnum)
    
    return np.array(coords), np.array(residues)


def get_molecule_centers(coords, residues):
    """Calculate center of each molecule."""
    unique_res = np.unique(residues)
    mol_centers = {}
    mol_atoms = {}
    
    for res in unique_res:
        mask = residues == res
        mol_centers[res] = coords[mask].mean(axis=0)
        mol_atoms[res] = np.where(mask)[0]
    
    return mol_centers, mol_atoms


def classify_bond(dist):
    """Classify bond type by distance."""
    for btype, info in BOND_TYPES.items():
        dmin, dmax = info['range']
        if dmin < dist < dmax:
            return btype
    return None


def find_central_pairs_by_type(pdb_file):
    """Find the most central pair for each bond type."""
    coords, residues = read_pdb_file(pdb_file)
    mol_centers, mol_atoms = get_molecule_centers(coords, residues)
    
    # System center
    system_center = coords.mean(axis=0)
    print(f"System center: ({system_center[0]:.2f}, {system_center[1]:.2f}, {system_center[2]:.2f}) Å")
    
    # Find all pairs and classify
    unique_res = list(mol_centers.keys())
    pairs_by_type = {1: [], 2: [], 3: [], 4: []}
    
    for i, r1 in enumerate(unique_res):
        for r2 in unique_res[i+1:]:
            dist = np.linalg.norm(mol_centers[r1] - mol_centers[r2])
            btype = classify_bond(dist)
            
            if btype is not None:
                pair_center = (mol_centers[r1] + mol_centers[r2]) / 2
                pair_dist_to_center = np.linalg.norm(pair_center - system_center)
                
                pairs_by_type[btype].append({
                    'mol1': r1,
                    'mol2': r2,
                    'center1': mol_centers[r1],
                    'center2': mol_centers[r2],
                    'atoms1': mol_atoms[r1],
                    'atoms2': mol_atoms[r2],
                    'distance': dist,
                    'pair_center': pair_center,
                    'dist_to_center': pair_dist_to_center
                })
    
    # Sort each type by distance to center and get the most central one
    central_pairs = {}
    for btype in [1, 2, 3, 4]:
        if pairs_by_type[btype]:
            pairs_by_type[btype].sort(key=lambda x: x['dist_to_center'])
            central_pairs[btype] = pairs_by_type[btype][0]
            info = BOND_TYPES[btype]
            pair = central_pairs[btype]
            print(f"Type {btype} ({info['name']}): Mol {pair['mol1']}-{pair['mol2']}, "
                  f"dist={pair['distance']:.2f} Å, from_center={pair['dist_to_center']:.2f} Å")
        else:
            print(f"Type {btype}: No pairs found")
    
    return central_pairs, coords, residues, system_center


def write_vmd_script(central_pairs, coords, residues, system_center, pdb_file, output_tcl, show_labels=True):
    """Write VMD script to visualize the central pairs only (excluding other molecules).
    
    Args:
        show_labels: If True, show text labels for bond types, axes, and legend.
    """
    
    # Collect all molecules and atoms to include
    selected_mols = set()
    mol_to_btype = {}  # Map molecule to its bond type for coloring
    for btype, pair in central_pairs.items():
        selected_mols.add(pair['mol1'])
        selected_mols.add(pair['mol2'])
        # Track which bond type each molecule belongs to (for coloring)
        if pair['mol1'] not in mol_to_btype:
            mol_to_btype[pair['mol1']] = btype
        if pair['mol2'] not in mol_to_btype:
            mol_to_btype[pair['mol2']] = btype
    
    # Read original PDB and write filtered PDB with only selected molecules
    filtered_pdb = output_tcl.replace('.tcl', '_selected.pdb')
    with open(pdb_file, 'r') as pf, open(filtered_pdb, 'w') as out:
        for line in pf:
            if line.startswith("TITLE") or line.startswith("REMARK") or line.startswith("CRYST1"):
                out.write(line)
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                resnum = int(line[22:26].strip())
                if resnum in selected_mols:
                    out.write(line)
            elif line.startswith("TER"):
                # Check if previous molecule was selected
                pass  # Skip TER records for simplicity
            elif line.startswith("END"):
                out.write(line)
        out.write("END\n")
    
    print(f"Saved {filtered_pdb} with {len(selected_mols)} molecules")
    
    with open(output_tcl, 'w') as f:
        f.write(f'''# VMD script to visualize central pairs of each bond type in pentacene
# Only selected molecules are shown (others excluded)
# Run with: vmd -e {Path(output_tcl).name}

# Load only selected molecules
mol new {filtered_pdb} type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Set background
color Display Background white
display projection Orthographic
axes location Off
display depthcue off

''')
        
        # Add each bond type pair with different colors
        for btype, pair in central_pairs.items():
            info = BOND_TYPES[btype]
            color = info['color']
            mol1, mol2 = pair['mol1'], pair['mol2']
            
            f.write(f'\n# Type {btype}: {info["name"]} ({pair["distance"]:.2f} Å)\n')
            f.write(f'# Molecules {mol1} and {mol2}\n')
            f.write(f'mol representation Licorice 0.15 12.0 12.0\n')
            f.write(f'mol color ColorID {_color_to_id(color)}\n')
            f.write(f'mol selection {{resid {mol1} {mol2}}}\n')
            f.write(f'mol material Opaque\n')
            f.write(f'mol addrep top\n')
        
        # Draw bond vectors between molecule centers
        f.write('\n# Draw bond vectors between molecule centers\n')
        for btype, pair in central_pairs.items():
            info = BOND_TYPES[btype]
            c1 = pair['center1']
            c2 = pair['center2']
            
            f.write(f'\n# Type {btype}: {info["name"]}\n')
            f.write(f'graphics top color {info["color"]}\n')
            f.write(f'graphics top cylinder {{ {c1[0]:.3f} {c1[1]:.3f} {c1[2]:.3f} }} '
                   f'{{ {c2[0]:.3f} {c2[1]:.3f} {c2[2]:.3f} }} radius 0.25 resolution 20 filled yes\n')
            
            # Label at midpoint
            if show_labels:
                mid = (c1 + c2) / 2
                f.write(f'graphics top text {{ {mid[0]:.3f} {mid[1]+1.5:.3f} {mid[2]:.3f} }} '
                       f'"T{btype}: {pair["distance"]:.2f}A" size 1.2\n')
        
        # Draw spheres at molecule centers (same radius as cylinder)
        f.write('\n# Draw molecule center markers\n')
        for btype, pair in central_pairs.items():
            info = BOND_TYPES[btype]
            f.write(f'graphics top color {info["color"]}\n')
            for center in [pair['center1'], pair['center2']]:
                f.write(f'graphics top sphere {{ {center[0]:.3f} {center[1]:.3f} {center[2]:.3f} }} '
                       f'radius 0.25 resolution 20\n')
        
        # Calculate center of selected molecules
        all_centers = []
        for pair in central_pairs.values():
            all_centers.extend([pair['center1'], pair['center2']])
        all_centers = np.array(all_centers)
        
        avg_center = all_centers.mean(axis=0)
        
        # Get unit cell vectors
        a_vec, b_vec, c_vec = get_unit_cell_vectors()
        
        # Convert molecule centers to fractional coordinates
        # Build transformation matrix (columns are a, b, c vectors)
        cell_matrix = np.column_stack([a_vec, b_vec, c_vec])
        inv_cell = np.linalg.inv(cell_matrix)
        
        # Convert centers to fractional coordinates
        frac_coords = np.array([inv_cell @ center for center in all_centers])
        
        # Find bounding box in fractional coordinates with margin
        frac_min = frac_coords.min(axis=0) - 0.3  # margin in fractional units
        frac_max = frac_coords.max(axis=0) + 0.3
        
        # Calculate box origin (corner) in Cartesian coordinates
        box_origin = frac_min[0] * a_vec + frac_min[1] * b_vec + frac_min[2] * c_vec
        
        # Calculate box dimensions in terms of unit cell vectors
        na = frac_max[0] - frac_min[0]
        nb = frac_max[1] - frac_min[1]
        nc = frac_max[2] - frac_min[2]
        
        # Draw parallelepiped box aligned with crystal axes (gray dashed lines)
        f.write('\n# Bounding box aligned with crystal axes\n')
        f.write('graphics top color gray\n')
        
        # 8 corners of the parallelepiped
        box_corners = [
            box_origin,                                           # 0
            box_origin + na * a_vec,                              # 1
            box_origin + nb * b_vec,                              # 2
            box_origin + na * a_vec + nb * b_vec,                 # 3
            box_origin + nc * c_vec,                              # 4
            box_origin + na * a_vec + nc * c_vec,                 # 5
            box_origin + nb * b_vec + nc * c_vec,                 # 6
            box_origin + na * a_vec + nb * b_vec + nc * c_vec,    # 7
        ]
        
        box_edges = [
            (0,1), (0,2), (1,3), (2,3),  # bottom face (along a and b)
            (4,5), (4,6), (5,7), (6,7),  # top face
            (0,4), (1,5), (2,6), (3,7)   # vertical edges (along c)
        ]
        
        for e1, e2 in box_edges:
            c1, c2 = box_corners[e1], box_corners[e2]
            f.write(f'graphics top line {{ {c1[0]:.3f} {c1[1]:.3f} {c1[2]:.3f} }} '
                   f'{{ {c2[0]:.3f} {c2[1]:.3f} {c2[2]:.3f} }} width 1 style dashed\n')
        
        # Draw unit cell vectors at the box origin
        axis_origin = box_origin - 0.2 * (a_vec + b_vec)
        axis_scale = 0.5
        
        f.write('\n# Unit cell vectors (a, b, c)\n')
        if show_labels:
            f.write('graphics top color black\n')
            f.write(f'graphics top text {{ {axis_origin[0]-1.5:.1f} {axis_origin[1]-1.5:.1f} {axis_origin[2]:.1f} }} "O" size 1.0\n')
        
        # a-axis (green)
        a_end = axis_origin + a_vec * axis_scale
        f.write('graphics top color green\n')
        f.write(f'graphics top cylinder {{ {axis_origin[0]:.3f} {axis_origin[1]:.3f} {axis_origin[2]:.3f} }} '
               f'{{ {a_end[0]:.3f} {a_end[1]:.3f} {a_end[2]:.3f} }} radius 0.12 resolution 20 filled yes\n')
        a_head = a_end + a_vec / np.linalg.norm(a_vec) * 0.6
        f.write(f'graphics top cone {{ {a_end[0]:.3f} {a_end[1]:.3f} {a_end[2]:.3f} }} '
               f'{{ {a_head[0]:.3f} {a_head[1]:.3f} {a_head[2]:.3f} }} radius 0.3 resolution 20\n')
        if show_labels:
            f.write(f'graphics top text {{ {a_head[0]+0.5:.3f} {a_head[1]:.3f} {a_head[2]:.3f} }} "a" size 1.2\n')
        
        # b-axis (green)
        b_end = axis_origin + b_vec * axis_scale
        f.write(f'graphics top cylinder {{ {axis_origin[0]:.3f} {axis_origin[1]:.3f} {axis_origin[2]:.3f} }} '
               f'{{ {b_end[0]:.3f} {b_end[1]:.3f} {b_end[2]:.3f} }} radius 0.12 resolution 20 filled yes\n')
        b_head = b_end + b_vec / np.linalg.norm(b_vec) * 0.6
        f.write(f'graphics top cone {{ {b_end[0]:.3f} {b_end[1]:.3f} {b_end[2]:.3f} }} '
               f'{{ {b_head[0]:.3f} {b_head[1]:.3f} {b_head[2]:.3f} }} radius 0.3 resolution 20\n')
        if show_labels:
            f.write(f'graphics top text {{ {b_head[0]:.3f} {b_head[1]+0.5:.3f} {b_head[2]:.3f} }} "b" size 1.2\n')
        
        # c-axis (blue)
        c_end = axis_origin + c_vec * axis_scale
        f.write('graphics top color blue\n')
        f.write(f'graphics top cylinder {{ {axis_origin[0]:.3f} {axis_origin[1]:.3f} {axis_origin[2]:.3f} }} '
               f'{{ {c_end[0]:.3f} {c_end[1]:.3f} {c_end[2]:.3f} }} radius 0.12 resolution 20 filled yes\n')
        c_head = c_end + c_vec / np.linalg.norm(c_vec) * 0.6
        f.write(f'graphics top cone {{ {c_end[0]:.3f} {c_end[1]:.3f} {c_end[2]:.3f} }} '
               f'{{ {c_head[0]:.3f} {c_head[1]:.3f} {c_head[2]:.3f} }} radius 0.3 resolution 20\n')
        if show_labels:
            f.write(f'graphics top text {{ {c_head[0]:.3f} {c_head[1]:.3f} {c_head[2]+0.5:.3f} }} "c" size 1.2\n')
        
        # Add legend at a fixed position relative to molecules
        if show_labels:
            f.write(f'''
# Legend
graphics top color black
graphics top text {{ {avg_center[0]+10:.1f} {avg_center[1]-15:.1f} {avg_center[2]:.1f} }} "Bond Types:" size 1.2
''')
            for i, (btype, info) in enumerate(BOND_TYPES.items()):
                y_offset = avg_center[1] - 15 - 2.5 * (i + 1)
                pair = central_pairs.get(btype)
                dist_str = f"{pair['distance']:.2f}" if pair else "N/A"
                f.write(f'graphics top color {info["color"]}\n')
                f.write(f'graphics top text {{ {avg_center[0]+10:.1f} {y_offset:.1f} {avg_center[2]:.1f} }} '
                       f'"Type {btype}: {info["name"]} ({dist_str} A)" size 1.0\n')
        
        f.write(f'''
# Set view centered on selected molecules
display resetview
rotate x by -60
rotate y by 10
scale by 1.2

puts "Visualization loaded!"
puts "Only selected molecules are shown."
puts "Bond types:"
''')
        for btype, pair in central_pairs.items():
            info = BOND_TYPES[btype]
            f.write(f'puts "  Type {btype} ({info["name"]}): Mol {pair["mol1"]}-{pair["mol2"]}, {pair["distance"]:.2f} A"\n')
    
    print(f"Saved {output_tcl}")


def _color_to_id(color_name):
    """Convert color name to VMD color ID."""
    color_map = {
        'red': 1,
        'orange': 3,
        'yellow': 4,
        'green': 7,
        'cyan': 10,
        'blue': 0,
        'purple': 11,
        'pink': 9,
        'white': 8,
        'black': 16,
        'gray': 2,
    }
    return color_map.get(color_name, 2)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Visualize central pairs of each bond type in pentacene crystal'
    )
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('-o', '--output', default='vmd_bond_types.tcl',
                        help='Output VMD script (default: vmd_bond_types.tcl)')
    parser.add_argument('--no-labels', action='store_true',
                        help='Disable text labels (bond types, axes, legend)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Visualizing central pairs of each bond type")
    print("=" * 70)
    
    central_pairs, coords, residues, system_center = find_central_pairs_by_type(args.pdb_file)
    
    write_vmd_script(central_pairs, coords, residues, system_center, 
                     args.pdb_file, args.output, show_labels=not args.no_labels)
    
    print("\n" + "=" * 70)
    print(f"To visualize: vmd -e {args.output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
