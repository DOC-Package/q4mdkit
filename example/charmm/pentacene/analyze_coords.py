#!/usr/bin/env python3
"""Calculate initial energy with proper PBC setup and check distances"""
import numpy as np

# Read CRD file
def read_crd(filename):
    atoms = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    # Skip header and count line
    idx = 0
    while idx < len(lines) and not lines[idx].strip().isdigit():
        idx += 1
    
    natoms = int(lines[idx].strip())
    idx += 1
    
    # Read coordinates (extended format: I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A4,F20.10)
    for i in range(natoms):
        line = lines[idx + i]
        if len(line) < 50:
            continue
        atom_id = int(line[0:10])
        resid = int(line[12:20])
        resname = line[22:30].strip()
        x = float(line[30:50])
        y = float(line[50:70])
        z = float(line[70:90])
        segid = line[92:100].strip() if len(line) > 92 else ""
        resid2 = line[102:106].strip() if len(line) > 102 else ""
        atoms.append({'id': atom_id, 'resid': resid, 'resname': resname,
                      'x': x, 'y': y, 'z': z, 'segid': segid})
    
    # Read box from last line
    last_line = lines[-1].strip().split()
    if len(last_line) == 6:
        box = [float(x) for x in last_line]
    else:
        box = None
    
    return atoms, box

print("="*70)
print("Analyzing pentacene_nowrap.crd")
print("="*70)
atoms_nowrap, box_nowrap = read_crd("pentacene_nowrap.crd")
print(f"Number of atoms: {len(atoms_nowrap)}")
if box_nowrap:
    print(f"Box: {box_nowrap[0]:.3f} {box_nowrap[1]:.3f} {box_nowrap[2]:.3f} Å")
    print(f"Angles: {box_nowrap[3]:.2f}° {box_nowrap[4]:.2f}° {box_nowrap[5]:.2f}°")

# Calculate minimum distances
coords = np.array([[a['x'], a['y'], a['z']] for a in atoms_nowrap])
print(f"\nCoordinate ranges:")
print(f"  X: {coords[:, 0].min():.3f} to {coords[:, 0].max():.3f} Å")
print(f"  Y: {coords[:, 1].min():.3f} to {coords[:, 1].max():.3f} Å")
print(f"  Z: {coords[:, 2].min():.3f} to {coords[:, 2].max():.3f} Å")

# Find minimum distance between atoms (excluding bonded neighbors)
min_dist = float('inf')
min_pair = None
for i in range(len(coords)):
    for j in range(i+5, len(coords)):  # Skip nearby atoms (likely bonded)
        dist = np.linalg.norm(coords[i] - coords[j])
        if dist < min_dist:
            min_dist = dist
            min_pair = (i, j)

print(f"\nMinimum distance between non-bonded atoms: {min_dist:.3f} Å")
if min_pair:
    i, j = min_pair
    print(f"  Atom {atoms_nowrap[i]['id']} ({atoms_nowrap[i]['resname']}) <-> Atom {atoms_nowrap[j]['id']} ({atoms_nowrap[j]['resname']})")

# Check for very short distances
short_pairs = []
for i in range(len(coords)):
    for j in range(i+5, len(coords)):
        dist = np.linalg.norm(coords[i] - coords[j])
        if dist < 2.0:
            short_pairs.append((i, j, dist))

print(f"\nNumber of atom pairs with distance < 2.0 Å: {len(short_pairs)}")
if short_pairs:
    print("First 10 short distances:")
    for i, j, d in sorted(short_pairs, key=lambda x: x[2])[:10]:
        print(f"  {atoms_nowrap[i]['id']:4d} - {atoms_nowrap[j]['id']:4d}: {d:.3f} Å")

print("\n" + "="*70)
print("Analyzing pentacene.crd (wrapped)")
print("="*70)
atoms_wrap, box_wrap = read_crd("pentacene.crd")
print(f"Number of atoms: {len(atoms_wrap)}")
if box_wrap:
    print(f"Box: {box_wrap[0]:.3f} {box_wrap[1]:.3f} {box_wrap[2]:.3f} Å")
    print(f"Angles: {box_wrap[3]:.2f}° {box_wrap[4]:.2f}° {box_wrap[5]:.2f}°")

coords_wrap = np.array([[a['x'], a['y'], a['z']] for a in atoms_wrap])
print(f"\nCoordinate ranges:")
print(f"  X: {coords_wrap[:, 0].min():.3f} to {coords_wrap[:, 0].max():.3f} Å")
print(f"  Y: {coords_wrap[:, 1].min():.3f} to {coords_wrap[:, 1].max():.3f} Å")
print(f"  Z: {coords_wrap[:, 2].min():.3f} to {coords_wrap[:, 2].max():.3f} Å")

min_dist = float('inf')
min_pair = None
for i in range(len(coords_wrap)):
    for j in range(i+5, len(coords_wrap)):
        dist = np.linalg.norm(coords_wrap[i] - coords_wrap[j])
        if dist < min_dist:
            min_dist = dist
            min_pair = (i, j)

print(f"\nMinimum distance between non-bonded atoms: {min_dist:.3f} Å")
if min_pair:
    i, j = min_pair
    print(f"  Atom {atoms_wrap[i]['id']} ({atoms_wrap[i]['resname']}) <-> Atom {atoms_wrap[j]['id']} ({atoms_wrap[j]['resname']})")

short_pairs_wrap = []
for i in range(len(coords_wrap)):
    for j in range(i+5, len(coords_wrap)):
        dist = np.linalg.norm(coords_wrap[i] - coords_wrap[j])
        if dist < 2.0:
            short_pairs_wrap.append((i, j, dist))

print(f"\nNumber of atom pairs with distance < 2.0 Å: {len(short_pairs_wrap)}")
if short_pairs_wrap:
    print("First 10 short distances:")
    for i, j, d in sorted(short_pairs_wrap, key=lambda x: x[2])[:10]:
        print(f"  {atoms_wrap[i]['id']:4d} - {atoms_wrap[j]['id']:4d}: {d:.3f} Å")

print("\n" + "="*70)
print("CONCLUSION:")
print("="*70)
if len(short_pairs) < len(short_pairs_wrap):
    print("✓ pentacene_nowrap.crd has FEWER overlapping atoms")
    print(f"  Use: COOR = 'pentacene_nowrap.crd' in run_md.py")
else:
    print("✗ Both files have atom overlaps - need to regenerate structure")
    print(f"  Suggestion: Run cif2crd with different options")
