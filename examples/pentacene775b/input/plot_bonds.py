#!/usr/bin/env python
"""
Plot pentacene crystal 2D in-plane bonds
Uses 3D center-to-center distances for bond type classification
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import numpy as np
from ase.io import read
import warnings
warnings.filterwarnings('ignore')

# Load unit cell
atoms = read('pentacene.cif')
positions = atoms.get_positions()
cell = atoms.get_cell()

# Get molecular centers (even=A, odd=B)
mol_a_idx = list(range(0, 72, 2))
mol_b_idx = list(range(1, 72, 2))

center_a = positions[mol_a_idx].mean(axis=0)  # 3D center
center_b = positions[mol_b_idx].mean(axis=0)  # 3D center

# Cell vectors (3D)
a3d = cell[0]  # a-axis
b3d = cell[1]  # b-axis

# Create extended neighbor positions (3D for distance calculation)
neighbor_offsets = [
    (0, 0), (1, 0), (-1, 0), (0, 1), (0, -1),
    (1, 1), (-1, -1), (1, -1), (-1, 1),
]

all_centers_3d = []
all_centers_2d = []
all_types = []
all_cell_idx = []

for di, dj in neighbor_offsets:
    offset_3d = di * a3d + dj * b3d
    offset_2d = di * a3d[:2] + dj * b3d[:2]
    
    all_centers_3d.append(center_a + offset_3d)
    all_centers_2d.append(center_a[:2] + offset_2d)
    all_types.append('A')
    all_cell_idx.append((di, dj))
    
    all_centers_3d.append(center_b + offset_3d)
    all_centers_2d.append(center_b[:2] + offset_2d)
    all_types.append('B')
    all_cell_idx.append((di, dj))

all_centers_3d = np.array(all_centers_3d)
all_centers_2d = np.array(all_centers_2d)
all_types = np.array(all_types)

from scipy.spatial.distance import cdist
dist_matrix = cdist(all_centers_3d, all_centers_3d)  # Use 3D distances

bonds = {'type1': [], 'type2': [], 'type3': [], 'type4': []}
central_idx = [0, 1]

for i in central_idx:
    for j in range(len(all_centers_3d)):
        if i == j:
            continue
        d = dist_matrix[i, j]
        same_type = (all_types[i] == all_types[j])
        
        if 4.5 < d < 5.2 and not same_type:
            bonds['type1'].append((i, j, d))
        elif 6.0 < d < 6.5 and same_type:
            bonds['type2'].append((i, j, d))
        elif 7.0 < d < 7.5 and not same_type:
            bonds['type3'].append((i, j, d))
        elif 7.5 < d < 8.0 and same_type:
            bonds['type4'].append((i, j, d))

print('Bonds from central unit cell:')
for btype, blist in bonds.items():
    print(f'  {btype}: {len(blist)} bonds')

fig, ax = plt.subplots(figsize=(10, 10))

# Cell vectors for drawing (2D)
a = cell[0, :2]
b = cell[1, :2]

# Draw unit cell
cell_corners = np.array([[0, 0], a, a + b, b, [0, 0]])
ax.plot(cell_corners[:, 0], cell_corners[:, 1], 'k-', linewidth=2, zorder=5)
ax.fill(cell_corners[:-1, 0], cell_corners[:-1, 1], alpha=0.1, color='yellow', zorder=0)

# Draw bonds
bond_styles = {
    'type1': ('red', 'Type 1: A-B herringbone (4.89 A)'),
    'type2': ('blue', 'Type 2: A-A/B-B a-axis (6.28 A)'),
    'type3': ('green', 'Type 3: A-B diagonal (7.29 A)'),
    'type4': ('orange', 'Type 4: A-A/B-B b-axis (7.71 A)'),
}

for btype, blist in bonds.items():
    color, label = bond_styles[btype]
    for i, j, d in blist:
        c1, c2 = all_centers_2d[i], all_centers_2d[j]  # Use 2D for plotting
        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], color=color, linewidth=4, alpha=0.7, zorder=2)

# Draw molecules - use 2D positions for plotting
for i, (center, mtype) in enumerate(zip(all_centers_2d, all_types)):
    if all_cell_idx[i] == (0, 0):
        color = 'royalblue' if mtype == 'A' else 'coral'
        circle = Circle(center, 0.4, color=color, alpha=1.0, zorder=4)
        ax.add_patch(circle)
        ax.text(center[0], center[1], mtype, ha='center', va='center', fontsize=12, fontweight='bold', color='white', zorder=5)
    else:
        color = 'royalblue' if mtype == 'A' else 'coral'
        circle = Circle(center, 0.35, facecolor='white', edgecolor=color, linewidth=2, alpha=0.8, zorder=3)
        ax.add_patch(circle)
        ax.text(center[0], center[1], mtype, ha='center', va='center', fontsize=10, color=color, zorder=4)

# Cell vectors
ax.annotate('', xy=a, xytext=[0, 0], arrowprops=dict(arrowstyle='->', color='black', lw=2))
ax.text(a[0]/2, a[1]/2 - 0.8, 'a = 6.28 A', fontsize=11, ha='center')
ax.annotate('', xy=b, xytext=[0, 0], arrowprops=dict(arrowstyle='->', color='black', lw=2))
ax.text(b[0]/2 - 1.2, b[1]/2, 'b = 7.71 A', fontsize=11, ha='center', rotation=80)

# Legend
legend_elements = [
    Line2D([0], [0], color='red', linewidth=4, label='Type 1: A-B herringbone (4.89 A)'),
    Line2D([0], [0], color='blue', linewidth=4, label='Type 2: A-A/B-B a-axis (6.28 A)'),
    Line2D([0], [0], color='green', linewidth=4, label='Type 3: A-B diagonal (7.29 A)'),
    Line2D([0], [0], color='orange', linewidth=4, label='Type 4: A-A/B-B b-axis (7.71 A)'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=11, framealpha=0.9)

ax.set_xlabel('x (Angstrom)', fontsize=12)
ax.set_ylabel('y (Angstrom)', fontsize=12)
ax.set_title('Pentacene Crystal: 4 Types of 2D In-Plane Bonds\n(Unit cell with neighboring molecules)', fontsize=14)
ax.set_aspect('equal')
ax.set_xlim(-8, 14)
ax.set_ylim(-8, 16)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('pentacene_unitcell_bonds.png', dpi=200, bbox_inches='tight')
print('Saved pentacene_unitcell_bonds.png')
