# mol_utils.py - Molecular Graph and Coordinate Utilities

## Overview

Core utility functions for molecular crystal structure manipulation, including graph-based molecule detection, coordinate transformations, and graph isomorphism matching. These functions support the workflow in `cif2gro.py` and `cif2crd.py`.

## Functions

### 1. `build_graph(atoms, scale=1.10)`

Build an undirected covalent bond graph from an ASE Atoms object.

**Parameters:**
- `atoms` (Atoms): ASE Atoms object
- `scale` (float): Scale factor for covalent radii (default: 1.10)

**Returns:**
- `nx.Graph`: NetworkX graph with node attributes:
  - `Z`: Atomic number
  - `deg`: Degree (number of bonds)
  - `hnb`: Number of hydrogen neighbors

**Algorithm:**
1. Calculate cutoff distances: `cutoff[i] = covalent_radii[Z[i]] * scale`
2. Use ASE NeighborList to find bonds under periodic boundaries
3. Two atoms i and j are bonded if: `distance(i,j) <= cutoff[i] + cutoff[j]`
4. Annotate nodes with invariants for isomorphism matching

**Example:**
```python
from mol_utils import build_graph
from ase.io import read

atoms = read("pentacene.cif")
G = build_graph(atoms, scale=1.10)

print(f"Atoms: {len(G.nodes)}")
print(f"Bonds: {len(G.edges)}")
print(f"Node 0: Z={G.nodes[0]['Z']}, degree={G.nodes[0]['deg']}, H-neighbors={G.nodes[0]['hnb']}")
```

---

### 2. `split_molecules_pbc(atoms, scale=1.10)`

Identify molecules as connected components under periodic boundary conditions.

**Parameters:**
- `atoms` (Atoms): ASE Atoms object with PBC
- `scale` (float): Scale factor for covalent radii (default: 1.10)

**Returns:**
- `List[List[int]]`: List of lists, each containing atom indices for one molecule

**Algorithm:**
1. Build covalent bond graph with PBC
2. Find connected components using NetworkX
3. Each component is one molecule
4. Sort atom indices within each component

**Example:**
```python
from mol_utils import split_molecules_pbc
from ase.io import read

atoms = read("pentacene_crystal.cif")
molecules = split_molecules_pbc(atoms, scale=1.10)

print(f"Detected {len(molecules)} molecules")
print(f"First molecule has {len(molecules[0])} atoms")
print(f"Atom indices: {molecules[0]}")
```

**Use Cases:**
- Molecular crystal detection
- Quality check (should detect separate molecules, not a single network)
- Preparing for per-molecule processing

---

### 3. `unwrap_to_single_image(atoms, idxs)`

Unwrap a molecule so all its atoms are in the same periodic image (no bonds crossing boundaries).

**Parameters:**
- `atoms` (Atoms): ASE Atoms object with PBC
- `idxs` (List[int]): Indices of atoms in the molecule

**Returns:**
- `np.ndarray`: Nx3 array of unwrapped Cartesian positions (Å)

**Algorithm:**
1. Get fractional coordinates of molecule atoms
2. Choose first atom as reference
3. For each other atom:
   - Calculate fractional distance to reference
   - If distance > 0.5 in any direction, shift by -1 or +1 cells
   - Result: all atoms within ±0.5 of reference
4. Convert to Cartesian coordinates

**Key Insight:**
```python
# Before unwrap (molecule split across boundary):
fractional = [[0.05, 0.1, 0.2],   # Atom 0
              [0.95, 0.1, 0.2]]   # Atom 1 (looks far, but bonded via PBC)

# After unwrap (continuous):
fractional = [[0.05, 0.1, 0.2],   # Atom 0 (reference)
              [-0.05, 0.1, 0.2]]  # Atom 1 (shifted to be near reference)
```

**Example:**
```python
from mol_utils import split_molecules_pbc, unwrap_to_single_image
from ase.io import read

atoms = read("crystal.cif")
molecules = split_molecules_pbc(atoms)

# Unwrap first molecule
unwrapped_pos = unwrap_to_single_image(atoms, molecules[0])
print(f"Unwrapped positions shape: {unwrapped_pos.shape}")
# Output: (36, 3) for 36 atoms
```

**Use Cases:**
- Preparing coordinates for PDB/GRO output (bonds should not cross box)
- Visualization (molecules appear as continuous units)
- RMSD calculations (avoid artificial distances)

---

### 4. `recenter_system(unwrapped_positions, cell)`

Recenter all unwrapped molecules to form a compact cluster around the box center.

**Parameters:**
- `unwrapped_positions` (List[np.ndarray]): List of unwrapped positions for each molecule (Å)
- `cell` (np.ndarray): 3x3 cell matrix (Å)

**Returns:**
- `List[np.ndarray]`: Recentered positions for each molecule

**Algorithm:**
1. Calculate center of mass (COM) of entire system
2. Calculate box center: `box_center = 0.5 * (a + b + c)`
3. Shift all molecules: `shift = box_center - system_COM`
4. Apply shift to all molecule positions

**Example:**
```python
from mol_utils import recenter_system
import numpy as np

# Unwrapped positions for 3 molecules
unwrapped = [np.array([[1, 2, 3], [1.5, 2.5, 3.5]]),  # Mol 1
             np.array([[10, 20, 30], [10.5, 20.5, 30.5]]),  # Mol 2
             np.array([[-5, -10, -15], [-5.5, -10.5, -15.5]])]  # Mol 3

cell = np.diag([50, 50, 50])  # 50x50x50 Å box

recentered = recenter_system(unwrapped, cell)
# Now system COM is at [25, 25, 25] (box center)
```

**Use Cases:**
- Visualization (compact cluster easier to view)
- Preparing for analysis tools that expect centered systems
- Reducing coordinate range for numerical stability

---

### 5. `wrap_to_box(unwrapped_positions, cell)`

Wrap molecule positions back into the primary simulation box [0, 1) while maintaining molecular integrity.

**Parameters:**
- `unwrapped_positions` (List[np.ndarray]): List of unwrapped positions (Å)
- `cell` (np.ndarray): 3x3 cell matrix (Å)

**Returns:**
- `List[np.ndarray]`: Wrapped positions for each molecule

**Algorithm:**
1. Convert each molecule to fractional coordinates
2. Calculate COM in fractional space
3. Find integer shift to bring COM into [0, 1)
4. Apply shift to all atoms in molecule (keeps molecule together)
5. Convert back to Cartesian

**Example:**
```python
from mol_utils import wrap_to_box
import numpy as np

# Molecule with COM at fractional [-0.3, 0.5, 1.2]
unwrapped = [np.array([[-3, 5, 12], [-2.5, 5.5, 12.5]])]
cell = np.diag([10, 10, 10])

wrapped = wrap_to_box(unwrapped, cell)
# COM now at [0.7, 0.5, 0.2] (wrapped into [0,1))
```

**Important:**
- Each molecule moves as a rigid unit (no internal distortion)
- Bonds within molecules remain intact
- Different from per-atom wrapping

---

### 6. `deterministic_template_order(G, atoms, idxs)`

Choose a stable internal atom order for a template molecule using BFS from a unique root.

**Parameters:**
- `G` (nx.Graph): Molecular graph
- `atoms` (Atoms): ASE Atoms object
- `idxs` (List[int]): Global indices of atoms in this molecule

**Returns:**
- `List[int]`: List of local indices (0..n-1) in the chosen order

**Root Selection Criteria (priority order):**
1. Atomic number (Z) - descending (heavy atoms first)
2. Degree (number of bonds) - descending (well-connected first)
3. H-neighbor count - ascending (fewer H neighbors first)
4. Index - ascending (tie-breaker)

**Algorithm:**
1. Rank all atoms by selection criteria
2. Choose highest-ranked atom as root
3. Perform BFS from root to get traversal order
4. Convert global indices to local (0..n-1)

**Example:**
```python
from mol_utils import build_graph, deterministic_template_order
from ase.io import read

atoms = read("pentacene.cif")
mol_indices = list(range(36))  # First molecule
G = build_graph(atoms)

order = deterministic_template_order(G, atoms, mol_indices)
print(f"Atom order: {order}")
# [0, 5, 10, 15, 1, 6, ...] (BFS from root)
```

**Why Deterministic?**
- Same molecule always produces same order
- Different runs produce identical results
- Essential for consistent topology application

---

### 7. `best_isomorphism(src, ref, src_idxs, ref_idxs, scale=1.10)`

Find graph isomorphism mapping from source to reference molecule, minimizing RMSD among all valid mappings.

**Parameters:**
- `src` (Atoms): Source molecule
- `ref` (Atoms): Reference (template) molecule
- `src_idxs` (List[int]): Local indices for source (0..n-1)
- `ref_idxs` (List[int]): Local indices for reference (0..n-1)
- `scale` (float): Scale factor for covalent radii (default: 1.10)

**Returns:**
- `Dict[int, int]`: Mapping from src local indices to ref local indices

**Algorithm:**
1. Build graphs for both molecules
2. Find all graph isomorphisms matching by:
   - Atomic number (Z)
   - Degree (number of bonds)
   - H-neighbor count
3. For each isomorphism:
   - Apply Kabsch algorithm (optimal rotation)
   - Calculate RMSD
4. Return mapping with minimum RMSD

**Node Matching:**
```python
def node_match(n1, n2):
    return (n1["Z"] == n2["Z"]) and \
           (n1["deg"] == n2["deg"]) and \
           (n1["hnb"] == n2["hnb"])
```

**Kabsch Algorithm:**
1. Center both structures at origin
2. Calculate cross-covariance matrix: `H = src^T @ ref`
3. SVD decomposition: `H = U @ S @ V^T`
4. Optimal rotation: `R = U @ V^T`
5. Calculate RMSD after rotation

**Example:**
```python
from mol_utils import best_isomorphism
from ase import Atoms
import numpy as np

# Two pentacene molecules
mol1 = Atoms('C22H14', positions=pos1)
mol2 = Atoms('C22H14', positions=pos2)

mapping = best_isomorphism(
    mol1, mol2,
    src_idxs=list(range(36)),
    ref_idxs=list(range(36)),
    scale=1.10
)

# mapping: {0: 5, 1: 10, 2: 15, ...}
# Atom 0 in mol1 corresponds to atom 5 in mol2
```

**Handling Symmetry:**
- Symmetric molecules have multiple valid isomorphisms
- RMSD-based selection ensures deterministic choice
- Physically equivalent but computationally consistent

**Error Handling:**
```python
if best_map is None:
    raise RuntimeError(
        "Graph isomorphism failed; "
        "try adjusting --bond-scale or ensure H atoms are present."
    )
```

---

## Typical Workflow

### Complete Pipeline Example

```python
from mol_utils import (
    split_molecules_pbc,
    unwrap_to_single_image,
    recenter_system,
    wrap_to_box,
    build_graph,
    deterministic_template_order,
    best_isomorphism
)
from ase.io import read
import numpy as np

# 1. Load crystal structure
atoms = read("pentacene_crystal.cif")

# 2. Detect molecules
molecules = split_molecules_pbc(atoms, scale=1.10)
print(f"Detected {len(molecules)} molecules")

# 3. Unwrap each molecule
unwrapped = []
for mol_idxs in molecules:
    pos = unwrap_to_single_image(atoms, mol_idxs)
    unwrapped.append(pos)

# 4. Recenter system
unwrapped = recenter_system(unwrapped, np.array(atoms.get_cell()))

# 5. Wrap back into box
unwrapped = wrap_to_box(unwrapped, np.array(atoms.get_cell()))

# 6. Define template order (first molecule)
G = build_graph(atoms, scale=1.10)
template_order = deterministic_template_order(G, atoms, molecules[0])

# 7. Map other molecules to template
mappings = {}
for i, mol_idxs in enumerate(molecules):
    if i == 0:
        mappings[i] = template_order
    else:
        mapping = best_isomorphism(
            atoms[mol_idxs], atoms[molecules[0]],
            list(range(len(mol_idxs))),
            list(range(len(molecules[0]))),
            scale=1.10
        )
        # Convert mapping to order
        inv = {v: k for k, v in mapping.items()}
        order = [inv[j] for j in template_order]
        mappings[i] = order

# Now all molecules have consistent atom ordering
```

## Performance Considerations

### NeighborList Efficiency
- Use `skin=0.0` for static structures
- PBC calculations are optimized in ASE

### Graph Isomorphism
- Complexity: exponential in worst case
- Typically fast for small molecules (<100 atoms)
- Symmetric molecules may have many isomorphisms
- RMSD calculation is O(n) per isomorphism

### Memory Usage
- Unwrapped positions stored per molecule
- Graph isomorphism requires temporary coordinate arrays
- Typical usage: <100 MB for 1000-molecule systems

## Common Issues

### "Single connected component" Warning
- Molecules appear as one large network
- **Solution**: Decrease `--bond-scale` (e.g., 1.05)
- **Cause**: Intermolecular distances too small

### "Graph isomorphism failed" Error
- Cannot find valid atom mapping
- **Solutions**:
  1. Adjust `--bond-scale`
  2. Ensure all H atoms present
  3. Check molecular integrity
  4. Verify consistent molecule types

### Molecules Split Incorrectly
- Some atoms missing from detected molecules
- **Solution**: Increase `--bond-scale` (e.g., 1.15)
- **Cause**: Bonds not detected (loose structure)

### High RMSD in Isomorphism
- Warning if RMSD > 0.5 Å after matching
- **Possible causes**:
  1. Different molecular conformations
  2. Numerical precision issues
  3. Wrong molecule type comparison

## Dependencies

- **NumPy**: Array operations, linear algebra
- **ASE**: Atomic structure representation, NeighborList
- **NetworkX**: Graph algorithms, isomorphism matching

## Related Files

- **cif2gro.py**: Uses these functions for GRO file generation
- **cif2crd.py**: Uses these functions for CHARMM CRD format
- **CIF2GRO_USAGE.md**: Higher-level documentation for end users
