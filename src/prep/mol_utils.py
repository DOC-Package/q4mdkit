import numpy as np
from typing import List, Dict
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.data import covalent_radii
import networkx as nx

def build_graph(atoms: Atoms, scale: float = 1.10) -> nx.Graph:
    """
    Build covalent bond graph from ASE Atoms.
    See MOL_UTILS_USAGE.md for details.
    """
    Z = atoms.get_atomic_numbers()
    cutoffs = [covalent_radii[z] * scale for z in Z]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    G = nx.Graph()
    for i, z in enumerate(Z):
        G.add_node(i, Z=int(z))
    for i in range(len(Z)):
        idxs, _ = nl.get_neighbors(i)
        for j in idxs:
            if i < j:
                G.add_edge(i, j)
    
    # Annotate with simple invariants for matching
    deg = dict(G.degree())
    for i in G.nodes:
        G.nodes[i]["deg"] = int(deg[i])
        hnb = sum(1 for j in G.neighbors(i) if atoms.numbers[j] == 1)
        G.nodes[i]["hnb"] = int(hnb)
    
    return G


def split_molecules_pbc(atoms: Atoms, scale: float = 1.10) -> List[List[int]]:
    """
    Identify molecules as connected components under PBC.
    See MOL_UTILS_USAGE.md for details.
    """
    Z = atoms.get_atomic_numbers()
    cutoffs = [covalent_radii[z] * scale for z in Z]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
    nl.update(atoms)

    G = nx.Graph()
    G.add_nodes_from(range(len(atoms)))
    for i in range(len(atoms)):
        idxs, _ = nl.get_neighbors(i)
        for j in idxs:
            if i < j:
                G.add_edge(i, j)
    
    components = [sorted(list(c)) for c in nx.connected_components(G)]
    return components


def unwrap_to_single_image(atoms: Atoms, idxs: List[int]) -> np.ndarray:
    """
    Unwrap molecule to single periodic image (no bonds across boundaries).
    See MOL_UTILS_USAGE.md for details.
    """
    spos = atoms.get_scaled_positions()
    cell = atoms.get_cell()
    s = spos[idxs].copy()
    ref = s[0]
    
    # Shift each atom to be near the reference atom in fractional space
    for i in range(len(s)):
        d = s[i] - ref
        s[i] -= np.round(d)  # bring within [-0.5, 0.5)
    
    cart = s @ cell  # Convert to Cartesian (Å)
    return np.asarray(cart)


def recenter_system(unwrapped_positions: List[np.ndarray], cell: np.ndarray) -> List[np.ndarray]:
    """
    Recenter system to form compact cluster around box center.
    See MOL_UTILS_USAGE.md for details.
    """
    # Calculate system COM
    all_positions = np.vstack(unwrapped_positions)
    system_com = np.mean(all_positions, axis=0)
    
    # Box center in Cartesian coordinates
    box_center = 0.5 * (cell[0] + cell[1] + cell[2])
    
    # Calculate shift to center the system
    shift = box_center - system_com
    
    # Apply shift to all molecules
    recentered = []
    for mol_pos in unwrapped_positions:
        new_pos = mol_pos + shift
        recentered.append(new_pos)
    
    return recentered


def wrap_to_box(unwrapped_positions: List[np.ndarray], cell: np.ndarray) -> List[np.ndarray]:
    """
    Wrap molecules into simulation box (maintains molecular integrity).
    See MOL_UTILS_USAGE.md for details.
    """
    cell_inv = np.linalg.inv(cell)
    wrapped = []
    
    for mol_pos in unwrapped_positions:
        # Convert to fractional coordinates
        frac = mol_pos @ cell_inv
        
        # Find COM in fractional coordinates
        com_frac = np.mean(frac, axis=0)
        
        # Calculate integer shift to bring COM into [0, 1)
        shift_frac = -np.floor(com_frac)
        
        # Apply shift to all atoms in molecule
        frac_shifted = frac + shift_frac
        
        # Convert back to Cartesian
        cart_shifted = frac_shifted @ cell
        wrapped.append(cart_shifted)
    
    return wrapped


def deterministic_template_order(G: nx.Graph, atoms: Atoms, idxs: List[int]) -> List[int]:
    """
    Choose stable atom order for template molecule using BFS from unique root.
    See MOL_UTILS_USAGE.md for details.
    """
    # Build mapping: local index -> global index
    local_to_global = {li: gi for li, gi in enumerate(idxs)}
    
    # Choose root
    scores = []
    for li, gi in enumerate(idxs):
        z = int(atoms.numbers[gi])
        deg = int(G.degree[gi])
        hnb = sum(1 for j in G.neighbors(gi) if atoms.numbers[j] == 1)
        scores.append((-z, -deg, hnb, li))  # sort ascending
    
    root_local = sorted(range(len(idxs)), key=lambda k: scores[k])[0]
    root_global = idxs[root_local]

    # BFS from root
    order_global = list(nx.bfs_tree(G, source=root_global).nodes())
    # Filter to this molecule only and stabilize order
    order_global = [g for g in order_global if g in idxs]
    order_global.sort(key=lambda g: (0,))  # keep BFS order

    # Map to local 0..n-1 order
    inv = {gi: li for li, gi in enumerate(idxs)}
    order_local = [inv[g] for g in order_global]
    
    return order_local




def best_isomorphism(src: Atoms, ref: Atoms, src_idxs: List[int], ref_idxs: List[int], 
                    scale: float = 1.10) -> Dict[int, int]:
    """
    Find graph isomorphism mapping minimizing RMSD (Kabsch alignment).
    See MOL_UTILS_USAGE.md for details.
    
    WARNING: This is slow for large molecules (>50 atoms). Use best_isomorphism_fast instead.
    """
    # Build graphs *restricted to the molecule* (reindexed 0..n-1)
    def sub_atoms(atoms, idxs):
        return Atoms(numbers=atoms.numbers[idxs], positions=atoms.positions[idxs])

    sub_s = sub_atoms(src, src_idxs)
    sub_r = sub_atoms(ref, ref_idxs)

    Gs = build_graph(sub_s, scale)
    Gr = build_graph(sub_r, scale)

    def node_match(n1, n2):
        return (n1["Z"] == n2["Z"]) and (n1["deg"] == n2["deg"]) and (n1["hnb"] == n2["hnb"])

    GM = nx.algorithms.isomorphism.GraphMatcher(Gs, Gr, node_match=node_match)

    Ps = sub_s.get_positions()
    Pr = sub_r.get_positions()
    best_map, best_rmsd = None, 1e9
    
    for mapping in GM.isomorphisms_iter():
        # mapping: src_local -> ref_local
        # Kabsch RMSD in template order 0..n-1
        src_ordered = Ps[np.array(sorted(mapping.keys()), dtype=int)]
        ref_ordered = Pr[np.array([mapping[k] for k in sorted(mapping.keys())], dtype=int)]
        
        # Center
        src_c = src_ordered - src_ordered.mean(0)
        ref_c = ref_ordered - ref_ordered.mean(0)
        
        # Kabsch algorithm
        H = src_c.T @ ref_c
        U, S, Vt = np.linalg.svd(H)
        R = U @ Vt
        rot_src = src_c @ R
        
        rmsd = np.sqrt(((rot_src - ref_c)**2).sum() / len(src_c))
        
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_map = mapping
    
    if best_map is None:
        raise RuntimeError("Graph isomorphism failed; try adjusting --bond-scale or ensure H atoms are present.")
    
    return best_map


def weisfeiler_lehman_hash(G: nx.Graph, iterations: int = 3) -> Dict[int, str]:
    """
    Compute Weisfeiler-Lehman graph hash for each node.
    
    This creates a unique "fingerprint" for each atom based on its local 
    graph neighborhood, which can be used for fast matching.
    
    Args:
        G: NetworkX graph with node attributes 'Z', 'deg', 'hnb'
        iterations: Number of WL iterations (more = more discriminating)
    
    Returns:
        Dictionary mapping node index to hash string
    """
    # Initialize labels with node attributes
    labels = {}
    for node in G.nodes():
        attrs = G.nodes[node]
        labels[node] = f"{attrs['Z']}_{attrs['deg']}_{attrs['hnb']}"
    
    for _ in range(iterations):
        new_labels = {}
        for node in G.nodes():
            # Get sorted neighbor labels
            neighbor_labels = sorted([labels[n] for n in G.neighbors(node)])
            # Combine own label with neighbor labels
            new_labels[node] = labels[node] + "_[" + ",".join(neighbor_labels) + "]"
        labels = new_labels
    
    return labels


def best_isomorphism_fast(src: Atoms, ref: Atoms, src_idxs: List[int], ref_idxs: List[int], 
                          scale: float = 1.10) -> Dict[int, int]:
    """
    Fast graph isomorphism using Weisfeiler-Lehman hashing + coordinate matching.
    
    Algorithm:
    1. Compute WL hash for each atom (captures local graph structure)
    2. Group atoms by (element, WL hash)
    3. For atoms with same hash, match by nearest distance after alignment
    
    This is O(n log n) instead of O(n!) for the naive approach.
    Accurate for molecular crystals where molecules have identical topology.
    
    Args:
        src: Source molecule Atoms object
        ref: Reference (template) molecule Atoms object
        src_idxs: Indices of source atoms (should be 0..n-1 for extracted molecule)
        ref_idxs: Indices of reference atoms (should be 0..n-1 for extracted molecule)
        scale: Bond scale factor for graph building
    
    Returns:
        Dictionary mapping src index -> ref index
    """
    def sub_atoms(atoms, idxs):
        return Atoms(numbers=atoms.numbers[idxs], positions=atoms.positions[idxs])

    sub_s = sub_atoms(src, src_idxs)
    sub_r = sub_atoms(ref, ref_idxs)

    Gs = build_graph(sub_s, scale)
    Gr = build_graph(sub_r, scale)
    
    # Compute WL hashes
    hash_s = weisfeiler_lehman_hash(Gs)
    hash_r = weisfeiler_lehman_hash(Gr)
    
    Ps = sub_s.get_positions()
    Pr = sub_r.get_positions()
    
    # First, do Kabsch alignment of src to ref (using all atoms)
    Ps_centered = Ps - Ps.mean(axis=0)
    Pr_centered = Pr - Pr.mean(axis=0)
    
    # Kabsch rotation matrix
    H = Ps_centered.T @ Pr_centered
    U, S, Vt = np.linalg.svd(H)
    R = U @ Vt
    # Handle reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = U @ Vt
    
    # Align source positions to reference frame
    Ps_aligned = Ps_centered @ R
    
    # Group reference atoms by their WL hash
    ref_by_hash: Dict[str, List[int]] = {}
    for i, h in hash_r.items():
        if h not in ref_by_hash:
            ref_by_hash[h] = []
        ref_by_hash[h].append(i)
    
    # Match each source atom to reference atom with same hash and nearest position
    mapping = {}
    used_ref = set()
    
    for src_i in range(len(sub_s)):
        src_hash = hash_s[src_i]
        
        if src_hash not in ref_by_hash:
            raise RuntimeError(f"No matching hash found for atom {src_i}. "
                             "Molecules may have different topology.")
        
        # Find nearest unmatched reference atom with same hash
        candidates = [r for r in ref_by_hash[src_hash] if r not in used_ref]
        
        if not candidates:
            raise RuntimeError(f"No available match for atom {src_i} with hash {src_hash}. "
                             "This shouldn't happen if molecules are identical.")
        
        # Find nearest by position
        src_pos = Ps_aligned[src_i]
        best_ref = None
        best_dist = float('inf')
        
        for ref_i in candidates:
            dist = np.linalg.norm(src_pos - Pr_centered[ref_i])
            if dist < best_dist:
                best_dist = dist
                best_ref = ref_i
        
        mapping[src_i] = best_ref
        used_ref.add(best_ref)
    
    return mapping

