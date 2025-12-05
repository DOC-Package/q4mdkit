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
