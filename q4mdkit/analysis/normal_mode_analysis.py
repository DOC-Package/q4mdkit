#!/usr/bin/env python3
"""
Normal Mode Coordinate Analysis from MD Trajectory

This module implements the axis-switching and normal mode projection method 
for analyzing molecular dynamics trajectories in terms of reference normal modes.

The method consists of:
1. Axis switching (rigid-body alignment via QCP algorithm using MDAnalysis.lib.qcprot)
2. Mass-weighted displacement calculation
3. Projection onto reference normal modes

All coordinates are handled in Angstrom units internally.

Reference:
- Horne & Li (1998) IDJ
- Coutsias et al. (2004) JCC
- Theobald (2005) ACA
- Liu et al. (2010) JCC

Author: q4mdkit
"""

import numpy as np
from typing import Tuple, Optional, List, Union
import warnings

# MDAnalysis QCP rotation
from MDAnalysis.lib.qcprot import CalcRMSDRotationalMatrix

# Unit conversion constants
ANGSTROM_TO_BOHR = 1.8897259886
BOHR_TO_ANGSTROM = 1.0 / ANGSTROM_TO_BOHR


# Atomic masses (element -> mass [amu])
# Using values from DFTB+ SK files (3ob)
ATOMIC_MASSES = {
    'H': 1.008,
    'He': 4.002602,
    'Li': 6.941,
    'Be': 9.012182,
    'B': 10.811,
    'C': 12.01,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984032,
    'Ne': 20.1797,
    'Na': 22.98976928,
    'Mg': 24.3050,
    'Al': 26.9815386,
    'Si': 28.0855,
    'P': 30.973762,
    'S': 32.065,
    'Cl': 35.453,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
}


def read_gen_file(filename: str, unit: str = 'angstrom') -> Tuple[np.ndarray, List[str]]:
    """
    Read DFTB+ .gen file format.
    
    Parameters
    ----------
    filename : str
        Path to .gen file
    unit : str
        Output unit: 'angstrom' (default) or 'bohr'
        Note: The .gen file is assumed to be in Angstrom units.
    
    Returns
    -------
    coords : np.ndarray
        Atomic coordinates (shape: (N, 3))
    atom_types : list
        List of element symbols for each atom
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Line 1: Number of atoms and coordinate type
    first_line = lines[0].split()
    n_atoms = int(first_line[0])
    coord_type = first_line[1].upper()  # 'C' for cluster, 'S' for supercell, 'F' for fractional
    
    # Line 2: List of element types
    element_types = lines[1].split()
    
    # Lines 3 onwards: Atom data (index, type_index, x, y, z)
    coords = np.zeros((n_atoms, 3))
    atom_types = []
    
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        type_idx = int(parts[1]) - 1  # 1-indexed -> 0-indexed
        atom_types.append(element_types[type_idx])
        coords[i-2, 0] = float(parts[2])
        coords[i-2, 1] = float(parts[3])
        coords[i-2, 2] = float(parts[4])
    
    # Convert units if needed
    if unit.lower() == 'bohr':
        coords = coords * ANGSTROM_TO_BOHR
    
    return coords, atom_types


def read_hessian_eigenvectors(filename: str) -> np.ndarray:
    """
    Read eigenvectors from hessian_eigenvectors.txt file.
    
    The eigenvectors are mass-weighted normal modes u_k.
    Each column is one eigenvector (ascending eigenvalue order).
    
    Parameters
    ----------
    filename : str
        Path to eigenvectors file
    
    Returns
    -------
    eigenvectors : np.ndarray
        Eigenvector matrix (3N x 3N), columns are eigenvectors
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            values = [float(x) for x in line.split()]
            data.append(values)
    
    eigenvectors = np.array(data)
    return eigenvectors


def read_hessian_eigenvalues(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read eigenvalues from hessian_eigenvalues.txt file.
    
    Parameters
    ----------
    filename : str
        Path to eigenvalues file
    
    Returns
    -------
    eigenvalues : np.ndarray
        Eigenvalues (3N)
    frequencies_cm : np.ndarray
        Frequencies in cm^-1 (3N)
    """
    eigenvalues = []
    frequencies_cm = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            # Format: index, eigenvalue, freq_hartree, freq_cm, type
            eigenvalues.append(float(parts[1]))
            frequencies_cm.append(float(parts[3]))
    
    return np.array(eigenvalues), np.array(frequencies_cm)


def get_masses(atom_types: List[str]) -> np.ndarray:
    """
    Get atomic masses from atom types.
    
    Parameters
    ----------
    atom_types : list
        List of element symbols
    
    Returns
    -------
    masses : np.ndarray
        Atomic masses in amu (shape: (N,))
    """
    return np.array([ATOMIC_MASSES[atom] for atom in atom_types])


def center_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Calculate center of mass.
    
    Parameters
    ----------
    coords : np.ndarray
        Atomic coordinates (N, 3)
    masses : np.ndarray
        Atomic masses (N,)
    
    Returns
    -------
    com : np.ndarray
        Center of mass (3,)
    """
    M = np.sum(masses)
    return np.sum(masses[:, np.newaxis] * coords, axis=0) / M


def remove_center_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Remove rigid translation by shifting to center of mass.
    
    r_bar_a = R_a - (1/M) * sum_a(m_a * R_a)
    
    Parameters
    ----------
    coords : np.ndarray
        Atomic coordinates (N, 3) or (n_frames, N, 3)
    masses : np.ndarray
        Atomic masses (N,)
    
    Returns
    -------
    coords_centered : np.ndarray
        Centered coordinates (same shape as input)
    """
    if coords.ndim == 2:
        com = center_of_mass(coords, masses)
        return coords - com
    elif coords.ndim == 3:
        # Multiple frames
        n_frames = coords.shape[0]
        centered = np.zeros_like(coords)
        for i in range(n_frames):
            com = center_of_mass(coords[i], masses)
            centered[i] = coords[i] - com
        return centered
    else:
        raise ValueError("coords must be 2D (N, 3) or 3D (n_frames, N, 3)")


def qcp_rotation(P: np.ndarray, Q: np.ndarray, 
                 weights: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
    """
    Calculate optimal rotation matrix using QCP algorithm via MDAnalysis.lib.qcprot.
    
    Uses the Quaternion Characteristic Polynomial (QCP) algorithm which is 
    faster and more numerically stable than SVD-based Kabsch.
    
    Finds rotation R that minimizes:
        sum_a w_a ||R @ P_a - Q_a||^2
    
    This aligns P onto Q (rotates P into Q's frame).
    
    Parameters
    ----------
    P : np.ndarray
        Source coordinates to be rotated (N, 3), should be centered
    Q : np.ndarray
        Target coordinates (N, 3), should be centered
    weights : np.ndarray, optional
        Weights for each atom (N,), default is uniform weights
    
    Returns
    -------
    R : np.ndarray
        Rotation matrix (3, 3) in SO(3)
    rmsd : float
        RMSD after optimal alignment
    """
    n_atoms = P.shape[0]
    
    # QCP requires contiguous float64 arrays
    P_c = np.ascontiguousarray(P, dtype=np.float64)
    Q_c = np.ascontiguousarray(Q, dtype=np.float64)
    
    # Prepare rotation matrix output (flat array of 9 elements)
    rot_flat = np.zeros(9, dtype=np.float64)
    
    # Prepare weights (None for uniform weights)
    if weights is not None:
        weights_c = np.ascontiguousarray(weights, dtype=np.float64)
    else:
        weights_c = None
    
    # Call QCP algorithm
    # CalcRMSDRotationalMatrix(ref, conf, n_atoms, rot, weights)
    # Returns RMSD, rot is modified in place
    rmsd = CalcRMSDRotationalMatrix(Q_c, P_c, n_atoms, rot_flat, weights_c)
    
    # Reshape to 3x3 rotation matrix
    R = rot_flat.reshape(3, 3)
    
    return R, rmsd








def kabsch_rotation(P: np.ndarray, Q: np.ndarray, 
                    weights: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Calculate optimal rotation matrix using Kabsch/QCP algorithm.
    
    Uses MDAnalysis.lib.qcprot if available (QCP algorithm), otherwise falls back to SVD.
    
    Finds rotation R that minimizes:
        sum_a w_a ||R @ P_a - Q_a||^2
    
    This aligns P onto Q (rotates P into Q's frame).
    
    Parameters
    ----------
    P : np.ndarray
        Source coordinates to be rotated (N, 3), should be centered
    Q : np.ndarray
        Target coordinates (N, 3), should be centered
    weights : np.ndarray, optional
        Weights for each atom (N,), default is uniform weights
    
    Returns
    -------
    R : np.ndarray
        Rotation matrix (3, 3) in SO(3)
    """
    R, _ = qcp_rotation(P, Q, weights)
    return R


def align_to_reference(coords: np.ndarray, ref_coords: np.ndarray, 
                       masses: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Align MD frame to reference structure using mass-weighted Kabsch algorithm.
    
    This implements the axis-switching procedure:
    R(t) = argmin_{R in SO(3)} sum_a m_a ||R @ r_bar_a(t) - r_bar_a^(0)||^2
    
    Parameters
    ----------
    coords : np.ndarray
        Current frame coordinates (N, 3), centered at COM
    ref_coords : np.ndarray
        Reference coordinates (N, 3), centered at COM
    masses : np.ndarray
        Atomic masses (N,)
    
    Returns
    -------
    aligned_coords : np.ndarray
        Aligned coordinates R @ coords (N, 3)
    R : np.ndarray
        Rotation matrix (3, 3)
    """
    # Use masses as weights
    R = kabsch_rotation(coords, ref_coords, weights=masses)
    aligned_coords = (R @ coords.T).T
    return aligned_coords, R


def calculate_displacements(coords: np.ndarray, ref_coords: np.ndarray,
                            masses: np.ndarray, 
                            align: bool = True) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Calculate displacements from reference structure.
    
    Implements:
    1. Remove rigid translation (center of mass)
    2. Remove rigid rotation (Kabsch alignment) if align=True
    3. Calculate displacement: delta_r_a(t) = R(t) @ r_bar_a(t) - r_bar_a^(0)
    
    Parameters
    ----------
    coords : np.ndarray
        Trajectory coordinates (n_frames, N, 3) or single frame (N, 3)
    ref_coords : np.ndarray
        Reference coordinates (N, 3)
    masses : np.ndarray
        Atomic masses (N,)
    align : bool
        Whether to perform Kabsch alignment (default: True)
    
    Returns
    -------
    displacements : np.ndarray
        Displacements (n_frames, N, 3) or (N, 3)
    rotations : np.ndarray or None
        Rotation matrices (n_frames, 3, 3) or (3, 3) if align=True
    """
    single_frame = coords.ndim == 2
    if single_frame:
        coords = coords[np.newaxis, :, :]
    
    n_frames = coords.shape[0]
    n_atoms = coords.shape[1]
    
    # Center reference structure
    ref_centered = remove_center_of_mass(ref_coords, masses)
    
    # Center trajectory frames
    coords_centered = remove_center_of_mass(coords, masses)
    
    displacements = np.zeros((n_frames, n_atoms, 3))
    rotations = np.zeros((n_frames, 3, 3)) if align else None
    
    for i in range(n_frames):
        if align:
            aligned, R = align_to_reference(coords_centered[i], ref_centered, masses)
            displacements[i] = aligned - ref_centered
            rotations[i] = R
        else:
            displacements[i] = coords_centered[i] - ref_centered
    
    if single_frame:
        displacements = displacements[0]
        if rotations is not None:
            rotations = rotations[0]
    
    return displacements, rotations


def mass_weight_displacements(displacements: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Convert displacements to mass-weighted displacements.
    
    delta_x(t) = M^(1/2) @ delta_r(t)
    
    Parameters
    ----------
    displacements : np.ndarray
        Cartesian displacements (n_frames, N, 3) or (N, 3)
    masses : np.ndarray
        Atomic masses (N,)
    
    Returns
    -------
    mw_displacements : np.ndarray
        Mass-weighted displacements (n_frames, 3N) or (3N,)
    """
    single_frame = displacements.ndim == 2
    if single_frame:
        displacements = displacements[np.newaxis, :, :]
    
    n_frames = displacements.shape[0]
    n_atoms = displacements.shape[1]
    
    # Expand masses to 3N
    sqrt_masses_3n = np.sqrt(np.repeat(masses, 3))
    
    # Flatten displacements to (n_frames, 3N)
    disp_flat = displacements.reshape(n_frames, 3 * n_atoms)
    
    # Apply mass weighting
    mw_disp = disp_flat * sqrt_masses_3n
    
    if single_frame:
        return mw_disp[0]
    return mw_disp


def project_onto_normal_modes(mw_displacements: np.ndarray, 
                               eigenvectors: np.ndarray,
                               mode_indices: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Project mass-weighted displacements onto normal modes.
    
    d_k(t) = u_k^(0)^T @ delta_x(t)
    
    Parameters
    ----------
    mw_displacements : np.ndarray
        Mass-weighted displacements (n_frames, 3N) or (3N,)
    eigenvectors : np.ndarray
        Normal mode eigenvectors (3N, 3N), columns are modes
    mode_indices : np.ndarray, optional
        Indices of modes to project onto (default: all modes except first 6)
    
    Returns
    -------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes) or (n_modes,)
    """
    n_total = eigenvectors.shape[0]
    
    if mode_indices is None:
        # Skip first 6 modes (rotational/translational)
        mode_indices = np.arange(6, n_total)
    
    single_frame = mw_displacements.ndim == 1
    if single_frame:
        mw_displacements = mw_displacements[np.newaxis, :]
    
    n_frames = mw_displacements.shape[0]
    n_modes = len(mode_indices)
    
    # Select modes
    modes = eigenvectors[:, mode_indices]  # (3N, n_modes)
    
    # Project: d_k = u_k^T @ delta_x
    mode_coords = mw_displacements @ modes  # (n_frames, n_modes)
    
    if single_frame:
        return mode_coords[0]
    return mode_coords


def analyze_trajectory(trajectory_coords: np.ndarray,
                       ref_coords: np.ndarray,
                       masses: np.ndarray,
                       eigenvectors: np.ndarray,
                       mode_indices: Optional[np.ndarray] = None,
                       align: bool = True) -> dict:
    """
    Full normal mode analysis of MD trajectory.
    
    Parameters
    ----------
    trajectory_coords : np.ndarray
        Trajectory coordinates (n_frames, N, 3)
    ref_coords : np.ndarray
        Reference (optimized) coordinates (N, 3)
    masses : np.ndarray
        Atomic masses (N,)
    eigenvectors : np.ndarray
        Normal mode eigenvectors (3N, 3N)
    mode_indices : np.ndarray, optional
        Indices of modes to analyze (default: skip first 6)
    align : bool
        Whether to perform Kabsch alignment
    
    Returns
    -------
    results : dict
        Dictionary containing:
        - 'mode_coords': Normal mode coordinates (n_frames, n_modes)
        - 'displacements': Cartesian displacements (n_frames, N, 3)
        - 'rotations': Rotation matrices (n_frames, 3, 3) if align=True
        - 'rmsd': RMSD from reference for each frame
    """
    n_total = eigenvectors.shape[0]
    
    if mode_indices is None:
        mode_indices = np.arange(6, n_total)
    
    # Calculate displacements
    displacements, rotations = calculate_displacements(
        trajectory_coords, ref_coords, masses, align=align
    )
    
    # Mass-weight displacements
    mw_displacements = mass_weight_displacements(displacements, masses)
    
    # Project onto normal modes
    mode_coords = project_onto_normal_modes(mw_displacements, eigenvectors, mode_indices)
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(displacements**2, axis=(1, 2)))
    
    results = {
        'mode_coords': mode_coords,
        'displacements': displacements,
        'rotations': rotations,
        'rmsd': rmsd,
        'mode_indices': mode_indices,
    }
    
    return results


def calculate_mode_statistics(mode_coords: np.ndarray, 
                               frequencies_cm: np.ndarray,
                               temperature: float = 300.0) -> dict:
    """
    Calculate statistical properties of normal mode coordinates.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    frequencies_cm : np.ndarray
        Frequencies in cm^-1 for the selected modes (n_modes,)
    temperature : float
        Temperature in Kelvin (default: 300 K)
    
    Returns
    -------
    stats : dict
        Dictionary containing:
        - 'mean': Mean of mode coordinates
        - 'std': Standard deviation
        - 'variance': Variance <d_k^2>
        - 'expected_variance': Classical expected variance kT/(m*omega^2)
    """
    # Constants
    kB = 3.1668114e-6  # Hartree/K
    cm_to_hartree = 1.0 / 219474.63137
    
    mean = np.mean(mode_coords, axis=0)
    std = np.std(mode_coords, axis=0)
    variance = np.var(mode_coords, axis=0)
    
    # Expected classical variance: <d_k^2> = kT / omega_k^2
    # (in mass-weighted coordinates, omega is angular frequency)
    omega_hartree = frequencies_cm * cm_to_hartree
    expected_variance = np.zeros_like(frequencies_cm)
    nonzero = np.abs(omega_hartree) > 1e-10
    expected_variance[nonzero] = kB * temperature / (omega_hartree[nonzero]**2)
    
    stats = {
        'mean': mean,
        'std': std,
        'variance': variance,
        'expected_variance': expected_variance,
    }
    
    return stats


def calculate_mode_correlation(mode_coords: np.ndarray,
                                mode_indices: Optional[np.ndarray] = None) -> dict:
    """
    Calculate correlation matrix of normal mode coordinates.
    
    Computes Pearson correlation coefficient r_kl = Cov(d_k, d_l) / (std_k * std_l)
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    mode_indices : np.ndarray, optional
        Mode indices for labeling (default: 7, 8, 9, ...)
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'correlation_matrix': Full correlation matrix (n_modes, n_modes)
        - 'mode_indices': Mode indices used for labeling
        - 'strong_correlations': List of (mode_i, mode_j, r) for |r| > 0.5
        - 'statistics': Dict with max, mean, median of off-diagonal |r|
    """
    n_modes = mode_coords.shape[1]
    
    if mode_indices is None:
        mode_indices = np.arange(7, 7 + n_modes)
    
    # Pearson correlation matrix
    corr_matrix = np.corrcoef(mode_coords.T)
    
    # Extract off-diagonal elements
    upper_tri_idx = np.triu_indices(n_modes, k=1)
    off_diag = corr_matrix[upper_tri_idx]
    
    # Find strong correlations (|r| > 0.5)
    strong_correlations = []
    for idx in range(len(off_diag)):
        r = off_diag[idx]
        if np.abs(r) > 0.5:
            i = upper_tri_idx[0][idx]
            j = upper_tri_idx[1][idx]
            strong_correlations.append((mode_indices[i], mode_indices[j], r))
    
    # Sort by absolute correlation
    strong_correlations.sort(key=lambda x: -np.abs(x[2]))
    
    # Statistics
    statistics = {
        'max_abs': np.max(np.abs(off_diag)),
        'mean_abs': np.mean(np.abs(off_diag)),
        'median_abs': np.median(np.abs(off_diag)),
        'std_abs': np.std(np.abs(off_diag)),
    }
    
    return {
        'correlation_matrix': corr_matrix,
        'mode_indices': mode_indices,
        'strong_correlations': strong_correlations,
        'statistics': statistics,
    }


def plot_mode_correlation(mode_coords: np.ndarray,
                          mode_indices: Optional[np.ndarray] = None,
                          output_file: Optional[str] = None,
                          figsize: Tuple[int, int] = (10, 8),
                          cmap: str = 'RdBu_r',
                          vmin: float = -1.0,
                          vmax: float = 1.0,
                          show_values: bool = False,
                          title: Optional[str] = None) -> 'matplotlib.figure.Figure':
    """
    Plot correlation matrix of normal mode coordinates as a heatmap.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    mode_indices : np.ndarray, optional
        Mode indices for axis labels (default: 7, 8, 9, ...)
    output_file : str, optional
        Path to save the figure (PNG, PDF, etc.)
    figsize : tuple
        Figure size (width, height) in inches
    cmap : str
        Colormap name (default: 'RdBu_r' for diverging red-blue)
    vmin, vmax : float
        Color scale limits (default: -1 to 1)
    show_values : bool
        Whether to annotate cells with correlation values
    title : str, optional
        Plot title
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    
    # Calculate correlation
    result = calculate_mode_correlation(mode_coords, mode_indices)
    corr_matrix = result['correlation_matrix']
    mode_idx = result['mode_indices']
    stats = result['statistics']
    
    n_modes = len(mode_idx)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    im = ax.imshow(corr_matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Pearson correlation $r$', fontsize=12)
    
    # Axis labels
    if n_modes <= 30:
        # Show all tick labels for small matrices
        ax.set_xticks(np.arange(n_modes))
        ax.set_yticks(np.arange(n_modes))
        ax.set_xticklabels(mode_idx, fontsize=8)
        ax.set_yticklabels(mode_idx, fontsize=8)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    else:
        # Show subset of tick labels for large matrices
        step = max(1, n_modes // 10)
        tick_positions = np.arange(0, n_modes, step)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels(mode_idx[tick_positions], fontsize=10)
        ax.set_yticklabels(mode_idx[tick_positions], fontsize=10)
    
    ax.set_xlabel('Mode index', fontsize=12)
    ax.set_ylabel('Mode index', fontsize=12)
    
    # Title
    if title is None:
        title = f'Mode Correlation Matrix\nmax|r|={stats["max_abs"]:.3f}, mean|r|={stats["mean_abs"]:.3f}'
    ax.set_title(title, fontsize=14)
    
    # Annotate values if requested
    if show_values and n_modes <= 20:
        for i in range(n_modes):
            for j in range(n_modes):
                val = corr_matrix[i, j]
                color = 'white' if np.abs(val) > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                       color=color, fontsize=6)
    
    plt.tight_layout()
    
    if output_file is not None:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved correlation plot to {output_file}")
    
    return fig


def plot_strong_correlations(mode_coords: np.ndarray,
                              mode_indices: Optional[np.ndarray] = None,
                              threshold: float = 0.5,
                              max_pairs: int = 12,
                              output_file: Optional[str] = None,
                              figsize: Optional[Tuple[int, int]] = None) -> 'matplotlib.figure.Figure':
    """
    Plot scatter plots of strongly correlated mode pairs.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    mode_indices : np.ndarray, optional
        Mode indices for labeling
    threshold : float
        Minimum |r| to include (default: 0.5)
    max_pairs : int
        Maximum number of pairs to plot (default: 12)
    output_file : str, optional
        Path to save the figure
    figsize : tuple, optional
        Figure size (default: auto-calculated)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    
    # Calculate correlations
    result = calculate_mode_correlation(mode_coords, mode_indices)
    strong_corr = result['strong_correlations']
    mode_idx = result['mode_indices']
    
    # Filter by threshold
    pairs = [(i, j, r) for (i, j, r) in strong_corr if np.abs(r) >= threshold][:max_pairs]
    
    if len(pairs) == 0:
        print(f"No correlations with |r| >= {threshold}")
        return None
    
    n_pairs = len(pairs)
    n_cols = min(4, n_pairs)
    n_rows = (n_pairs + n_cols - 1) // n_cols
    
    if figsize is None:
        figsize = (4 * n_cols, 4 * n_rows)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_pairs == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Map mode index to column index
    mode_to_col = {m: i for i, m in enumerate(mode_idx)}
    
    for idx, (mode_i, mode_j, r) in enumerate(pairs):
        ax = axes[idx]
        col_i = mode_to_col[mode_i]
        col_j = mode_to_col[mode_j]
        
        x = mode_coords[:, col_i]
        y = mode_coords[:, col_j]
        
        ax.scatter(x, y, alpha=0.5, s=10, c='steelblue')
        ax.set_xlabel(f'$d_{{{mode_i}}}$', fontsize=11)
        ax.set_ylabel(f'$d_{{{mode_j}}}$', fontsize=11)
        ax.set_title(f'Mode {mode_i}-{mode_j}: r={r:.3f}', fontsize=11)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
        
        # Add regression line
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, p(x_line), 'r-', linewidth=1.5, alpha=0.7)
    
    # Hide empty subplots
    for idx in range(n_pairs, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(f'Strongly Correlated Mode Pairs (|r| ≥ {threshold})', fontsize=14, y=1.02)
    plt.tight_layout()
    
    if output_file is not None:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved scatter plots to {output_file}")
    
    return fig


def print_correlation_summary(mode_coords: np.ndarray,
                               mode_indices: Optional[np.ndarray] = None,
                               top_n: int = 10) -> None:
    """
    Print summary of mode correlations to console.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    mode_indices : np.ndarray, optional
        Mode indices for labeling
    top_n : int
        Number of top correlations to show
    """
    result = calculate_mode_correlation(mode_coords, mode_indices)
    stats = result['statistics']
    strong_corr = result['strong_correlations']
    
    print("=" * 60)
    print("Normal Mode Correlation Summary")
    print("=" * 60)
    print(f"Number of modes: {len(result['mode_indices'])}")
    print(f"Number of mode pairs: {len(result['mode_indices']) * (len(result['mode_indices'])-1) // 2}")
    print()
    print("Off-diagonal correlation statistics:")
    print(f"  Max |r|:    {stats['max_abs']:.6f}")
    print(f"  Mean |r|:   {stats['mean_abs']:.6f}")
    print(f"  Median |r|: {stats['median_abs']:.6f}")
    print(f"  Std |r|:    {stats['std_abs']:.6f}")
    print()
    
    if strong_corr:
        print(f"Top {min(top_n, len(strong_corr))} correlated mode pairs:")
        print("-" * 40)
        print(f"{'Mode i':>8}  {'Mode j':>8}  {'r':>10}")
        print("-" * 40)
        for mode_i, mode_j, r in strong_corr[:top_n]:
            print(f"{mode_i:>8}  {mode_j:>8}  {r:>10.6f}")
    else:
        print("No strong correlations (|r| > 0.5) found.")
    print("=" * 60)


def plot_mode_time_evolution(mode_coords: np.ndarray,
                              selected_modes: Union[List[int], np.ndarray],
                              mode_indices: Optional[np.ndarray] = None,
                              dt: float = 1.0,
                              time_unit: str = 'frame',
                              output_file: Optional[str] = None,
                              figsize: Optional[Tuple[int, int]] = None,
                              colors: Optional[List[str]] = None,
                              title: Optional[str] = None,
                              show_mean: bool = True,
                              share_y: bool = False) -> 'matplotlib.figure.Figure':
    """
    Plot time evolution of selected normal mode coordinates.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    selected_modes : list or array
        Mode indices to plot (1-indexed if mode_indices starts from 7)
    mode_indices : np.ndarray, optional
        Mode indices array for column mapping (default: 7, 8, 9, ...)
    dt : float
        Time step between frames (default: 1.0)
    time_unit : str
        Time unit label for x-axis (default: 'frame')
    output_file : str, optional
        Path to save the figure
    figsize : tuple, optional
        Figure size (width, height) in inches
    colors : list, optional
        Colors for each mode line
    title : str, optional
        Plot title
    show_mean : bool
        Whether to show mean value as horizontal line (default: True)
    share_y : bool
        Whether to share y-axis across subplots (default: False)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    
    n_frames = mode_coords.shape[0]
    n_modes_total = mode_coords.shape[1]
    
    if mode_indices is None:
        mode_indices = np.arange(7, 7 + n_modes_total)
    
    # Map mode labels to column indices
    mode_to_col = {m: i for i, m in enumerate(mode_indices)}
    
    # Validate selected modes
    valid_modes = []
    for m in selected_modes:
        if m in mode_to_col:
            valid_modes.append(m)
        else:
            print(f"Warning: Mode {m} not found in mode_indices, skipping")
    
    if len(valid_modes) == 0:
        print("Error: No valid modes to plot")
        return None
    
    n_selected = len(valid_modes)
    
    # Time axis
    time = np.arange(n_frames) * dt
    
    # Figure setup
    if figsize is None:
        figsize = (12, 3 * n_selected)
    
    fig, axes = plt.subplots(n_selected, 1, figsize=figsize, sharex=True, sharey=share_y)
    if n_selected == 1:
        axes = [axes]
    
    # Colors
    if colors is None:
        cmap = plt.cm.tab10
        colors = [cmap(i % 10) for i in range(n_selected)]
    
    for idx, mode in enumerate(valid_modes):
        ax = axes[idx]
        col = mode_to_col[mode]
        data = mode_coords[:, col]
        
        ax.plot(time, data, color=colors[idx], linewidth=0.8, alpha=0.8)
        
        if show_mean:
            mean_val = np.mean(data)
            ax.axhline(mean_val, color=colors[idx], linestyle='--', 
                      linewidth=1.5, alpha=0.7, label=f'mean={mean_val:.3f}')
            ax.axhline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        
        ax.set_ylabel(f'$d_{{{mode}}}$', fontsize=12)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # Add statistics annotation
        std_val = np.std(data)
        ax.text(0.02, 0.95, f'std={std_val:.3f}', transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    axes[-1].set_xlabel(f'Time ({time_unit})', fontsize=12)
    
    if title is None:
        title = f'Normal Mode Time Evolution (modes: {", ".join(map(str, valid_modes))})'
    fig.suptitle(title, fontsize=14, y=1.01)
    
    plt.tight_layout()
    
    if output_file is not None:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved time evolution plot to {output_file}")
    
    return fig


def plot_mode_time_evolution_combined(mode_coords: np.ndarray,
                                       selected_modes: Union[List[int], np.ndarray],
                                       mode_indices: Optional[np.ndarray] = None,
                                       dt: float = 1.0,
                                       time_unit: str = 'frame',
                                       output_file: Optional[str] = None,
                                       figsize: Tuple[int, int] = (12, 6),
                                       normalize: bool = False,
                                       title: Optional[str] = None) -> 'matplotlib.figure.Figure':
    """
    Plot time evolution of multiple modes on a single plot.
    
    Parameters
    ----------
    mode_coords : np.ndarray
        Normal mode coordinates (n_frames, n_modes)
    selected_modes : list or array
        Mode indices to plot
    mode_indices : np.ndarray, optional
        Mode indices array for column mapping
    dt : float
        Time step between frames
    time_unit : str
        Time unit label for x-axis
    output_file : str, optional
        Path to save the figure
    figsize : tuple
        Figure size
    normalize : bool
        Whether to normalize each mode by its std (default: False)
    title : str, optional
        Plot title
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    
    n_frames = mode_coords.shape[0]
    n_modes_total = mode_coords.shape[1]
    
    if mode_indices is None:
        mode_indices = np.arange(7, 7 + n_modes_total)
    
    mode_to_col = {m: i for i, m in enumerate(mode_indices)}
    
    # Validate selected modes
    valid_modes = [m for m in selected_modes if m in mode_to_col]
    if len(valid_modes) == 0:
        print("Error: No valid modes to plot")
        return None
    
    # Time axis
    time = np.arange(n_frames) * dt
    
    fig, ax = plt.subplots(figsize=figsize)
    
    cmap = plt.cm.tab10
    for idx, mode in enumerate(valid_modes):
        col = mode_to_col[mode]
        data = mode_coords[:, col]
        
        if normalize:
            data = (data - np.mean(data)) / np.std(data)
            ylabel = 'Normalized $d_k$'
        else:
            ylabel = '$d_k$ (mass-weighted)'
        
        ax.plot(time, data, label=f'Mode {mode}', 
               color=cmap(idx % 10), linewidth=0.8, alpha=0.8)
    
    ax.axhline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    ax.set_xlabel(f'Time ({time_unit})', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(loc='upper right', fontsize=10, ncol=min(4, len(valid_modes)))
    ax.grid(True, alpha=0.3)
    
    if title is None:
        title = 'Normal Mode Time Evolution'
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    
    if output_file is not None:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved combined time evolution plot to {output_file}")
    
    return fig


class NormalModeAnalyzer:
    """
    Class for normal mode coordinate analysis of MD trajectories.
    
    This implements the axis-switching and normal mode projection procedure.
    
    Parameters
    ----------
    ref_coords : np.ndarray
        Reference structure coordinates (N, 3)
    atom_types : list
        List of element symbols
    eigenvectors : np.ndarray
        Normal mode eigenvectors (3N, 3N)
    eigenvalues : np.ndarray, optional
        Eigenvalues (3N,)
    frequencies_cm : np.ndarray, optional
        Frequencies in cm^-1 (3N,)
    mode_range : tuple, optional
        Range of modes to analyze (start, end), default (6, 3N)
    
    Attributes
    ----------
    masses : np.ndarray
        Atomic masses
    n_atoms : int
        Number of atoms
    n_modes : int
        Number of vibrational modes
    """
    
    def __init__(self, 
                 ref_coords: np.ndarray,
                 atom_types: List[str],
                 eigenvectors: np.ndarray,
                 eigenvalues: Optional[np.ndarray] = None,
                 frequencies_cm: Optional[np.ndarray] = None,
                 mode_range: Optional[Tuple[int, int]] = None):
        
        self.ref_coords = ref_coords
        self.atom_types = atom_types
        self.eigenvectors = eigenvectors
        self.eigenvalues = eigenvalues
        self.frequencies_cm = frequencies_cm
        
        self.masses = get_masses(atom_types)
        self.n_atoms = len(atom_types)
        self.n_dof = 3 * self.n_atoms
        
        if mode_range is None:
            mode_range = (6, self.n_dof)
        self.mode_range = mode_range
        self.mode_indices = np.arange(mode_range[0], mode_range[1])
        self.n_modes = len(self.mode_indices)
        
        # Center reference structure
        self._ref_centered = remove_center_of_mass(ref_coords, self.masses)
    
    @classmethod
    def from_files(cls, 
                   gen_file: str,
                   eigenvectors_file: str,
                   eigenvalues_file: Optional[str] = None,
                   mode_range: Optional[Tuple[int, int]] = None) -> 'NormalModeAnalyzer':
        """
        Create analyzer from DFTB+ output files.
        
        Parameters
        ----------
        gen_file : str
            Path to optimized structure .gen file
        eigenvectors_file : str
            Path to hessian_eigenvectors.txt
        eigenvalues_file : str, optional
            Path to hessian_eigenvalues.txt
        mode_range : tuple, optional
            Range of modes to analyze
        
        Returns
        -------
        analyzer : NormalModeAnalyzer
        """
        ref_coords, atom_types = read_gen_file(gen_file)
        eigenvectors = read_hessian_eigenvectors(eigenvectors_file)
        
        eigenvalues = None
        frequencies_cm = None
        if eigenvalues_file is not None:
            eigenvalues, frequencies_cm = read_hessian_eigenvalues(eigenvalues_file)
        
        return cls(ref_coords, atom_types, eigenvectors, 
                   eigenvalues, frequencies_cm, mode_range)
    
    def analyze_frame(self, coords: np.ndarray, align: bool = True) -> dict:
        """
        Analyze a single MD frame.
        
        Parameters
        ----------
        coords : np.ndarray
            Frame coordinates (N, 3)
        align : bool
            Whether to perform Kabsch alignment
        
        Returns
        -------
        result : dict
            Analysis results for single frame
        """
        displacements, rotation = calculate_displacements(
            coords, self.ref_coords, self.masses, align=align
        )
        mw_disp = mass_weight_displacements(displacements, self.masses)
        mode_coords = project_onto_normal_modes(mw_disp, self.eigenvectors, self.mode_indices)
        
        rmsd = np.sqrt(np.mean(displacements**2))
        
        return {
            'mode_coords': mode_coords,
            'displacement': displacements,
            'rotation': rotation,
            'rmsd': rmsd,
        }
    
    def analyze_trajectory(self, trajectory_coords: np.ndarray, 
                           align: bool = True) -> dict:
        """
        Analyze full MD trajectory.
        
        Parameters
        ----------
        trajectory_coords : np.ndarray
            Trajectory coordinates (n_frames, N, 3)
        align : bool
            Whether to perform Kabsch alignment
        
        Returns
        -------
        results : dict
            Analysis results
        """
        return analyze_trajectory(
            trajectory_coords, self.ref_coords, self.masses,
            self.eigenvectors, self.mode_indices, align
        )
    
    def get_mode_frequencies(self) -> np.ndarray:
        """Get frequencies for the selected modes."""
        if self.frequencies_cm is None:
            return None
        return self.frequencies_cm[self.mode_indices]
    
    def calculate_statistics(self, mode_coords: np.ndarray, 
                              temperature: float = 300.0) -> dict:
        """Calculate statistics for mode coordinates."""
        freq = self.get_mode_frequencies()
        if freq is None:
            freq = np.ones(self.n_modes)  # Dummy values
        return calculate_mode_statistics(mode_coords, freq, temperature)


def load_trajectory_mdtraj(dcd_file: str, top_file: str,
                            atom_indices: Optional[np.ndarray] = None,
                            unit: str = 'angstrom') -> Tuple[np.ndarray, int]:
    """
    Load trajectory using MDTraj.
    
    Parameters
    ----------
    dcd_file : str
        Path to DCD trajectory file
    top_file : str
        Path to topology file (PDB, PSF, etc.)
    atom_indices : np.ndarray, optional
        Indices of atoms to extract
    unit : str
        Output unit: 'angstrom' (default) or 'bohr'
    
    Returns
    -------
    coords : np.ndarray
        Coordinates (n_frames, N, 3)
    n_frames : int
        Number of frames
    """
    try:
        import mdtraj as md
    except ImportError:
        raise ImportError("MDTraj is required for loading DCD files. "
                         "Install with: pip install mdtraj")
    
    traj = md.load(dcd_file, top=top_file)
    
    if atom_indices is not None:
        traj = traj.atom_slice(atom_indices)
    
    # MDTraj uses nm, convert to Angstrom (1 nm = 10 Å)
    coords = traj.xyz * 10.0  # nm to Angstrom
    
    # Convert to Bohr if requested
    if unit.lower() == 'bohr':
        coords = coords * ANGSTROM_TO_BOHR
    
    return coords, traj.n_frames


def main():
    """Example usage of normal mode analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Normal mode coordinate analysis')
    parser.add_argument('--dcd', '-d', required=True, help='DCD trajectory file')
    parser.add_argument('--top', '-t', required=True, help='Topology file (PDB)')
    parser.add_argument('--gen', '-g', required=True, help='Reference structure (.gen)')
    parser.add_argument('--eigvec', '-e', required=True, help='Eigenvectors file')
    parser.add_argument('--eigval', '-v', help='Eigenvalues file')
    parser.add_argument('--output', '-o', default='mode_coords', help='Output prefix')
    parser.add_argument('--mode-start', type=int, default=6, help='First mode index')
    parser.add_argument('--mode-end', type=int, help='Last mode index')
    parser.add_argument('--no-align', action='store_true', help='Skip Kabsch alignment')
    
    args = parser.parse_args()
    
    # Load reference and eigenvectors
    print(f"Loading reference structure: {args.gen}")
    ref_coords, atom_types = read_gen_file(args.gen)
    
    print(f"Loading eigenvectors: {args.eigvec}")
    eigenvectors = read_hessian_eigenvectors(args.eigvec)
    
    eigenvalues, frequencies_cm = None, None
    if args.eigval:
        print(f"Loading eigenvalues: {args.eigval}")
        eigenvalues, frequencies_cm = read_hessian_eigenvalues(args.eigval)
    
    # Create analyzer
    n_dof = 3 * len(atom_types)
    mode_end = args.mode_end if args.mode_end else n_dof
    mode_range = (args.mode_start, mode_end)
    
    analyzer = NormalModeAnalyzer(
        ref_coords, atom_types, eigenvectors,
        eigenvalues, frequencies_cm, mode_range
    )
    
    print(f"\nReference structure: {analyzer.n_atoms} atoms")
    print(f"Mode range: {mode_range[0]} to {mode_range[1]} ({analyzer.n_modes} modes)")
    
    # Load trajectory
    print(f"\nLoading trajectory: {args.dcd}")
    traj_coords, n_frames = load_trajectory_mdtraj(args.dcd, args.top)
    print(f"  Frames: {n_frames}")
    print(f"  Atoms in trajectory: {traj_coords.shape[1]}")
    
    if traj_coords.shape[1] != analyzer.n_atoms:
        print(f"\nWarning: Atom count mismatch!")
        print(f"  Trajectory: {traj_coords.shape[1]}")
        print(f"  Reference: {analyzer.n_atoms}")
        return
    
    # Analyze
    print(f"\nPerforming analysis (align={not args.no_align})...")
    results = analyzer.analyze_trajectory(traj_coords, align=not args.no_align)
    
    # Statistics
    stats = analyzer.calculate_statistics(results['mode_coords'])
    
    print(f"\nRMSD: mean = {np.mean(results['rmsd']):.4f} Bohr")
    print(f"       std = {np.std(results['rmsd']):.4f} Bohr")
    
    # Save results
    output_file = f"{args.output}.npz"
    print(f"\nSaving results to: {output_file}")
    
    np.savez(output_file,
             mode_coords=results['mode_coords'],
             rmsd=results['rmsd'],
             mode_indices=results['mode_indices'],
             mode_mean=stats['mean'],
             mode_std=stats['std'],
             mode_variance=stats['variance'],
             frequencies_cm=analyzer.get_mode_frequencies())
    
    # Also save as text
    txt_file = f"{args.output}.txt"
    with open(txt_file, 'w') as f:
        f.write(f"# Normal mode coordinates\n")
        f.write(f"# Frames: {n_frames}, Modes: {analyzer.n_modes}\n")
        f.write(f"# Mode indices: {mode_range[0]} to {mode_range[1]}\n")
        f.write(f"# Columns: frame_index, mode_1, mode_2, ..., mode_n\n")
        for i in range(n_frames):
            line = f"{i:8d}"
            for d in results['mode_coords'][i]:
                line += f" {d:15.8e}"
            f.write(line + "\n")
    
    print(f"Saved text output: {txt_file}")
    print("\nDone!")


if __name__ == "__main__":
    main()
