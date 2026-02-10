#!/usr/bin/env python
"""
DFTB+ Hessian matrix loading, symmetry check, and diagonalization script

Output format:
- 4 values per line
- Each value is 16 characters wide with 10 decimal places (f16.10)
- Hessian matrix output column by column

Note:
- DFTB+ hessian.out is a regular Hessian (not mass-weighted)
- To calculate frequencies, conversion to mass-weighted Hessian is required
"""

import numpy as np
import os

# Atomic masses (element -> mass [amu])
# Using values from DFTB+ SK files (3ob)
ATOMIC_MASSES = {
    'H': 1.008,       # SKファイル: 1.008
    'He': 4.002602,
    'Li': 6.941,
    'Be': 9.012182,
    'B': 10.811,
    'C': 12.01,       # SKファイル: 12.01
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


def read_gen_file(filename):
    """
    Read atom types from DFTB+ .gen file
    
    Parameters
    ----------
    filename : str
        Path to .gen file
    
    Returns
    -------
    atom_types : list
        List of element symbols for each atom
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Line 1: Number of atoms and coordinate type
    first_line = lines[0].split()
    n_atoms = int(first_line[0])
    
    # Line 2: List of element types
    element_types = lines[1].split()
    
    # Line 3 onwards: Atom data
    atom_types = []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        type_idx = int(parts[1]) - 1  # 1-indexed -> 0-indexed
        atom_types.append(element_types[type_idx])
    
    return atom_types


def read_hessian(filename):
    """
    Read DFTB+ hessian.out file
    
    Parameters
    ----------
    filename : str
        Path to hessian.out file
    
    Returns
    -------
    hessian : ndarray
        Hessian matrix (N x N)
    """
    # Read all values from file
    values = []
    with open(filename, 'r') as f:
        for line in f:
            # Fixed 16-character width format, split by whitespace
            parts = line.split()
            for val in parts:
                values.append(float(val))
    
    values = np.array(values)
    
    # Estimate matrix size (N x N square matrix)
    n = int(np.sqrt(len(values)))
    if n * n != len(values):
        raise ValueError(f"Number of values {len(values)} does not form a square matrix")
    
    # Output is column-by-column, so reshape in Fortran order (column-major)
    hessian = values.reshape((n, n), order='F')
    
    return hessian


def check_symmetry(hessian, tol=1e-8):
    """
    Check symmetry of Hessian matrix
    
    Parameters
    ----------
    hessian : ndarray
        Hessian matrix
    tol : float
        Tolerance for symmetry
    
    Returns
    -------
    is_symmetric : bool
        Whether the matrix is symmetric
    max_diff : float
        Maximum difference
    mean_diff : float
        Mean difference
    """
    diff = hessian - hessian.T
    max_diff = np.max(np.abs(diff))
    mean_diff = np.mean(np.abs(diff))
    
    is_symmetric = max_diff < tol
    
    return is_symmetric, max_diff, mean_diff


def symmetrize(hessian):
    """Symmetrize Hessian matrix"""
    return (hessian + hessian.T) / 2


def mass_weight_hessian(hessian, atom_types, use_amu=True):
    """
    Convert Hessian to mass-weighted Hessian
    
    H_mw[i,j] = H[i,j] / sqrt(m_i * m_j)
    
    Parameters
    ----------
    hessian : ndarray
        Regular Hessian matrix (3N x 3N)
    atom_types : list
        List of element symbols for each atom (length N)
    use_amu : bool
        True: Calculate using amu like DFTB+ (DFTB+ compatible mode)
        False: Convert to atomic units (me) for calculation
    
    Returns
    -------
    hessian_mw : ndarray
        Mass-weighted Hessian matrix
    masses_3n : ndarray
        Mass for each degree of freedom (amu or atomic units)
    """
    n_atoms = len(atom_types)
    n = 3 * n_atoms
    
    # Mass of each atom (amu) - using values from DFTB+ SK files
    # Note: ATOMIC_MASSES uses standard atomic weights, may differ slightly from SK files
    masses = np.array([ATOMIC_MASSES[atom] for atom in atom_types])
    
    # Expand to 3N dimensions (same mass for x, y, z)
    masses_3n = np.repeat(masses, 3)
    
    if use_amu:
        # DFTB+ modes compatible mode: calculate using amu
        # Eigenvalue units are Hartree / (Bohr^2 * amu)
        sqrt_masses = np.sqrt(masses_3n)
    else:
        # Atomic units mode: convert amu -> electron mass
        # 1 amu = 1822.888486 m_e
        amu_to_me = 1822.888486209
        masses_3n_au = masses_3n * amu_to_me
        sqrt_masses = np.sqrt(masses_3n_au)
    
    # Mass-weighted Hessian: H_mw = M^(-1/2) H M^(-1/2)
    mass_matrix_inv_sqrt = np.diag(1.0 / sqrt_masses)
    hessian_mw = mass_matrix_inv_sqrt @ hessian @ mass_matrix_inv_sqrt
    
    return hessian_mw, masses_3n


def diagonalize(hessian):
    """
    Diagonalize Hessian matrix
    
    Parameters
    ----------
    hessian : ndarray
        Hessian matrix
    
    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues (ascending order)
    eigenvectors : ndarray
        Eigenvectors
    """
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    return eigenvalues, eigenvectors


def convert_to_frequencies(eigenvalues, use_amu=True):
    """
    Convert eigenvalues to frequencies
    
    DFTB+ modes compatible mode (use_amu=True):
        Calculate frequency from eigenvalues when mass is in amu
        ω = sqrt(λ) [units: sqrt(Ha/(Bohr^2*amu))]
        frequency(cm^-1) = ω * Hartree_cm / sqrt(amu_to_me)
        
    Atomic units mode (use_amu=False):
        Calculate frequency from eigenvalues when mass is in atomic units (me)
        ω = sqrt(λ) [atomic units: Ha]
        frequency(cm^-1) = ω * Hartree_cm
    
    Parameters
    ----------
    eigenvalues : ndarray
        Eigenvalues of mass-weighted Hessian
    use_amu : bool
        True: DFTB+ modes compatible mode
        False: Atomic units mode
    
    Returns
    -------
    frequencies_cm : ndarray
        Frequencies (cm^-1)
    frequencies_hartree : ndarray
        Frequencies (Hartree units, for comparison with DFTB+ output)
    """
    # Conversion constants
    hartree_to_cm = 219474.63137  # 1 Hartree = 219474.63 cm^-1
    amu_to_me = 1822.888486209    # 1 amu = 1822.89 * m_e
    
    frequencies_hartree = np.zeros_like(eigenvalues)
    frequencies_cm = np.zeros_like(eigenvalues)
    
    for i, ev in enumerate(eigenvalues):
        if ev >= 0:
            omega = np.sqrt(ev)
        else:
            omega = -np.sqrt(-ev)  # Negative eigenvalues displayed as negative (imaginary frequency)
        
        if use_amu:
            # DFTB+ compatible: ω units are sqrt(Ha/(Bohr^2*amu))
            # Convert to atomic units: ω_au = ω / sqrt(amu_to_me)
            omega_au = omega / np.sqrt(amu_to_me)
            freq_hartree = omega_au  # DFTB+ outputs this as "frequencies"
            freq_cm = omega_au * hartree_to_cm
        else:
            # Atomic units mode: ω is already in atomic units (Ha)
            freq_hartree = omega
            freq_cm = omega * hartree_to_cm
        
        frequencies_hartree[i] = freq_hartree
        frequencies_cm[i] = freq_cm
    
    return frequencies_cm, frequencies_hartree


def save_results(filename, hessian_orig, hessian_sym, hessian_mw, eigenvalues, eigenvectors, 
                 frequencies_cm, frequencies_hartree, atom_types):
    """
    Save analysis results to files
    
    Parameters
    ----------
    filename : str
        Output filename (without extension)
    hessian_orig : ndarray
        Original Hessian matrix
    hessian_sym : ndarray
        Symmetrized Hessian matrix
    hessian_mw : ndarray
        Mass-weighted Hessian matrix
    eigenvalues : ndarray
        Eigenvalues
    eigenvectors : ndarray
        Eigenvectors
    frequencies_cm : ndarray
        Frequencies (cm^-1)
    frequencies_hartree : ndarray
        Frequencies (Hartree units)
    atom_types : list
        List of atom types
    """
    n = hessian_orig.shape[0]
    n_atoms = n // 3
    
    # 1. Analysis report
    report_file = f"{filename}_analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("DFTB+ Hessian Analysis Report\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("## Matrix Information\n")
        f.write(f"  Size: {n} x {n}\n")
        f.write(f"  Number of atoms: {n_atoms}\n")
        f.write(f"  Degrees of freedom: {n}\n\n")
        
        # Symmetry
        is_sym, max_diff, mean_diff = check_symmetry(hessian_orig)
        f.write("## Symmetry Check (Original Hessian)\n")
        f.write(f"  Max |H - H^T|: {max_diff:.6e}\n")
        f.write(f"  Mean |H - H^T|: {mean_diff:.6e}\n")
        f.write(f"  Symmetry: {'symmetric' if is_sym else 'asymmetric (symmetrization applied)'}\n\n")
        
        # Eigenvalue analysis
        n_zero_1e4 = np.sum(np.abs(eigenvalues) < 1e-4)
        n_zero_1e5 = np.sum(np.abs(eigenvalues) < 1e-5)
        n_negative = np.sum(eigenvalues < -1e-4)
        n_positive = np.sum(eigenvalues > 1e-4)
        
        f.write("## Eigenvalue Analysis\n")
        f.write(f"  Total eigenvalues: {len(eigenvalues)}\n")
        f.write(f"  |lambda| < 1e-4: {n_zero_1e4} (rotational/translational modes)\n")
        f.write(f"  |lambda| < 1e-5: {n_zero_1e5}\n")
        f.write(f"  lambda < -1e-4: {n_negative}\n")
        f.write(f"  lambda > 1e-4: {n_positive} (vibrational modes)\n\n")
        
        f.write("## Rotational/Translational Modes (6 smallest eigenvalues)\n")
        for i in range(min(6, len(eigenvalues))):
            f.write(f"  lambda_{i+1:3d} = {eigenvalues[i]:+15.10e}\n")
        f.write("\n")
        
        f.write("## Conclusion\n")
        if n_zero_1e4 == 6:
            f.write("  OK: 6 small eigenvalues (rotational/translational modes) correctly detected\n")
            f.write("  OK: Normal result for nonlinear molecule\n")
        elif n_zero_1e4 == 5:
            f.write("  OK: 5 small eigenvalues detected (possibly linear molecule)\n")
        elif n_zero_1e4 == 3:
            f.write("  Note: 3 small eigenvalues (translation only, rotation may be constrained)\n")
        else:
            f.write(f"  Note: {n_zero_1e4} small eigenvalues detected\n")
        
        f.write("\n" + "=" * 70 + "\n")
    
    print(f"  -> {report_file}")
    
    # 2. All eigenvalues
    eigenvalues_file = f"{filename}_eigenvalues.txt"
    with open(eigenvalues_file, 'w') as f:
        f.write("# DFTB+ Mass-weighted Hessian eigenvalues and frequencies\n")
        f.write(f"# Matrix size: {n} x {n}, Number of atoms: {n_atoms}\n")
        f.write(f"# Atom types: {' '.join(atom_types)}\n")
        f.write("# Index    Eigenvalue           Freq(Hartree)      Freq(cm^-1)    Type\n")
        f.write("#" + "-" * 80 + "\n")
        for i, (ev, freq_h, freq_cm) in enumerate(zip(eigenvalues, frequencies_hartree, frequencies_cm)):
            if abs(ev) < 1e-7:
                mode_type = "rot/trans"
            elif ev < 0:
                mode_type = "imaginary"
            else:
                mode_type = "vibration"
            f.write(f"{i+1:6d}  {ev:+20.12e}  {freq_h:+15.10e}  {freq_cm:+12.4f}     {mode_type}\n")
    
    print(f"  -> {eigenvalues_file}")
    
    # 3. Hessian matrix (before symmetrization)
    hessian_orig_file = f"{filename}_original.txt"
    with open(hessian_orig_file, 'w') as f:
        f.write(f"# DFTB+ Hessian matrix (original data)\n")
        f.write(f"# Size: {n} x {n}\n")
        f.write(f"# Number of atoms: {n_atoms}\n\n")
        for i in range(n):
            row_str = " ".join([f"{hessian_orig[i, j]:+15.10e}" for j in range(n)])
            f.write(row_str + "\n")
    
    print(f"  -> {hessian_orig_file}")
    
    # 4. Hessian matrix (after symmetrization)
    hessian_sym_file = f"{filename}_symmetrized.txt"
    with open(hessian_sym_file, 'w') as f:
        f.write(f"# DFTB+ Hessian matrix (symmetrized)\n")
        f.write(f"# Size: {n} x {n}\n")
        f.write(f"# Number of atoms: {n_atoms}\n")
        f.write(f"# Symmetrization: (H + H^T) / 2\n\n")
        for i in range(n):
            row_str = " ".join([f"{hessian_sym[i, j]:+15.10e}" for j in range(n)])
            f.write(row_str + "\n")
    
    print(f"  -> {hessian_sym_file}")
    
    # 4b. Mass-weighted Hessian matrix
    hessian_mw_file = f"{filename}_massweighted.txt"
    with open(hessian_mw_file, 'w') as f:
        f.write(f"# DFTB+ Mass-weighted Hessian matrix\n")
        f.write(f"# Size: {n} x {n}\n")
        f.write(f"# Number of atoms: {n_atoms}\n")
        f.write(f"# H_mw = M^(-1/2) H M^(-1/2)\n")
        f.write(f"# Atom types: {' '.join(atom_types)}\n\n")
        for i in range(n):
            row_str = " ".join([f"{hessian_mw[i, j]:+15.10e}" for j in range(n)])
            f.write(row_str + "\n")
    
    print(f"  -> {hessian_mw_file}")
    
    # 5. Eigenvectors
    eigenvectors_file = f"{filename}_eigenvectors.txt"
    with open(eigenvectors_file, 'w') as f:
        f.write(f"# DFTB+ Mass-weighted Hessian eigenvectors\n")
        f.write(f"# Size: {n} x {n}\n")
        f.write(f"# Each column is one eigenvector (ascending eigenvalue order)\n\n")
        for i in range(n):
            row_str = " ".join([f"{eigenvectors[i, j]:+15.10e}" for j in range(n)])
            f.write(row_str + "\n")
    
    print(f"  -> {eigenvectors_file}")
    
    # 6. Save in NumPy format (easy to load later)
    npz_file = f"{filename}_data.npz"
    np.savez(npz_file,
             hessian_original=hessian_orig,
             hessian_symmetrized=hessian_sym,
             hessian_mass_weighted=hessian_mw,
             eigenvalues=eigenvalues,
             eigenvectors=eigenvectors,
             frequencies_cm=frequencies_cm,
             frequencies_hartree=frequencies_hartree,
             atom_types=atom_types)
    print(f"  -> {npz_file}")


def main():
    import sys
    
    # Filename
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "hessian.out"
    
    # .gen file (atom type information)
    if len(sys.argv) > 2:
        gen_file = sys.argv[2]
    else:
        # Default: look for geom.out.gen in the same directory
        gen_file = os.path.join(os.path.dirname(filename) or '.', 'opt.gen')
        if not os.path.exists(gen_file):
            gen_file = os.path.join(os.path.dirname(filename) or '.', 'geo.gen')
    
    # Output file basename
    basename = os.path.splitext(filename)[0]
    
    print(f"Input file: {filename}")
    print(f"Structure file: {gen_file}")
    print("=" * 60)
    
    # Read atom types
    if os.path.exists(gen_file):
        atom_types = read_gen_file(gen_file)
        print(f"Atom types: {set(atom_types)}")
        print(f"Number of atoms: {len(atom_types)} ({', '.join([f'{a}:{atom_types.count(a)}' for a in set(atom_types)])})")
    else:
        print(f"Warning: {gen_file} not found")
        print("  Assuming pentacene (C22H14)")
        atom_types = ['C'] * 22 + ['H'] * 14
    
    # Read Hessian
    hessian_orig = read_hessian(filename)
    n = hessian_orig.shape[0]
    n_atoms = n // 3
    
    print(f"\nHessian matrix size: {n} x {n}")
    
    if len(atom_types) != n_atoms:
        print(f"Warning: Number of atoms mismatch (Hessian: {n_atoms}, gen: {len(atom_types)})")
    
    # Symmetry check
    is_symmetric, max_diff, mean_diff = check_symmetry(hessian_orig)
    
    print(f"\n=== Symmetry Check ===")
    print(f"Max H - H^T: {max_diff:.2e}")
    print(f"Mean H - H^T: {mean_diff:.2e}")
    
    if is_symmetric:
        print(f"OK: Matrix is symmetric")
        hessian_sym = hessian_orig.copy()
    else:
        print(f"NG: Matrix is asymmetric")
        print("  -> Applying symmetrization (H + H^T) / 2")
        hessian_sym = symmetrize(hessian_orig)
    
    # Convert to mass-weighted Hessian (DFTB+ modes compatible mode: calculate using amu)
    print(f"\n=== Mass-weighted Hessian (DFTB+ modes compatible) ===")
    hessian_mw, masses = mass_weight_hessian(hessian_sym, atom_types, use_amu=True)
    print(f"Calculating H_mw = M^(-1/2) H M^(-1/2) (mass in amu)")
    
    # Diagonalize mass-weighted Hessian
    eigenvalues, eigenvectors = diagonalize(hessian_mw)
    
    # Eigenvalue analysis (screen output)
    print(f"\n=== Eigenvalue Analysis (Mass-weighted Hessian) ===")
    print(f"Number of eigenvalues: {len(eigenvalues)}")
    
    n_zero = np.sum(np.abs(eigenvalues) < 1e-8)
    print(f"|lambda| < 1e-8: {n_zero} (rotational/translational modes)")
    
    print(f"\n10 smallest eigenvalues:")
    for i in range(min(10, len(eigenvalues))):
        ev = eigenvalues[i]
        label = ""
        if abs(ev) < 1e-8:
            label = " <- rotational/translational mode"
        elif ev < 0:
            label = " <- negative eigenvalue"
        print(f"  {i+1:3d}: {ev:+15.10e}{label}")
    
    # Frequency conversion (DFTB+ modes compatible mode)
    frequencies_cm, frequencies_hartree = convert_to_frequencies(eigenvalues, use_amu=True)
    
    print(f"\n12 smallest frequencies:")
    print(f"{'Index':>6} {'Hartree':>15} {'cm^-1':>12}  Note")
    print("-" * 50)
    for i in range(min(12, len(frequencies_cm))):
        freq_h = frequencies_hartree[i]
        freq_cm = frequencies_cm[i]
        label = ""
        if abs(freq_cm) < 1:
            label = "rot/trans"
        elif freq_cm < 0:
            label = "imaginary"
        print(f"  {i+1:3d}: {freq_h:+15.10e} {freq_cm:+12.2f}  {label}")
    
    # Save results to files
    print(f"\n=== File Output ===")
    save_results(basename, hessian_orig, hessian_sym, hessian_mw, eigenvalues, eigenvectors, 
                 frequencies_cm, frequencies_hartree, atom_types)
    
    print("\n" + "=" * 60)
    print("Analysis complete")


if __name__ == "__main__":
    main()
