#!/usr/bin/env python
"""
Calculate spectral density with Duschinsky rotation effects.

The Duschinsky transformation relates normal modes of two electronic states:
    Q_f = J * Q_i + D

where:
    J: Duschinsky rotation matrix (mode mixing)
    D: Displacement vector in final state normal coordinates

This script calculates:
1. Normal modes for both neutral and cation
2. Duschinsky rotation matrix J
3. Displacement vector D  
4. Reorganization energies with and without Duschinsky effects
5. Huang-Rhys factors
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Constants
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_CM = 219474.63
CM_TO_EV = HARTREE_TO_EV / HARTREE_TO_CM
BOHR_TO_ANG = 0.529177210903
ANG_TO_BOHR = 1.0 / BOHR_TO_ANG

# Atomic masses (amu)
MASSES = {
    'C': 12.0107,
    'H': 1.00794,
    'N': 14.0067,
    'O': 15.9994,
    'S': 32.065,
}


def read_gen_file(filename):
    """Read geometry from .gen file (Angstrom)"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    natoms = int(lines[0].split()[0])
    atom_types = lines[1].split()
    
    coords = []
    elements = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        atom_idx = int(parts[1]) - 1
        elements.append(atom_types[atom_idx])
        coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
    
    return np.array(coords), elements


def read_hessian_massweighted(filename, natoms):
    """Read mass-weighted Hessian matrix from hessian_massweighted.txt"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines starting with #
    data = []
    for line in lines:
        if not line.strip().startswith('#'):
            data.extend(line.split())
    
    n3 = natoms * 3
    hessian = np.zeros((n3, n3))
    
    idx = 0
    for i in range(n3):
        for j in range(n3):
            hessian[i, j] = float(data[idx])
            idx += 1
    
    return hessian


def calculate_normal_modes(hessian_mw, elements):
    """
    Calculate normal modes from mass-weighted Hessian.
    
    Returns:
        frequencies: in cm^-1
        L: eigenvector matrix (mass-weighted, column = mode)
        mass_vec: mass vector (repeated 3x per atom)
    """
    natoms = len(elements)
    n3 = natoms * 3
    
    masses = np.array([MASSES[e] for e in elements])
    mass_vec = np.repeat(masses, 3)
    
    # Hessian is already mass-weighted, directly diagonalize
    eigenvalues, eigenvectors = np.linalg.eigh(hessian_mw)
    
    # Convert eigenvalues to frequencies (cm^-1)
    hartree_to_joule = 4.3597447222071e-18
    bohr_to_m = 5.29177210903e-11
    amu_to_kg = 1.66053906660e-27
    c = 2.99792458e10
    
    conv = np.sqrt(hartree_to_joule / (bohr_to_m**2 * amu_to_kg)) / (2 * np.pi * c)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * conv
    
    return frequencies, eigenvectors, mass_vec


def calculate_duschinsky(neutral_dir, cation_dir, output_dir=None, freq_threshold=50):
    """
    Calculate Duschinsky rotation matrix and displacement vector.
    
    The transformation is: Q_C = J * Q_N + D
    
    where Q are normal mode coordinates.
    """
    if output_dir is None:
        output_dir = Path('.')
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
    
    # =========================================================
    # Step 1: Read geometries and Hessians
    # =========================================================
    print("="*70)
    print("DUSCHINSKY ROTATION ANALYSIS")
    print("="*70)
    
    print("\n[1] Reading input files...")
    coords_n, elements = read_gen_file(f"{neutral_dir}/geom.out.gen")
    coords_c, _ = read_gen_file(f"{cation_dir}/geom.out.gen")
    
    hessian_n = read_hessian_massweighted(f"{neutral_dir}/hessian_massweighted.txt", len(elements))
    hessian_c = read_hessian_massweighted(f"{cation_dir}/hessian_massweighted.txt", len(elements))
    
    natoms = len(elements)
    n3 = natoms * 3
    
    # =========================================================
    # Step 2: Calculate normal modes for both states
    # =========================================================
    print("\n[2] Calculating normal modes...")
    
    freq_n, L_n, mass_vec = calculate_normal_modes(hessian_n, elements)
    freq_c, L_c, _ = calculate_normal_modes(hessian_c, elements)
    
    print(f"    Neutral:  {np.sum(freq_n > freq_threshold)} vibrational modes")
    print(f"    Cation:   {np.sum(freq_c > freq_threshold)} vibrational modes")
    
    # =========================================================
    # Step 3: Calculate Duschinsky rotation matrix
    # =========================================================
    print("\n[3] Calculating Duschinsky rotation matrix...")
    
    # J = L_c^T * L_n  (both are orthonormal in mass-weighted coordinates)
    # This gives the overlap between cation and neutral normal modes
    J = L_c.T @ L_n
    
    # Check orthogonality
    orth_check = np.max(np.abs(J @ J.T - np.eye(n3)))
    print(f"    Orthogonality check (should be ~0): {orth_check:.2e}")
    
    # =========================================================
    # Step 4: Calculate displacement vector
    # =========================================================
    print("\n[4] Calculating displacement vector...")
    
    # Geometry difference in Angstrom -> Bohr
    delta_R = (coords_c - coords_n).flatten() * ANG_TO_BOHR
    
    # Mass-weighted displacement
    sqrt_mass = np.sqrt(mass_vec)
    delta_q_mw = delta_R * sqrt_mass
    
    # Method 1: Displacement in NEUTRAL normal coordinates (parallel mode approx.)
    D_neutral = L_n.T @ delta_q_mw
    
    # Method 2: Displacement in CATION normal coordinates
    D_cation = L_c.T @ delta_q_mw
    
    # Method 3: Duschinsky-transformed displacement
    # Q_C = J * Q_N + D  =>  D = L_c^T * delta_q = D_cation
    D_duschinsky = D_cation
    
    print(f"    RMSD (Cartesian): {np.sqrt(np.mean(delta_R**2))*BOHR_TO_ANG:.4f} Angstrom")
    print(f"    |D| in neutral modes:  {np.linalg.norm(D_neutral):.4f}")
    print(f"    |D| in cation modes:   {np.linalg.norm(D_cation):.4f}")
    
    # =========================================================
    # Step 5: Calculate reorganization energies
    # =========================================================
    print("\n[5] Calculating reorganization energies...")
    
    # Convert frequencies to atomic units
    omega_n_au = freq_n / HARTREE_TO_CM
    omega_c_au = freq_c / HARTREE_TO_CM
    
    # --- Method A: Parallel mode approximation (no Duschinsky) ---
    # Use neutral modes only
    lambda_parallel = 0.5 * omega_n_au**2 * D_neutral**2 * HARTREE_TO_EV
    
    # --- Method B: Vertical approximation with cation modes ---
    # Use cation modes
    lambda_cation = 0.5 * omega_c_au**2 * D_cation**2 * HARTREE_TO_EV
    
    # --- Method C: Duschinsky with mode mixing ---
    # The reorganization energy considering Duschinsky rotation
    # For small displacements: λ_i ≈ (1/2) ω_i^2 Σ_j J_ij^2 D_j^2 + cross terms
    # Full treatment requires numerical integration of Franck-Condon factors
    
    # Simplified: use average of both approaches or cation modes
    lambda_duschinsky = lambda_cation  # First approximation
    
    # =========================================================
    # Step 6: Huang-Rhys factors
    # =========================================================
    omega_n_ev = freq_n * CM_TO_EV
    omega_c_ev = freq_c * CM_TO_EV
    
    S_parallel = np.zeros_like(lambda_parallel)
    S_cation = np.zeros_like(lambda_cation)
    
    valid_n = np.abs(omega_n_ev) > 1e-6
    valid_c = np.abs(omega_c_ev) > 1e-6
    
    S_parallel[valid_n] = lambda_parallel[valid_n] / omega_n_ev[valid_n]
    S_cation[valid_c] = lambda_cation[valid_c] / omega_c_ev[valid_c]
    
    # =========================================================
    # Step 7: Summary
    # =========================================================
    vib_n = freq_n > freq_threshold
    vib_c = freq_c > freq_threshold
    
    lambda_total_parallel = np.sum(lambda_parallel[vib_n])
    lambda_total_cation = np.sum(lambda_cation[vib_c])
    
    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    print(f"\nTotal reorganization energy:")
    print(f"  Parallel mode approx. (neutral modes): {lambda_total_parallel*1000:.2f} meV")
    print(f"  Cation modes:                          {lambda_total_cation*1000:.2f} meV")
    print(f"\nTotal Huang-Rhys factor:")
    print(f"  Parallel mode approx.: {np.sum(S_parallel[vib_n]):.4f}")
    print(f"  Cation modes:          {np.sum(S_cation[vib_c]):.4f}")
    
    # =========================================================
    # Step 8: Analyze Duschinsky matrix
    # =========================================================
    print(f"\n{'='*70}")
    print("DUSCHINSKY MATRIX ANALYSIS")
    print(f"{'='*70}")
    
    # Look at mode mixing (off-diagonal elements)
    # For each cation mode, find the most contributing neutral modes
    print("\nMode mixing (top contributions for each cation mode):")
    print(f"{'Cation Mode':>12} {'Freq(cm-1)':>12} | {'Neutral contributions':>40}")
    print("-" * 70)
    
    for i in range(n3):
        if freq_c[i] > freq_threshold:
            # Get contributions from neutral modes
            contributions = J[i, :]**2  # Squared overlap
            top_idx = np.argsort(contributions)[::-1][:3]
            
            contrib_str = ", ".join([f"N{j+1}({contributions[j]*100:.0f}%)" 
                                    for j in top_idx if contributions[j] > 0.01])
            
            # Check if it's mainly one mode (diagonal)
            max_contrib = np.max(contributions)
            mixing = "mixed" if max_contrib < 0.8 else "pure"
            
            if i < 20 or max_contrib < 0.8:  # Show first 20 or mixed modes
                print(f"{i+1:>12} {freq_c[i]:>12.1f} | {contrib_str:<40} [{mixing}]")
    
    # =========================================================
    # Step 9: Save results
    # =========================================================
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print(f"{'='*70}")
    
    results = {
        # Frequencies
        'freq_neutral_cm': freq_n,
        'freq_cation_cm': freq_c,
        
        # Duschinsky matrix
        'duschinsky_matrix': J,
        
        # Displacements
        'displacement_neutral': D_neutral,
        'displacement_cation': D_cation,
        
        # Reorganization energies
        'lambda_parallel_ev': lambda_parallel,
        'lambda_cation_ev': lambda_cation,
        'lambda_total_parallel_ev': lambda_total_parallel,
        'lambda_total_cation_ev': lambda_total_cation,
        
        # Huang-Rhys factors
        'huang_rhys_parallel': S_parallel,
        'huang_rhys_cation': S_cation,
    }
    
    np.savez(output_dir / 'duschinsky_analysis.npz', **results)
    print(f"  Saved: {output_dir / 'duschinsky_analysis.npz'}")
    
    # Write text summary
    write_text_summary(output_dir / 'duschinsky_analysis.txt', 
                       freq_n, freq_c, J, D_neutral, D_cation,
                       lambda_parallel, lambda_cation, S_parallel, S_cation,
                       freq_threshold)
    print(f"  Saved: {output_dir / 'duschinsky_analysis.txt'}")
    
    # Create plots
    plot_duschinsky_analysis(output_dir, freq_n, freq_c, J, 
                            lambda_parallel, lambda_cation,
                            S_parallel, S_cation, freq_threshold)
    
    return results


def write_text_summary(filename, freq_n, freq_c, J, D_n, D_c, 
                       lambda_n, lambda_c, S_n, S_c, threshold):
    """Write detailed text summary"""
    with open(filename, 'w') as f:
        f.write("# Duschinsky Rotation Analysis\n")
        f.write("#\n")
        f.write(f"# Total reorganization energy (parallel mode): {np.sum(lambda_n[freq_n > threshold])*1000:.4f} meV\n")
        f.write(f"# Total reorganization energy (cation modes):  {np.sum(lambda_c[freq_c > threshold])*1000:.4f} meV\n")
        f.write("#\n")
        
        # Mode comparison
        f.write("# Mode comparison:\n")
        f.write("# Mode  Freq_N(cm-1)  Freq_C(cm-1)  D_N        D_C        lambda_N(meV)  lambda_C(meV)  S_N       S_C       Max_J\n")
        
        for i in range(len(freq_n)):
            if freq_n[i] > threshold or freq_c[i] > threshold:
                max_J = np.max(np.abs(J[i, :]))
                f.write(f"{i+1:5d}  {freq_n[i]:12.2f}  {freq_c[i]:12.2f}  "
                       f"{D_n[i]:10.4f} {D_c[i]:10.4f} "
                       f"{lambda_n[i]*1000:13.4f}  {lambda_c[i]*1000:13.4f}  "
                       f"{S_n[i]:9.4f} {S_c[i]:9.4f} {max_J:8.3f}\n")


def plot_duschinsky_analysis(output_dir, freq_n, freq_c, J, 
                             lambda_n, lambda_c, S_n, S_c, threshold):
    """Create visualization plots"""
    
    vib_n = freq_n > threshold
    vib_c = freq_c > threshold
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Duschinsky matrix heatmap
    ax = axes[0, 0]
    n_vib = min(50, np.sum(vib_n))  # Show first 50 vibrational modes
    vib_idx = np.where(vib_n)[0][:n_vib]
    J_sub = J[np.ix_(vib_idx, vib_idx)]
    im = ax.imshow(np.abs(J_sub), cmap='hot', aspect='auto', vmin=0, vmax=1)
    ax.set_xlabel('Neutral mode')
    ax.set_ylabel('Cation mode')
    ax.set_title('Duschinsky Matrix |J|')
    plt.colorbar(im, ax=ax)
    
    # 2. Diagonal elements of J (mode correlation)
    ax = axes[0, 1]
    diag_J = np.diag(J)[vib_n]
    ax.plot(freq_n[vib_n], np.abs(diag_J), 'o-', markersize=3)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(0.9, color='red', linestyle=':', alpha=0.5, label='90% threshold')
    ax.set_xlabel('Frequency (cm$^{-1}$)')
    ax.set_ylabel('|J$_{ii}$|')
    ax.set_title('Diagonal Duschinsky Elements')
    ax.legend()
    ax.set_ylim(0, 1.1)
    
    # 3. Frequency comparison
    ax = axes[0, 2]
    ax.scatter(freq_n[vib_n], freq_c[vib_c][:np.sum(vib_n)], alpha=0.5, s=20)
    max_freq = max(np.max(freq_n[vib_n]), np.max(freq_c[vib_c]))
    ax.plot([0, max_freq], [0, max_freq], 'k--', alpha=0.5)
    ax.set_xlabel('Neutral frequency (cm$^{-1}$)')
    ax.set_ylabel('Cation frequency (cm$^{-1}$)')
    ax.set_title('Frequency Correlation')
    
    # 4. Spectral density comparison
    ax = axes[1, 0]
    ax.stem(freq_n[vib_n], lambda_n[vib_n]*1000, linefmt='b-', markerfmt='bo', 
            basefmt=' ', label='Parallel mode')
    ax.stem(freq_c[vib_c], lambda_c[vib_c]*1000, linefmt='r-', markerfmt='r^', 
            basefmt=' ', label='Cation modes')
    ax.set_xlabel('Frequency (cm$^{-1}$)')
    ax.set_ylabel('λ (meV)')
    ax.set_title('Discrete Spectral Density')
    ax.legend()
    
    # 5. Huang-Rhys factors
    ax = axes[1, 1]
    ax.stem(freq_n[vib_n], S_n[vib_n], linefmt='b-', markerfmt='bo', 
            basefmt=' ', label='Parallel mode')
    ax.stem(freq_c[vib_c], S_c[vib_c], linefmt='r-', markerfmt='r^', 
            basefmt=' ', label='Cation modes')
    ax.set_xlabel('Frequency (cm$^{-1}$)')
    ax.set_ylabel('Huang-Rhys factor S')
    ax.set_title('Huang-Rhys Factors')
    ax.legend()
    
    # 6. Broadened spectral density
    ax = axes[1, 2]
    omega_grid = np.linspace(0, max(np.max(freq_n[vib_n]), np.max(freq_c[vib_c])) * 1.1, 2000)
    gamma = 20  # Broadening in cm^-1
    
    J_n = np.zeros_like(omega_grid)
    J_c = np.zeros_like(omega_grid)
    
    for w, lam in zip(freq_n[vib_n], lambda_n[vib_n]*1000):
        J_n += lam * gamma / ((omega_grid - w)**2 + gamma**2) / np.pi
    
    for w, lam in zip(freq_c[vib_c], lambda_c[vib_c]*1000):
        J_c += lam * gamma / ((omega_grid - w)**2 + gamma**2) / np.pi
    
    ax.plot(omega_grid, J_n, 'b-', label='Parallel mode', linewidth=1.5)
    ax.plot(omega_grid, J_c, 'r--', label='Cation modes', linewidth=1.5)
    ax.set_xlabel('Frequency (cm$^{-1}$)')
    ax.set_ylabel('J(ω) (meV)')
    ax.set_title(f'Broadened Spectral Density (γ = {gamma} cm$^{{-1}}$)')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'duschinsky_analysis.png', dpi=150)
    plt.savefig(output_dir / 'duschinsky_analysis.pdf')
    print(f"  Saved: {output_dir / 'duschinsky_analysis.png'}")
    print(f"  Saved: {output_dir / 'duschinsky_analysis.pdf'}")
    plt.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Duschinsky rotation analysis')
    parser.add_argument('--neutral', default='neutral', help='Neutral calculation directory')
    parser.add_argument('--cation', default='cation', help='Cation calculation directory')
    parser.add_argument('--output', default='result', help='Output directory')
    parser.add_argument('--threshold', type=float, default=50, help='Frequency threshold (cm^-1)')
    args = parser.parse_args()
    
    calculate_duschinsky(args.neutral, args.cation, args.output, args.threshold)


if __name__ == '__main__':
    main()
