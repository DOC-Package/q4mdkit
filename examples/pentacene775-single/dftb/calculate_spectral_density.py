#!/usr/bin/env python
"""
Calculate spectral density from neutral and cation normal modes.

The spectral density J(ω) characterizes the electron-phonon coupling:
  J(ω) = Σ_i λ_i δ(ω - ω_i)

where λ_i is the reorganization energy for mode i:
  λ_i = (1/2) ω_i (ΔQ_i)^2

and ΔQ_i is the displacement along normal mode i between neutral and cation geometries.

The Huang-Rhys factor is:
  S_i = λ_i / (ℏω_i)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path

# Constants
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_CM = 219474.63  # cm^-1
CM_TO_EV = HARTREE_TO_EV / HARTREE_TO_CM
BOHR_TO_ANG = 0.529177210903

# Atomic masses (amu)
MASSES = {
    'C': 12.0107,
    'H': 1.00794,
}


def read_gen_file(filename):
    """Read geometry from .gen file"""
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
    """Calculate normal modes from mass-weighted Hessian
    
    DFTB+ outputs Hessian in Hartree/Bohr^2, and mass-weighted Hessian
    uses atomic mass units (amu). The eigenvalues have units of Hartree/(Bohr^2 * amu).
    
    Returns:
        frequencies: Vibrational frequencies in cm^-1
        eigenvectors: Normal mode eigenvectors
        mass_vec: Mass vector (amu, repeated 3x per atom)
        eigenvalues: Raw eigenvalues from diagonalization [Hartree/(Bohr^2 * amu)]
    """
    natoms = len(elements)
    n3 = natoms * 3
    
    masses = np.array([MASSES[e] for e in elements])
    mass_vec = np.repeat(masses, 3)
    
    # Hessian is already mass-weighted, directly diagonalize
    eigenvalues, eigenvectors = np.linalg.eigh(hessian_mw)
    
    # Convert eigenvalues [Hartree/(Bohr^2 * amu)] to frequencies [cm^-1]
    # ω = sqrt(k/m), where k is in Hartree/Bohr^2 and m is in amu
    # For mass-weighted eigenvalues: ω^2 = eigenvalue [Hartree/(Bohr^2 * amu)]
    # 
    # Unit conversion:
    # 1 Hartree = 4.3597447e-18 J
    # 1 Bohr = 5.29177e-11 m
    # 1 amu = 1.66054e-27 kg
    # c = 2.99792458e10 cm/s
    #
    # sqrt(Hartree/(Bohr^2 * amu)) -> rad/s -> cm^-1
    hartree_to_joule = 4.3597447222071e-18
    bohr_to_m = 5.29177210903e-11
    amu_to_kg = 1.66053906660e-27
    c = 2.99792458e10  # speed of light in cm/s
    
    # Conversion factor from sqrt(Hartree/(Bohr^2 * amu)) to cm^-1
    conv = np.sqrt(hartree_to_joule / (bohr_to_m**2 * amu_to_kg)) / (2 * np.pi * c)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * conv
    
    return frequencies, eigenvectors, mass_vec, eigenvalues


def calculate_spectral_density(neutral_dir, cation_dir, output_dir=None):
    """
    Calculate spectral density from neutral and cation calculations.
    
    Uses neutral normal modes as the basis (common approximation).
    """
    if output_dir is None:
        output_dir = Path('.')
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
    
    # Read neutral geometry and Hessian
    print("Reading neutral calculation...")
    coords_n, elements = read_gen_file(f"{neutral_dir}/geom.out.gen")
    hessian_n = read_hessian_massweighted(f"{neutral_dir}/hessian_massweighted.txt", len(elements))
    
    # Read cation geometry
    print("Reading cation calculation...")
    coords_c, _ = read_gen_file(f"{cation_dir}/geom.out.gen")
    
    # Calculate neutral normal modes
    print("Calculating normal modes...")
    frequencies, eigenvectors, mass_vec = calculate_normal_modes(hessian_n, elements)
    
    natoms = len(elements)
    n_modes = len(frequencies)
    
    # Calculate geometry difference
    # DFTB+ .gen files output coordinates in Bohr (atomic units)
    # No conversion needed as Hessian is also in atomic units (Hartree/Bohr^2)
    delta_coords = (coords_c - coords_n).flatten()  # Already in Bohr
    
    print(f"  Geometry difference: RMSD = {np.sqrt(np.mean(delta_coords**2))*BOHR_TO_ANG:.4f} Angstrom")
    
    # Mass-weight the displacement
    sqrt_mass = np.sqrt(mass_vec)
    delta_q_mw = delta_coords * sqrt_mass
    
    # Project onto normal modes: ΔQ_i = L_i^T · Δq (mass-weighted)
    # eigenvectors are already mass-weighted normal modes
    delta_Q = eigenvectors.T @ delta_q_mw
    
    # Calculate reorganization energies
    # λ_i = (1/2) * ω_i^2 * ΔQ_i^2 (in atomic units)
    # Convert frequency from cm^-1 to atomic units
    omega_au = frequencies / HARTREE_TO_CM  # in Hartree
    
    # Reorganization energy in eV
    # The eigenvalue from the Hessian is already in atomic units (Hartree/Bohr^2/amu)
    # ω_au is angular frequency in atomic units
    lambda_i = 0.5 * (omega_au ** 2) * (delta_Q ** 2) * HARTREE_TO_EV
    
    # Huang-Rhys factors: S_i = λ_i / (ℏω_i)
    omega_ev = frequencies * CM_TO_EV
    huang_rhys = np.zeros_like(lambda_i)
    valid_modes = np.abs(omega_ev) > 1e-6
    huang_rhys[valid_modes] = lambda_i[valid_modes] / np.abs(omega_ev[valid_modes])
    
    # Total reorganization energy
    lambda_total = np.sum(lambda_i[frequencies > 50])  # Skip translation/rotation
    
    print(f"\n{'='*60}")
    print("SPECTRAL DENSITY RESULTS")
    print(f"{'='*60}")
    print(f"Total reorganization energy: {lambda_total*1000:.2f} meV ({lambda_total:.4f} eV)")
    print(f"                           = {lambda_total/CM_TO_EV:.2f} cm^-1")
    print(f"{'='*60}")
    
    # Print mode-by-mode results
    print(f"\n{'Mode':>5} {'Freq (cm-1)':>12} {'λ (meV)':>12} {'S':>10} {'ΔQ':>12}")
    print("-" * 55)
    for i in range(n_modes):
        if frequencies[i] > 50:  # Skip translation/rotation
            print(f"{i+1:5d} {frequencies[i]:12.2f} {lambda_i[i]*1000:12.4f} {huang_rhys[i]:10.4f} {delta_Q[i]:12.6f}")
    
    # Save results
    results = {
        'frequencies_cm': frequencies,
        'frequencies_ev': omega_ev,
        'reorganization_energy_ev': lambda_i,
        'reorganization_energy_mev': lambda_i * 1000,
        'huang_rhys_factor': huang_rhys,
        'delta_Q': delta_Q,
        'total_reorganization_energy_ev': lambda_total,
    }
    
    np.savez(output_dir / 'spectral_density.npz', **results)
    print(f"\nSaved: {output_dir / 'spectral_density.npz'}")
    
    # Write text file
    with open(output_dir / 'spectral_density.txt', 'w') as f:
        f.write("# Spectral Density for Pentacene (Neutral -> Cation)\n")
        f.write(f"# Total reorganization energy: {lambda_total*1000:.4f} meV\n")
        f.write("#\n")
        f.write("# Mode  Freq(cm-1)  Lambda(meV)  Huang-Rhys   Delta_Q\n")
        for i in range(n_modes):
            if frequencies[i] > 50:
                f.write(f"{i+1:5d} {frequencies[i]:12.4f} {lambda_i[i]*1000:12.6f} {huang_rhys[i]:12.6f} {delta_Q[i]:12.6f}\n")
    
    print(f"Saved: {output_dir / 'spectral_density.txt'}")
    
    # Create plots
    plot_spectral_density(frequencies, lambda_i, huang_rhys, output_dir)
    
    return results


def plot_spectral_density(frequencies, lambda_i, huang_rhys, output_dir):
    """Create plots of spectral density"""
    
    # Filter valid modes (> 50 cm^-1)
    valid = frequencies > 50
    freq_valid = frequencies[valid]
    lambda_valid = lambda_i[valid] * 1000  # meV
    hr_valid = huang_rhys[valid]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Discrete spectral density (stick spectrum)
    ax1 = axes[0, 0]
    ax1.stem(freq_valid, lambda_valid, linefmt='b-', markerfmt='bo', basefmt=' ')
    ax1.set_xlabel('Frequency (cm$^{-1}$)')
    ax1.set_ylabel('Reorganization Energy (meV)')
    ax1.set_title('Discrete Spectral Density')
    ax1.set_xlim(0, max(freq_valid) * 1.1)
    
    # 2. Huang-Rhys factors
    ax2 = axes[0, 1]
    ax2.stem(freq_valid, hr_valid, linefmt='r-', markerfmt='ro', basefmt=' ')
    ax2.set_xlabel('Frequency (cm$^{-1}$)')
    ax2.set_ylabel('Huang-Rhys Factor S')
    ax2.set_title('Huang-Rhys Factors')
    ax2.set_xlim(0, max(freq_valid) * 1.1)
    
    # 3. Broadened spectral density (Lorentzian)
    ax3 = axes[1, 0]
    omega_grid = np.linspace(0, max(freq_valid) * 1.1, 2000)
    gamma = 10  # Broadening width in cm^-1
    
    J_omega = np.zeros_like(omega_grid)
    for w, lam in zip(freq_valid, lambda_valid):
        J_omega += lam * gamma / ((omega_grid - w)**2 + gamma**2) / np.pi
    
    ax3.plot(omega_grid, J_omega, 'b-', linewidth=1.5)
    ax3.fill_between(omega_grid, J_omega, alpha=0.3)
    ax3.set_xlabel('Frequency (cm$^{-1}$)')
    ax3.set_ylabel('J(ω) (meV)')
    ax3.set_title(f'Broadened Spectral Density (γ = {gamma} cm$^{{-1}}$)')
    ax3.set_xlim(0, max(freq_valid) * 1.1)
    
    # 4. Cumulative reorganization energy
    ax4 = axes[1, 1]
    sorted_idx = np.argsort(freq_valid)
    freq_sorted = freq_valid[sorted_idx]
    lambda_sorted = lambda_valid[sorted_idx]
    lambda_cumulative = np.cumsum(lambda_sorted)
    
    ax4.plot(freq_sorted, lambda_cumulative, 'g-', linewidth=2)
    ax4.axhline(lambda_cumulative[-1], color='gray', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Frequency (cm$^{-1}$)')
    ax4.set_ylabel('Cumulative λ (meV)')
    ax4.set_title(f'Cumulative Reorganization Energy (Total: {lambda_cumulative[-1]:.2f} meV)')
    ax4.set_xlim(0, max(freq_valid) * 1.1)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'spectral_density.png', dpi=150)
    plt.savefig(output_dir / 'spectral_density.pdf')
    print(f"Saved: {output_dir / 'spectral_density.png'}")
    print(f"Saved: {output_dir / 'spectral_density.pdf'}")
    plt.close()
    
    # Additional plot: Low-frequency region
    fig2, ax = plt.subplots(figsize=(10, 5))
    low_freq_mask = freq_valid < 500
    if np.any(low_freq_mask):
        ax.stem(freq_valid[low_freq_mask], lambda_valid[low_freq_mask], 
                linefmt='b-', markerfmt='bo', basefmt=' ')
        ax.set_xlabel('Frequency (cm$^{-1}$)')
        ax.set_ylabel('Reorganization Energy (meV)')
        ax.set_title('Low-Frequency Spectral Density (< 500 cm$^{-1}$)')
        plt.tight_layout()
        plt.savefig(output_dir / 'spectral_density_lowfreq.png', dpi=150)
        print(f"Saved: {output_dir / 'spectral_density_lowfreq.png'}")
        plt.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Calculate spectral density from neutral/cation normal modes')
    parser.add_argument('--neutral', default='neutral', help='Directory with neutral calculation')
    parser.add_argument('--cation', default='cation', help='Directory with cation calculation')
    parser.add_argument('--output', default='.', help='Output directory')
    args = parser.parse_args()
    
    calculate_spectral_density(args.neutral, args.cation, args.output)


if __name__ == '__main__':
    main()
