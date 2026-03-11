#!/usr/bin/env python3
"""
Compute spectral densities from energy_diff.dat with multiple molecule columns.

For each molecule column, calculates:
- Autocorrelation function
- Spectral density J(ω)
- Reorganization energy λ

Plots all spectral densities together for comparison.
"""

import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from pathlib import Path

# =============================================================================
# Physical constants
# =============================================================================
kB_cm = 0.695034800  # cm⁻¹/K (Boltzmann constant)
c_cm_fs = 2.99792458e-5  # cm/fs (speed of light)
two_pi_c = 2 * np.pi * c_cm_fs  # rad·cm/fs

# Unit conversions
Ha_to_cm = 219474.63  # cm⁻¹/Ha
Ha_to_eV = 27.211386245988  # eV/Ha
eV_to_cm = 8065.54  # cm⁻¹/eV
cm_to_eV = 1.23984e-4  # eV/cm⁻¹


# =============================================================================
# File I/O
# =============================================================================
def read_energy_diff_multi(filepath: str) -> dict:
    """Read energy_diff.dat file with multiple molecule columns.
    
    Returns:
        dict with keys:
            'frame': array of frame numbers
            'time': array of time values
            'energies': 2D array (n_frames, n_molecules) in eV
            'col_names': list of column names
    """
    frames = []
    times = []
    energies = []
    col_names = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            if not line_stripped:
                continue
            if line_stripped.startswith('#'):
                # Parse column names from header
                if 'dE_mol' in line or 'E_mol' in line:
                    parts = line.replace('#', '').split()
                    col_names = [p for p in parts if p.startswith('dE_mol') or p.startswith('E_mol')]
                continue
            
            parts = line_stripped.split()
            if len(parts) >= 3:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                energy_vals = [float(c) for c in parts[2:]]
                energies.append(energy_vals)
    
    energies = np.array(energies)
    n_cols = energies.shape[1] if len(energies.shape) > 1 else 1
    
    # Generate column names if not found
    if not col_names or len(col_names) != n_cols:
        col_names = [f'mol{i+1}' for i in range(n_cols)]
    else:
        # Clean up column names (remove (eV) suffix)
        col_names = [name.replace('(eV)', '').replace('dE_', '').replace('E_', '') for name in col_names]
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'energies': energies,  # eV
        'col_names': col_names,
    }


def interpolate_nan_multi(data: dict) -> dict:
    """Interpolate NaN values in the energies array (multi-column)."""
    energies = data['energies']
    time = data['time']
    n_cols = energies.shape[1]
    
    for j in range(n_cols):
        col = energies[:, j]
        nan_mask = np.isnan(col)
        n_nan = np.sum(nan_mask)
        
        if n_nan == 0:
            continue
        
        print(f"  Column {j}: Found {n_nan} NaN values, interpolating...")
        
        valid_mask = ~nan_mask
        if np.sum(valid_mask) < 2:
            raise ValueError(f"Column {j}: Not enough valid data points for interpolation")
        
        f_interp = interp1d(time[valid_mask], col[valid_mask], 
                            kind='linear', fill_value='extrapolate')
        energies[nan_mask, j] = f_interp(time[nan_mask])
    
    data['energies'] = energies
    return data


# =============================================================================
# Correlation function computation
# =============================================================================
def compute_correlation_fft(delta: np.ndarray, dt: float) -> tuple:
    """Compute autocorrelation function using FFT."""
    N = len(delta)
    delta_padded = np.concatenate([delta, np.zeros(N)])
    fft_delta = fft(delta_padded)
    power_spectrum = np.abs(fft_delta)**2
    autocorr_full = np.real(ifft(power_spectrum))[:N]
    norm = np.arange(N, 0, -1)
    C = autocorr_full / norm
    t_corr = np.arange(N) * dt
    return t_corr, C


def compute_cross_correlation_fft(delta1: np.ndarray, delta2: np.ndarray, dt: float) -> tuple:
    """
    Compute cross-correlation function using FFT.
    C_12(τ) = <δX_1(t) δX_2(t+τ)>
    """
    N = len(delta1)
    delta1_padded = np.concatenate([delta1, np.zeros(N)])
    delta2_padded = np.concatenate([delta2, np.zeros(N)])
    
    fft_delta1 = fft(delta1_padded)
    fft_delta2 = fft(delta2_padded)
    
    # Cross spectrum: conj(FFT1) * FFT2
    cross_spectrum = np.conj(fft_delta1) * fft_delta2
    cross_corr_full = np.real(ifft(cross_spectrum))[:N]
    
    norm = np.arange(N, 0, -1)
    C = cross_corr_full / norm
    t_corr = np.arange(N) * dt
    
    return t_corr, C


def compute_correlation_segmented(delta: np.ndarray, dt: float, 
                                   segment_length: int, corr_length: int) -> tuple:
    """Compute correlation function with segment averaging."""
    N = len(delta)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        segment = delta[start:end]
        segment = segment - np.mean(segment)
        _, C_seg = compute_correlation_fft(segment, dt)
        C_sum += C_seg[:corr_length]
    
    C_avg = C_sum / n_segments
    t_corr = np.arange(corr_length) * dt
    
    return t_corr, C_avg, n_segments


def compute_cross_correlation_segmented(delta1: np.ndarray, delta2: np.ndarray, dt: float, 
                                         segment_length: int, corr_length: int) -> tuple:
    """Compute cross-correlation function with segment averaging."""
    N = min(len(delta1), len(delta2))
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        seg1 = delta1[start:end] - np.mean(delta1[start:end])
        seg2 = delta2[start:end] - np.mean(delta2[start:end])
        _, C_seg = compute_cross_correlation_fft(seg1, seg2, dt)
        C_sum += C_seg[:corr_length]
    
    C_avg = C_sum / n_segments
    t_corr = np.arange(corr_length) * dt
    
    return t_corr, C_avg, n_segments


# =============================================================================
# Spectral density computation
# =============================================================================
def correlation_to_spectral_density(t_corr: np.ndarray, C: np.ndarray, 
                                     nu_out: np.ndarray, T: float,
                                     use_window: bool = True) -> np.ndarray:
    """
    Compute spectral density from correlation function.
    J(ν̃) = (2πc ν̃ / k_B T) ∫_0^∞ dt C_cl(t) cos(2πc ν̃ t)
    """
    dt = t_corr[1] - t_corr[0]
    
    if use_window:
        window = 0.5 * (1 + np.cos(np.pi * t_corr / t_corr[-1]))
        C_windowed = C * window
    else:
        C_windowed = C
    
    J = np.zeros(len(nu_out))
    
    for i, nu in enumerate(nu_out):
        omega = two_pi_c * nu
        integrand = C_windowed * np.cos(omega * t_corr)
        integral = simpson(integrand, dx=dt)
        J[i] = two_pi_c * nu / (2.0 * kB_cm * T) * integral
    
    return J


# =============================================================================
# Main analysis
# =============================================================================
def main():
    """Main analysis routine."""
    parser = argparse.ArgumentParser(description='Compute spectral densities from energy_diff.dat')
    parser.add_argument('--input', '-i', default='energy_diff.dat',
                        help='Input file (default: energy_diff.dat)')
    parser.add_argument('--dt', type=float, default=4.0,
                        help='Time step in fs (default: 4.0)')
    parser.add_argument('--temp', '-T', type=float, default=300.0,
                        help='Temperature in K (default: 300.0)')
    parser.add_argument('--segment', type=float, default=20.0,
                        help='Segment length in ps (default: 20.0)')
    parser.add_argument('--corr', type=float, default=10.0,
                        help='Correlation length in ps (default: 10.0)')
    parser.add_argument('--nu-max', type=float, default=2000,
                        help='Maximum frequency in cm^-1 (default: 2000)')
    parser.add_argument('--output', '-o', default='spectral_density',
                        help='Output prefix (default: spectral_density)')
    parser.add_argument('--molecules', '-m', type=str, default=None,
                        help='Comma-separated list of molecule IDs to plot (default: all)')
    args = parser.parse_args()
    
    dt = args.dt  # fs
    T = args.temp  # K
    segment_ps = args.segment
    corr_ps = args.corr
    nu_max = args.nu_max
    
    print(f"=== Spectral Density Analysis (Multi-column) ===")
    print(f"Input: {args.input}")
    print(f"Time step: {dt} fs, Temperature: {T} K")
    
    # Read data
    data = read_energy_diff_multi(args.input)
    data = interpolate_nan_multi(data)
    
    n_points = len(data['time'])
    n_cols = data['energies'].shape[1]
    col_names = data['col_names']
    
    # Filter molecules if specified
    if args.molecules:
        target_mols = [m.strip() for m in args.molecules.split(',')]
        target_mols_set = set(target_mols)
        
        # Find indices of matching molecules
        selected_indices = []
        selected_names = []
        for j, name in enumerate(col_names):
            # Extract molecule ID from name (e.g., 'mol128' -> '128')
            mol_id = name.replace('mol', '')
            if mol_id in target_mols_set or name in target_mols_set:
                selected_indices.append(j)
                selected_names.append(name)
        
        if not selected_indices:
            print(f"Warning: No matching molecules found for {args.molecules}")
            print(f"Available: {', '.join(col_names)}")
        else:
            # Filter data
            data['energies'] = data['energies'][:, selected_indices]
            col_names = selected_names
            n_cols = len(selected_indices)
            print(f"\nFiltered to {n_cols} molecules: {', '.join(col_names)}")
    
    print(f"\nUsing {n_points} frames, {n_cols} molecules")
    print(f"Molecules: {', '.join(col_names)}")
    
    # Convert to cm⁻¹ (input is in eV)
    energies_cm = data['energies'] * eV_to_cm
    
    # Compute fluctuations (subtract mean from each column)
    deltas = np.zeros_like(energies_cm)
    for j in range(n_cols):
        deltas[:, j] = energies_cm[:, j] - np.mean(energies_cm[:, j])
    
    # Analysis parameters
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    
    print(f"\nSegment length: {segment_length} frames ({segment_ps} ps)")
    print(f"Correlation length: {corr_length} frames ({corr_ps} ps)")
    
    # Compute correlation functions and spectral densities for each molecule
    nu_out = np.linspace(1, nu_max, 500)
    J_all = np.zeros((len(nu_out), n_cols))
    lambda_all = np.zeros(n_cols)
    kBT = kB_cm * T
    
    print(f"\nComputing spectral densities...")
    for j in range(n_cols):
        t_corr, C, n_seg = compute_correlation_segmented(deltas[:, j], dt, segment_length, corr_length)
        J = correlation_to_spectral_density(t_corr, C, nu_out, T)
        J_all[:, j] = J
        lambda_all[j] = C[0] / (2 * kBT)
    
    print(f"  {n_seg} segments used")
    
    # Print reorganization energies
    print(f"\nReorganization energies:")
    print(f"{'Molecule':<15} {'λ (cm⁻¹)':>12} {'λ (eV)':>12}")
    print("-" * 41)
    for j in range(n_cols):
        print(f"{col_names[j]:<15} {lambda_all[j]:>12.2f} {lambda_all[j]*cm_to_eV:>12.4f}")
    
    # Plot spectral densities
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Color map for distinguishing molecules
    colors = plt.cm.tab10(np.linspace(0, 1, min(n_cols, 10)))
    if n_cols > 10:
        colors = plt.cm.tab20(np.linspace(0, 1, n_cols))
    
    # Left: All spectral densities
    ax = axes[0]
    for j in range(n_cols):
        ax.plot(nu_out, J_all[:, j], linewidth=1.5, color=colors[j % len(colors)], 
                label=f'{col_names[j]} (λ={lambda_all[j]*cm_to_eV:.3f} eV)')
    ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=14)
    ax.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=14)
    ax.set_xlim(0, nu_max)
    ax.legend(fontsize=8, loc='upper right', ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_title('Spectral Densities', fontsize=13, fontweight='bold')
    
    # Right: Sum and average
    ax = axes[1]
    J_sum = np.sum(J_all, axis=1)
    J_avg = np.mean(J_all, axis=1)
    ax.plot(nu_out, J_sum, 'k-', linewidth=2, label=f'Sum (λ={np.sum(lambda_all)*cm_to_eV:.3f} eV)')
    ax.plot(nu_out, J_avg, 'r--', linewidth=2, label=f'Average (λ={np.mean(lambda_all)*cm_to_eV:.4f} eV)')
    ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=14)
    ax.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=14)
    ax.set_xlim(0, nu_max)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_title('Sum and Average', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{args.output}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.output}.pdf', bbox_inches='tight')
    print(f"\nFigures saved: {args.output}.png, .pdf")
    plt.close()
    
    # Plot individual spectral densities in subplots
    n_rows = int(np.ceil(n_cols / 4))
    fig, axes = plt.subplots(n_rows, 4, figsize=(16, 3*n_rows))
    axes = np.atleast_2d(axes)
    
    for j in range(n_cols):
        ax = axes[j // 4, j % 4]
        ax.plot(nu_out, J_all[:, j], 'b-', linewidth=1)
        ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=10)
        ax.set_ylabel(r'$J(\omega)$', fontsize=10)
        ax.set_xlim(0, nu_max)
        ax.set_title(f'{col_names[j]}\nλ={lambda_all[j]*cm_to_eV:.4f} eV', fontsize=10)
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for j in range(n_cols, n_rows * 4):
        axes[j // 4, j % 4].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{args.output}_individual.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.output}_individual.pdf', bbox_inches='tight')
    print(f"Figures saved: {args.output}_individual.png, .pdf")
    plt.close()
    
    # Save data
    with open(f'{args.output}.dat', 'w') as f:
        f.write("# Spectral density analysis\n")
        f.write(f"# Input: {args.input}\n")
        f.write(f"# Temperature: {T} K\n")
        f.write("# Reorganization energies:\n")
        for j in range(n_cols):
            f.write(f"#   {col_names[j]}: λ = {lambda_all[j]:.2f} cm⁻¹ = {lambda_all[j]*cm_to_eV:.4f} eV\n")
        f.write(f"# Sum: λ = {np.sum(lambda_all):.2f} cm⁻¹ = {np.sum(lambda_all)*cm_to_eV:.4f} eV\n")
        f.write(f"# Average: λ = {np.mean(lambda_all):.2f} cm⁻¹ = {np.mean(lambda_all)*cm_to_eV:.4f} eV\n")
        
        # Header line
        header = "# nu(cm-1)"
        for name in col_names:
            header += f"  J_{name}"
        header += "  J_sum  J_avg\n"
        f.write(header)
        
        for i in range(len(nu_out)):
            line = f"{nu_out[i]:10.2f}"
            for j in range(n_cols):
                line += f" {J_all[i, j]:14.6e}"
            line += f" {J_sum[i]:14.6e} {J_avg[i]:14.6e}\n"
            f.write(line)
    print(f"Data saved: {args.output}.dat")
    
    print("\n=== Analysis complete ===")


if __name__ == '__main__':
    main()
