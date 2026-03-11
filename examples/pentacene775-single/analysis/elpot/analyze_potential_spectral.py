#!/usr/bin/env python3
"""
Analyze displacement potential (φ_disp) fluctuations and compute spectral density
for each QM atom.

Read potential_phi_disp.dat and compute:
1. Autocorrelation function C(t) for φ_disp of each atom
2. Spectral density J(ω) from Fourier transform of C(t)

The relation between correlation function and spectral density:
    J(ν̃) = (2πc ν̃ / k_B T) ∫_0^∞ dt C_cl(t) cos(2πc ν̃ t)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.integrate import simpson
from pathlib import Path

# =============================================================================
# Physical constants
# =============================================================================
kB_cm = 0.695034800  # cm⁻¹/K (Boltzmann constant)
c_cm_fs = 2.99792458e-5  # cm/fs (speed of light)
two_pi_c = 2 * np.pi * c_cm_fs  # rad·cm/fs

# Unit conversions
eV_to_cm = 8065.5  # cm⁻¹/eV
cm_to_eV = 1.23984e-4  # eV/cm⁻¹
V_to_eV = 1.0  # V = eV/e, so potential in V is energy in eV when multiplied by e


# =============================================================================
# File I/O
# =============================================================================
def read_potential_data(filepath: str) -> tuple:
    """
    Read potential_phi_disp.dat file.
    
    Parameters
    ----------
    filepath : str
        Path to the potential data file
        
    Returns
    -------
    data : np.ndarray
        Potential data, shape (n_frames, n_atoms), in Volt
    n_atoms : int
        Number of QM atoms
    """
    data = np.loadtxt(filepath, comments='#')
    n_frames, n_atoms = data.shape
    print(f"  Loaded {n_frames} frames, {n_atoms} atoms")
    return data, n_atoms


# =============================================================================
# Correlation function computation
# =============================================================================
def compute_correlation_fft(delta: np.ndarray, dt: float) -> tuple:
    """
    Compute autocorrelation function using FFT (Wiener-Khinchin theorem).
    
    C(τ) = <δX(t) δX(t+τ)>
    
    Parameters
    ----------
    delta : np.ndarray
        Fluctuation time series (mean subtracted)
    dt : float
        Time step (fs)
        
    Returns
    -------
    t_corr : np.ndarray
        Correlation time array (fs)
    C : np.ndarray
        Correlation function (same units as δX²)
    """
    N = len(delta)
    
    # Zero-padding for linear (not circular) correlation
    delta_padded = np.concatenate([delta, np.zeros(N)])
    
    # FFT
    fft_delta = fft(delta_padded)
    
    # Power spectrum
    power_spectrum = np.abs(fft_delta)**2
    
    # Inverse FFT gives autocorrelation
    autocorr_full = np.real(ifft(power_spectrum))[:N]
    
    # Normalization: divide by number of pairs at each lag
    norm = np.arange(N, 0, -1)
    C = autocorr_full / norm
    
    t_corr = np.arange(N) * dt
    
    return t_corr, C


def compute_correlation_segmented(delta: np.ndarray, dt: float, 
                                   segment_length: int, corr_length: int) -> tuple:
    """
    Compute correlation function with segment averaging.
    """
    N = len(delta)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        segment = delta[start:end]
        
        # Subtract segment mean
        segment = segment - np.mean(segment)
        
        _, C_seg = compute_correlation_fft(segment, dt)
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
    
    Parameters
    ----------
    t_corr : np.ndarray
        Correlation time array (fs)
    C : np.ndarray
        Correlation function (cm⁻²)
    nu_out : np.ndarray
        Output wavenumber array (cm⁻¹)
    T : float
        Temperature (K)
    use_window : bool
        If True, apply Hanning window to reduce spectral leakage
        
    Returns
    -------
    J : np.ndarray
        Spectral density (cm⁻¹)
    """
    dt = t_corr[1] - t_corr[0]
    
    # Apply window function to reduce spectral leakage
    if use_window:
        window = 0.5 * (1 + np.cos(np.pi * t_corr / t_corr[-1]))
        C_windowed = C * window
    else:
        C_windowed = C
    
    J = np.zeros(len(nu_out))
    
    for i, nu in enumerate(nu_out):
        omega = two_pi_c * nu  # rad/fs
        integrand = C_windowed * np.cos(omega * t_corr)
        integral = simpson(integrand, dx=dt)
        J[i] = two_pi_c * nu / (2.0 * kB_cm * T) * integral
    
    return J


# =============================================================================
# Analysis
# =============================================================================
def analyze_potential_spectral(
    data: np.ndarray,
    dt: float,
    T: float,
    segment_ps: float = 40.0,
    corr_ps: float = 30.0
) -> dict:
    """
    Analyze potential fluctuations and compute spectral density for each atom.
    
    Parameters
    ----------
    data : np.ndarray
        Potential data (n_frames, n_atoms) in Volt
    dt : float
        Time step (fs)
    T : float
        Temperature (K)
    segment_ps : float
        Segment length for averaging (ps)
    corr_ps : float
        Correlation length (ps)
        
    Returns
    -------
    results : dict
        Analysis results
    """
    n_frames, n_atoms = data.shape
    total_time_ps = n_frames * dt / 1000
    
    print(f"\n=== Potential Statistics ===")
    print(f"  Frames: {n_frames}, Atoms: {n_atoms}")
    print(f"  Total time: {total_time_ps:.1f} ps")
    print(f"  Time step: {dt} fs")
    
    # Convert V to cm⁻¹ for spectral analysis
    # V = eV/e, so potential in V equals energy in eV when considering e
    # eV -> cm⁻¹: multiply by eV_to_cm
    data_cm = data * eV_to_cm  # V -> cm⁻¹
    
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    n_segments_expected = int(total_time_ps / segment_ps)
    
    print(f"\n=== Analysis Parameters ===")
    print(f"  Segment: {segment_ps} ps, Correlation length: {corr_ps} ps")
    print(f"  Expected segments: {n_segments_expected}")
    
    # Spectral density output grid
    nu_max = 2000  # cm⁻¹
    nu_out = np.linspace(1, nu_max, 500)
    
    # Storage
    t_corr = None
    C_all = []
    SD_all = []
    lambda_all = []
    
    kBT = kB_cm * T
    
    print(f"\n=== Processing {n_atoms} atoms ===")
    
    for atom_idx in range(n_atoms):
        potential_cm = data_cm[:, atom_idx]
        
        # Fluctuation
        delta = potential_cm - np.mean(potential_cm)
        
        if n_segments_expected >= 2:
            t_corr, C, _ = compute_correlation_segmented(delta, dt, segment_length, corr_length)
        else:
            t_corr, C = compute_correlation_fft(delta, dt)
            max_corr = min(len(t_corr), corr_length)
            t_corr = t_corr[:max_corr]
            C = C[:max_corr]
        
        # Spectral density
        SD = correlation_to_spectral_density(t_corr, C, nu_out, T)
        
        # Reorganization energy λ = C(0) / (2 k_B T)
        lambda_cm = C[0] / (2 * kBT)
        
        C_all.append(C)
        SD_all.append(SD)
        lambda_all.append(lambda_cm)
        
        if (atom_idx + 1) % 10 == 0:
            print(f"  Atom {atom_idx + 1}/{n_atoms} done")
    
    print(f"  All atoms processed.")
    
    C_all = np.array(C_all)  # (n_atoms, corr_length)
    SD_all = np.array(SD_all)  # (n_atoms, n_freq)
    lambda_all = np.array(lambda_all)  # (n_atoms,)
    
    # Sum and mean
    SD_sum = np.sum(SD_all, axis=0)
    SD_mean = np.mean(SD_all, axis=0)
    
    print(f"\n=== Reorganization Energy ===")
    print(f"  λ per atom (cm⁻¹): min={lambda_all.min():.2f}, max={lambda_all.max():.2f}, mean={lambda_all.mean():.2f}")
    print(f"  λ total (sum): {lambda_all.sum():.2f} cm⁻¹ = {lambda_all.sum() * cm_to_eV:.4f} eV")
    
    return {
        'n_atoms': n_atoms,
        'n_frames': n_frames,
        'dt': dt,
        'T': T,
        't_corr': t_corr,
        'nu': nu_out,
        'C_all': C_all,
        'SD_all': SD_all,
        'SD_sum': SD_sum,
        'SD_mean': SD_mean,
        'lambda_all': lambda_all,
        'lambda_total': lambda_all.sum(),
        'corr_ps': corr_ps,
    }


def plot_results(results: dict, output_dir: str = '.'):
    """Plot spectral densities for each atom."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    n_atoms = results['n_atoms']
    nu = results['nu']
    SD_all = results['SD_all']
    SD_sum = results['SD_sum']
    
    # 1. Individual atom spectra (grid)
    ncols = 6
    nrows = (n_atoms + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 3 * nrows))
    axes = axes.flatten()
    
    for i in range(n_atoms):
        ax = axes[i]
        ax.plot(nu, SD_all[i], 'b-', linewidth=0.8)
        ax.set_xlim(0, 2000)
        ax.set_ylim(0, None)
        ax.set_title(f'Atom {i}', fontsize=10)
        ax.set_xlabel(r'$\tilde{\nu}$ (cm⁻¹)', fontsize=8)
        ax.set_ylabel(r'$J(\tilde{\nu})$', fontsize=8)
        ax.tick_params(labelsize=7)
    
    # Hide unused axes
    for i in range(n_atoms, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    output_file = output_dir / 'spectral_density_atoms.png'
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()
    
    # 2. Sum and comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: All atoms overlaid
    ax = axes[0]
    for i in range(n_atoms):
        ax.plot(nu, SD_all[i], alpha=0.5, linewidth=0.8, label=f'{i}' if i < 5 else None)
    ax.set_xlim(0, 2000)
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel(r'$J(\tilde{\nu})$ (cm⁻¹)', fontsize=12)
    ax.set_title('Spectral Density (All Atoms)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Right: Sum over all atoms
    ax = axes[1]
    ax.plot(nu, SD_sum, 'r-', linewidth=1.5)
    ax.set_xlim(0, 2000)
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel(r'$J_{total}(\tilde{\nu})$ (cm⁻¹)', fontsize=12)
    ax.set_title('Total Spectral Density (Sum)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / 'spectral_density_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()
    
    # 3. Reorganization energy per atom
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(range(n_atoms), results['lambda_all'], color='steelblue')
    ax.set_xlabel('Atom Index', fontsize=12)
    ax.set_ylabel(r'$\lambda$ (cm⁻¹)', fontsize=12)
    ax.set_title('Reorganization Energy per Atom', fontsize=13, fontweight='bold')
    ax.set_xlim(-0.5, n_atoms - 0.5)
    plt.tight_layout()
    output_file = output_dir / 'reorganization_energy_atoms.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def save_results(results: dict, output_dir: str = '.'):
    """Save numerical results to files."""
    output_dir = Path(output_dir)
    
    # Spectral density for each atom
    header = f"# Spectral density J(nu) for {results['n_atoms']} atoms\n"
    header += "# Column 0: wavenumber (cm-1)\n"
    header += "# Columns 1-N: J(nu) for atoms 0 to N-1 (cm-1)\n"
    
    data_out = np.column_stack([results['nu'], results['SD_all'].T])
    np.savetxt(output_dir / 'spectral_density_atoms.dat', data_out, 
               header=header, fmt='%.6e')
    
    # Sum spectral density
    np.savetxt(output_dir / 'spectral_density_sum.dat', 
               np.column_stack([results['nu'], results['SD_sum']]),
               header="# Wavenumber (cm-1)  J_total (cm-1)", fmt='%.6e')
    
    # Reorganization energies
    np.savetxt(output_dir / 'reorganization_energy.dat',
               np.column_stack([np.arange(results['n_atoms']), 
                                results['lambda_all'],
                                results['lambda_all'] * cm_to_eV]),
               header="# Atom  lambda(cm-1)  lambda(eV)", fmt=['%d', '%.4f', '%.6f'])
    
    print(f"  Data saved to {output_dir}")


# =============================================================================
# Main
# =============================================================================
# =============================================================================
# Configuration
# =============================================================================
INPUT_FILE = 'potential_phi_disp.dat'
DT = 4.0           # Time step (fs)
TEMPERATURE = 300  # Temperature (K)
SEGMENT_PS = 40.0  # Segment length for averaging (ps)
CORR_PS = 30.0     # Correlation length (ps)
OUTPUT_DIR = '.'   # Output directory


def main():
    print("=" * 70)
    print("Potential Spectral Density Analysis")
    print("=" * 70)
    
    # Read data
    print(f"\nReading: {INPUT_FILE}")
    data, n_atoms = read_potential_data(INPUT_FILE)
    
    # Analyze
    results = analyze_potential_spectral(
        data, DT, TEMPERATURE,
        segment_ps=SEGMENT_PS,
        corr_ps=CORR_PS
    )
    
    # Plot
    print(f"\n=== Generating plots ===")
    plot_results(results, OUTPUT_DIR)
    
    # Save data
    save_results(results, OUTPUT_DIR)
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
