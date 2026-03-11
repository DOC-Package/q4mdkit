#!/usr/bin/env python3
"""
Analyze Mulliken charge difference fluctuations and compute spectral density
for each QM atom.

Read mulliken_charges_diff.dat and compute:
1. Autocorrelation function C(t) for charge fluctuations of each atom
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


# =============================================================================
# Configuration
# =============================================================================
INPUT_FILE = 'mulliken_charges_diff.dat'
DT = 4.0           # Time step (fs)
TEMPERATURE = 300  # Temperature (K)
SEGMENT_PS = 40.0  # Segment length for averaging (ps)
CORR_PS = 30.0     # Correlation length (ps)
OUTPUT_DIR = '.'   # Output directory


# =============================================================================
# File I/O
# =============================================================================
def read_charge_diff_data(filepath: str) -> tuple:
    """
    Read mulliken_charges_diff.dat file.
    
    Parameters
    ----------
    filepath : str
        Path to the charge difference data file
        
    Returns
    -------
    frames : np.ndarray
        Frame indices
    times : np.ndarray
        Time values (fs)
    charges : np.ndarray
        Charge difference data, shape (n_frames, n_atoms)
    atom_types : list or None
        Atom types if found in header
    """
    frames = []
    times = []
    charges = []
    atom_types = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line_strip = line.strip()
            
            # Parse header for atom types
            if line_strip.startswith('# Atom types:'):
                atom_types = line_strip.replace('# Atom types:', '').split()
                continue
            
            # Skip other comments
            if not line_strip or line_strip.startswith('#'):
                continue
            
            parts = line_strip.split()
            if len(parts) >= 3:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                # Columns 2:-1 are the per-atom charges, last column is Total
                charges.append([float(x) for x in parts[2:-1]])
    
    frames = np.array(frames)
    times = np.array(times)
    charges = np.array(charges)
    
    n_atoms = charges.shape[1]
    print(f"  Loaded {filepath}: {len(frames)} frames, {n_atoms} atoms")
    
    # Infer dt if not uniform
    if len(times) > 1:
        dt_inferred = times[1] - times[0]
        print(f"  Inferred dt: {dt_inferred} fs")
    
    return frames, times, charges, atom_types


def interpolate_nan(charges: np.ndarray, times: np.ndarray) -> np.ndarray:
    """
    Interpolate NaN values in the charges array using linear interpolation.
    
    Parameters
    ----------
    charges : np.ndarray
        Charge data (n_frames, n_atoms)
    times : np.ndarray
        Time array (n_frames,)
        
    Returns
    -------
    charges_interp : np.ndarray
        Charge data with NaN values interpolated
    """
    from scipy.interpolate import interp1d
    
    charges_interp = charges.copy()
    n_frames, n_atoms = charges.shape
    
    total_nan = np.sum(np.isnan(charges))
    if total_nan == 0:
        print("  No NaN values found in data.")
        return charges_interp
    
    print(f"  Found {total_nan} NaN values, interpolating...")
    
    for atom_idx in range(n_atoms):
        col = charges[:, atom_idx]
        nan_mask = np.isnan(col)
        n_nan = np.sum(nan_mask)
        
        if n_nan == 0:
            continue
        
        valid_mask = ~nan_mask
        if np.sum(valid_mask) < 2:
            print(f"    Warning: Atom {atom_idx} has insufficient valid data for interpolation")
            continue
        
        # Create interpolation function from valid data
        f_interp = interp1d(times[valid_mask], col[valid_mask], 
                            kind='linear', fill_value='extrapolate')
        
        # Interpolate NaN values
        charges_interp[nan_mask, atom_idx] = f_interp(times[nan_mask])
    
    print(f"  Interpolation complete.")
    return charges_interp


# =============================================================================
# Correlation function computation
# =============================================================================
def compute_correlation_fft(delta: np.ndarray, dt: float) -> tuple:
    """
    Compute autocorrelation function using FFT (Wiener-Khinchin theorem).
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
    
    Note: For charge fluctuations, we need to convert units appropriately.
    The correlation function C has units of e², and we want J in cm⁻¹.
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
def analyze_charge_spectral(
    charges: np.ndarray,
    dt: float,
    T: float,
    segment_ps: float = 40.0,
    corr_ps: float = 30.0
) -> dict:
    """
    Analyze charge fluctuations and compute spectral density for each atom.
    
    Parameters
    ----------
    charges : np.ndarray
        Charge difference data (n_frames, n_atoms) in elementary charge
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
    n_frames, n_atoms = charges.shape
    total_time_ps = n_frames * dt / 1000
    
    print(f"\n=== Charge Difference Statistics ===")
    print(f"  Frames: {n_frames}, Atoms: {n_atoms}")
    print(f"  Total time: {total_time_ps:.1f} ps")
    print(f"  Time step: {dt} fs")
    
    # For charge fluctuations, we don't convert units - 
    # the spectral density will be in units related to e²
    # To get energy-like spectral density, we'd need a coupling constant
    
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
    variance_all = []
    
    kBT = kB_cm * T
    
    print(f"\n=== Processing {n_atoms} atoms ===")
    
    for atom_idx in range(n_atoms):
        charge = charges[:, atom_idx]
        
        # Fluctuation
        delta = charge - np.mean(charge)
        variance = np.var(delta)
        variance_all.append(variance)
        
        if n_segments_expected >= 2:
            t_corr, C, _ = compute_correlation_segmented(delta, dt, segment_length, corr_length)
        else:
            t_corr, C = compute_correlation_fft(delta, dt)
            max_corr = min(len(t_corr), corr_length)
            t_corr = t_corr[:max_corr]
            C = C[:max_corr]
        
        # Spectral density
        SD = correlation_to_spectral_density(t_corr, C, nu_out, T)
        
        C_all.append(C)
        SD_all.append(SD)
        
        if (atom_idx + 1) % 10 == 0:
            print(f"  Atom {atom_idx + 1}/{n_atoms} done")
    
    print(f"  All atoms processed.")
    
    C_all = np.array(C_all)  # (n_atoms, corr_length)
    SD_all = np.array(SD_all)  # (n_atoms, n_freq)
    variance_all = np.array(variance_all)  # (n_atoms,)
    
    # Sum and mean
    SD_sum = np.sum(SD_all, axis=0)
    SD_mean = np.mean(SD_all, axis=0)
    
    print(f"\n=== Variance Statistics ===")
    print(f"  Variance per atom (e²): min={variance_all.min():.6f}, max={variance_all.max():.6f}, mean={variance_all.mean():.6f}")
    
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
        'variance_all': variance_all,
        'corr_ps': corr_ps,
    }


def plot_results(results: dict, output_dir: str = '.', atom_types: list = None):
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
        label = f'{atom_types[i]}{i}' if atom_types else f'Atom {i}'
        ax.set_title(label, fontsize=10)
        ax.set_xlabel(r'$\tilde{\nu}$ (cm⁻¹)', fontsize=8)
        ax.set_ylabel(r'$J(\tilde{\nu})$', fontsize=8)
        ax.tick_params(labelsize=7)
    
    # Hide unused axes
    for i in range(n_atoms, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    output_file = output_dir / 'charge_spectral_density_atoms.png'
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()
    
    # 2. Sum and comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: All atoms overlaid
    ax = axes[0]
    for i in range(n_atoms):
        ax.plot(nu, SD_all[i], alpha=0.5, linewidth=0.8)
    ax.set_xlim(0, 2000)
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel(r'$J(\tilde{\nu})$', fontsize=12)
    ax.set_title('Charge Spectral Density (All Atoms)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Right: Sum over all atoms
    ax = axes[1]
    ax.plot(nu, SD_sum, 'r-', linewidth=1.5)
    ax.set_xlim(0, 2000)
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel(r'$J_{total}(\tilde{\nu})$', fontsize=12)
    ax.set_title('Total Charge Spectral Density (Sum)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / 'charge_spectral_density_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()
    
    # 3. Variance per atom
    fig, ax = plt.subplots(figsize=(10, 4))
    colors = ['C0' if t == 'C' else 'C1' for t in atom_types] if atom_types else 'steelblue'
    ax.bar(range(n_atoms), results['variance_all'], color=colors)
    ax.set_xlabel('Atom Index', fontsize=12)
    ax.set_ylabel(r'Variance ($e^2$)', fontsize=12)
    ax.set_title('Charge Fluctuation Variance per Atom', fontsize=13, fontweight='bold')
    ax.set_xlim(-0.5, n_atoms - 0.5)
    plt.tight_layout()
    output_file = output_dir / 'charge_variance_atoms.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def save_results(results: dict, output_dir: str = '.'):
    """Save numerical results to files."""
    output_dir = Path(output_dir)
    
    # Spectral density for each atom
    header = f"# Charge spectral density J(nu) for {results['n_atoms']} atoms\n"
    header += "# Column 0: wavenumber (cm-1)\n"
    header += "# Columns 1-N: J(nu) for atoms 0 to N-1\n"
    
    data_out = np.column_stack([results['nu'], results['SD_all'].T])
    np.savetxt(output_dir / 'charge_spectral_density_atoms.dat', data_out, 
               header=header, fmt='%.6e')
    
    # Sum spectral density
    np.savetxt(output_dir / 'charge_spectral_density_sum.dat', 
               np.column_stack([results['nu'], results['SD_sum']]),
               header="# Wavenumber (cm-1)  J_total", fmt='%.6e')
    
    # Variance
    np.savetxt(output_dir / 'charge_variance.dat',
               np.column_stack([np.arange(results['n_atoms']), 
                                results['variance_all']]),
               header="# Atom  Variance(e^2)", fmt=['%d', '%.8f'])
    
    print(f"  Data saved to {output_dir}")


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 70)
    print("Charge Difference Spectral Density Analysis")
    print("=" * 70)
    
    # Read data
    print(f"\nReading: {INPUT_FILE}")
    frames, times, charges, atom_types = read_charge_diff_data(INPUT_FILE)
    
    # Interpolate NaN values
    charges = interpolate_nan(charges, times)
    
    # Analyze
    results = analyze_charge_spectral(
        charges, DT, TEMPERATURE,
        segment_ps=SEGMENT_PS,
        corr_ps=CORR_PS
    )
    
    # Plot
    print(f"\n=== Generating plots ===")
    plot_results(results, OUTPUT_DIR, atom_types)
    
    # Save data
    save_results(results, OUTPUT_DIR)
    
    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
