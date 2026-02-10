#!/usr/bin/env python3
"""
Compare spectral densities from different energy difference files.

Plot spectral densities from energy_diff_all.dat, energy_diff_500.dat, 
and energy_diff_1500.dat on the same figure.
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
Ha_to_cm = 219474.63  # cm⁻¹/Ha
Ha_to_eV = 27.211386245988  # eV/Ha
cm_to_eV = 1.23984e-4  # eV/cm⁻¹


# =============================================================================
# File I/O
# =============================================================================
def interpolate_nan(arr: np.ndarray) -> np.ndarray:
    """
    Interpolate NaN values using linear interpolation.
    
    Parameters
    ----------
    arr : np.ndarray
        Input array that may contain NaN values.
    
    Returns
    -------
    np.ndarray
        Array with NaN values replaced by interpolated values.
    """
    arr = arr.copy()
    nan_mask = np.isnan(arr)
    if not np.any(nan_mask):
        return arr
    
    n_nan = np.sum(nan_mask)
    print(f"  Warning: Found {n_nan} NaN values, interpolating...")
    
    # Get indices
    indices = np.arange(len(arr))
    valid_mask = ~nan_mask
    
    # Interpolate
    arr[nan_mask] = np.interp(indices[nan_mask], indices[valid_mask], arr[valid_mask])
    
    return arr


def read_energy_diff(filepath: str) -> dict:
    """
    Read energy_diff file.
    
    Expected format:
    # Energy Difference (Cation - Neutral)
    # Frame  Time(fs)        dE(a.u.)           dE(eV)
        0         0.000        0.2361048333      6.424740
        ...
    
    Returns
    -------
    data : dict
        Dictionary with keys: 'frame', 'time', 'delta_epsilon' (in Ha)
    """
    frames = []
    times = []
    delta_epsilon = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                delta_epsilon.append(float(parts[2]))  # dE(a.u.) column
    
    delta_epsilon_arr = np.array(delta_epsilon)
    
    # Interpolate NaN values if present
    delta_epsilon_arr = interpolate_nan(delta_epsilon_arr)
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'delta_epsilon': delta_epsilon_arr,  # Ha
    }


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
def compute_spectral_density(filepath: str, dt: float, T: float,
                              segment_ps: float = 20.0,
                              corr_ps: float = 10.0) -> dict:
    """
    Compute spectral density from energy difference file.
    """
    # Read data
    data = read_energy_diff(filepath)
    
    # Convert units to cm⁻¹
    delta_epsilon_cm = data['delta_epsilon'] * Ha_to_cm
    
    # Compute fluctuations (subtract mean)
    d_delta_epsilon = delta_epsilon_cm - np.mean(delta_epsilon_cm)
    
    N = len(delta_epsilon_cm)
    total_time_ps = N * dt / 1000
    
    # Determine analysis method based on data length
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    n_segments_expected = int(total_time_ps / segment_ps)
    
    if n_segments_expected >= 2:
        t_corr, C, n_seg = compute_correlation_segmented(
            d_delta_epsilon, dt, segment_length, corr_length)
    else:
        t_corr, C = compute_correlation_fft(d_delta_epsilon, dt)
        max_corr = min(len(t_corr), int(corr_ps * 1000 / dt))
        t_corr = t_corr[:max_corr]
        C = C[:max_corr]
    
    # Spectral density
    nu_max = 2000  # cm⁻¹
    nu_out = np.linspace(1, nu_max, 500)
    SD = correlation_to_spectral_density(t_corr, C, nu_out, T)
    
    # Reorganization energy
    kBT = kB_cm * T
    lambda_cm = C[0] / (2 * kBT)
    lambda_eV = lambda_cm * cm_to_eV
    
    return {
        'nu': nu_out,
        'SD': SD,
        'lambda_cm': lambda_cm,
        'lambda_eV': lambda_eV,
        'N_frames': N,
        'total_time_ps': total_time_ps,
    }


def main():
    """Main routine to compare spectral densities."""
    # Input parameters
    dt = 4.0  # fs (time step between frames)
    T = 300.0  # K (temperature)
    segment_ps = 20.0  # ps
    corr_ps = 10.0  # ps
    
    # Files to compare
    files = {
        'all': 'energy_diff.dat',
        '750-1000': 'energy_diff_750-1000.dat',
    }
    
    # Colors and labels for plotting
    colors = {
        'all': 'black',
        '750-1000': 'cornflowerblue',
    }
    labels = {
        'all': 'Total',
        '750-1000': '750-1000 cm⁻¹ filtered',
    }
    
    print("=== Comparing Spectral Densities ===")
    print(f"Time step: {dt} fs")
    print(f"Temperature: {T} K")
    print()
    
    # Compute spectral density for each file
    results = {}
    for key, filepath in files.items():
        if not Path(filepath).exists():
            print(f"Warning: {filepath} not found, skipping...")
            continue
        
        print(f"Processing {filepath}...")
        results[key] = compute_spectral_density(
            filepath, dt, T, segment_ps=segment_ps, corr_ps=corr_ps
        )
        print(f"  Frames: {results[key]['N_frames']}, "
              f"Time: {results[key]['total_time_ps']:.1f} ps, "
              f"λ = {results[key]['lambda_cm']:.1f} cm⁻¹ = {results[key]['lambda_eV']:.4f} eV")
    
    if not results:
        print("No files found!")
        return
    
    # Plot comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot 750-1000 first (lines), then all on top (dashed)
    for key in ['750-1000']:
        if key not in results:
            continue
        ax.plot(results[key]['nu'], results[key]['SD'], 
                color=colors[key], linewidth=2.0, label=labels[key])
    
    # Plot 'all' on top with dashed line
    if 'all' in results:
        ax.plot(results['all']['nu'], results['all']['SD'], 
                '--', color=colors['all'], linewidth=2.0, label=labels['all'])
    
    ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=15)
    ax.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=15)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlim(0, 2000)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figures
    output_dir = Path('.')
    plt.savefig(output_dir / 'spectral_density_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'spectral_density_comparison.pdf', bbox_inches='tight')
    print()
    print(f"Figure saved: spectral_density_comparison.png")
    print(f"Figure saved: spectral_density_comparison.pdf")
    plt.close()
    
    # Save numerical data
    data_file = output_dir / 'spectral_density_comparison.dat'
    with open(data_file, 'w') as f:
        f.write("# Spectral density comparison\n")
        f.write(f"# Temperature: {T:.1f} K\n")
        for key in results:
            f.write(f"# {labels[key]}: λ = {results[key]['lambda_cm']:.2f} cm⁻¹ = {results[key]['lambda_eV']:.4f} eV\n")
        
        # Header
        header = "# Wavenumber(cm⁻¹)"
        for key in ['all', '750-1000']:
            if key in results:
                header += f"  J_{key}(cm⁻¹)"
        f.write(header + "\n")
        
        # Data
        nu = results[list(results.keys())[0]]['nu']
        for i, nu_val in enumerate(nu):
            line = f"{nu_val:12.2f}"
            for key in ['all', '750-1000']:
                if key in results:
                    line += f"  {results[key]['SD'][i]:16.6e}"
            f.write(line + "\n")
    
    print(f"Data saved: {data_file}")
    print()
    print("=== Done ===")


if __name__ == '__main__':
    main()
