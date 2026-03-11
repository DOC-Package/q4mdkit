#!/usr/bin/env python3
"""
Analyze energy difference (Δε) fluctuations and compute spectral density.

Read energy_diff_all.dat and compute:
1. Autocorrelation function C(t) for Δε
2. Spectral density J(ω) from Fourier transform of C(t)

The relation between correlation function and spectral density:
    C_cl(t) = (2/π) k_B T ∫_0^∞ dν̃ (J(ν̃)/ν̃) cos(2πc ν̃ t)
    
Inverse relation:
    J(ν̃) = (2π ν̃ / k_B T) ∫_0^∞ dt C_cl(t) cos(2πc ν̃ t)
"""

import numpy as np
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
cm_to_eV = 1.23984e-4  # eV/cm⁻¹
meV_to_cm = 8.0655  # cm⁻¹/meV


# =============================================================================
# File I/O
# =============================================================================
def read_energy_diff(filepath: str) -> dict:
    """
    Read energy_diff_all.dat file.
    
    Expected format:
    # Energy differences: energies*.dat - qm_energy_sampled.dat
    # Frame  Time(fs)        Energy_orig(a.u.)  Energy_qm(a.u.)    Difference(a.u.)
         0         0.000      -42.6202870676      -42.8561405051        0.2358534375
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
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'delta_epsilon': np.array(delta_epsilon),  # Ha
    }


def interpolate_nan(data: dict) -> dict:
    """
    Interpolate NaN values in the delta_epsilon array using linear interpolation.
    
    Parameters
    ----------
    data : dict
        Dictionary with keys: 'frame', 'time', 'delta_epsilon'
        
    Returns
    -------
    data : dict
        Dictionary with NaN values interpolated
    """
    delta_epsilon = data['delta_epsilon']
    time = data['time']
    
    # Find NaN indices
    nan_mask = np.isnan(delta_epsilon)
    n_nan = np.sum(nan_mask)
    
    if n_nan == 0:
        print("  No NaN values found in data.")
        return data
    
    print(f"  Found {n_nan} NaN values ({100*n_nan/len(delta_epsilon):.2f}%), interpolating...")
    
    # Get valid indices
    valid_mask = ~nan_mask
    
    if np.sum(valid_mask) < 2:
        raise ValueError("Not enough valid data points for interpolation")
    
    # Create interpolation function from valid data
    f_interp = interp1d(time[valid_mask], delta_epsilon[valid_mask], 
                        kind='linear', fill_value='extrapolate')
    
    # Interpolate NaN values
    delta_epsilon_interp = delta_epsilon.copy()
    delta_epsilon_interp[nan_mask] = f_interp(time[nan_mask])
    
    data['delta_epsilon'] = delta_epsilon_interp
    
    return data


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
    
    Split the time series into segments and average the correlation
    functions. This reduces noise for long trajectories.
    
    Parameters
    ----------
    delta : np.ndarray
        Fluctuation time series
    dt : float
        Time step (fs)
    segment_length : int
        Number of points per segment
    corr_length : int
        Number of correlation points to compute
        
    Returns
    -------
    t_corr : np.ndarray
        Correlation time array
    C_avg : np.ndarray
        Averaged correlation function
    n_segments : int
        Number of segments used
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
        # Use Hanning window (raised cosine)
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


def spectral_density_fft(t_corr: np.ndarray, C: np.ndarray, T: float) -> tuple:
    """
    Compute spectral density using FFT.
    
    Returns wavenumber and spectral density arrays.
    """
    dt = t_corr[1] - t_corr[0]
    N = len(C)
    
    # FFT of correlation function
    C_fft = fft(C)
    
    # Frequency array
    freq = fftfreq(N, d=dt)  # 1/fs
    
    # Convert to wavenumber
    nu = freq / c_cm_fs  # cm⁻¹
    
    # Power spectral density: S(ν) = 2 Re[C̃(ν)] dt
    S = 2 * np.real(C_fft) * dt
    
    # Spectral density: J(ν) = π ν S(ν) / (2 k_B T)
    with np.errstate(divide='ignore', invalid='ignore'):
        J = np.pi * np.abs(nu) * S / (2 * kB_cm * T)
        J = np.nan_to_num(J, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Return only positive frequencies
    pos_mask = nu >= 0
    return nu[pos_mask], J[pos_mask]


# =============================================================================
# Analysis
# =============================================================================
def analyze_energy_diff_fluctuations(data: dict, dt: float, T: float, 
                                      segment_ps: float = 40.0,
                                      corr_ps: float = 30.0) -> dict:
    """
    Analyze energy difference (Δε) fluctuations.
    
    Parameters
    ----------
    data : dict
        Output from read_energy_diff()
    dt : float
        Time step (fs)
    T : float
        Temperature (K)
    segment_ps : float
        Segment length in ps
    corr_ps : float
        Correlation length (ps)
        
    Returns
    -------
    results : dict
        Analysis results
    """
    # Convert units to cm⁻¹
    delta_epsilon_cm = data['delta_epsilon'] * Ha_to_cm   # Ha -> cm⁻¹
    
    # Compute fluctuations (subtract mean)
    d_delta_epsilon = delta_epsilon_cm - np.mean(delta_epsilon_cm)
    
    N = len(delta_epsilon_cm)
    total_time_ps = N * dt / 1000
    
    results = {
        'dt': dt,
        'T': T,
        'N_frames': N,
        'total_time_ps': total_time_ps,
        # Statistics (cm⁻¹)
        'delta_epsilon_mean_cm': np.mean(delta_epsilon_cm),
        'delta_epsilon_std_cm': np.std(d_delta_epsilon),
        # Statistics (Ha, eV)
        'delta_epsilon_mean_Ha': np.mean(data['delta_epsilon']),
        'delta_epsilon_std_Ha': np.std(data['delta_epsilon'] - np.mean(data['delta_epsilon'])),
        'delta_epsilon_mean_eV': np.mean(data['delta_epsilon']) * Ha_to_eV,
        'delta_epsilon_std_eV': np.std(data['delta_epsilon'] - np.mean(data['delta_epsilon'])) * Ha_to_eV,
    }
    
    print(f"\n=== Energy Difference Statistics ===")
    print(f"  Mean Δε: {results['delta_epsilon_mean_Ha']:.6f} Ha = {results['delta_epsilon_mean_eV']:.4f} eV")
    print(f"  Std Δε:  {results['delta_epsilon_std_Ha']:.6f} Ha = {results['delta_epsilon_std_eV']:.4f} eV")
    print(f"  Total trajectory: {total_time_ps:.1f} ps ({N} frames)")
    
    # Determine analysis method based on data length
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    n_segments_expected = int(total_time_ps / segment_ps)
    
    if n_segments_expected >= 2:
        print(f"\n=== Using segmented analysis ===")
        print(f"  Segment length: {segment_ps} ps ({segment_length} frames)")
        print(f"  Correlation length: {corr_ps} ps ({corr_length} frames)")
        print(f"  Expected segments: {n_segments_expected}")
        
        t_corr, C, n_seg = compute_correlation_segmented(d_delta_epsilon, dt, segment_length, corr_length)
        results['n_segments'] = n_seg
        print(f"  Actual segments used: {n_seg}")
    else:
        print(f"\n=== Using full trajectory correlation ===")
        t_corr, C = compute_correlation_fft(d_delta_epsilon, dt)
        # Truncate to reasonable length
        max_corr = min(len(t_corr), int(corr_ps * 1000 / dt))
        t_corr = t_corr[:max_corr]
        C = C[:max_corr]
    
    results['t_corr'] = t_corr
    results['C'] = C
    
    # Spectral density
    nu_max = 2000  # cm⁻¹
    nu_out = np.linspace(1, nu_max, 500)
    
    results['nu'] = nu_out
    results['SD'] = correlation_to_spectral_density(t_corr, C, nu_out, T)
    
    # Reorganization energy from C(0) / (2 k_B T)
    kBT = kB_cm * T
    results['lambda_cm'] = C[0] / (2 * kBT)
    results['lambda_eV'] = results['lambda_cm'] * cm_to_eV
    
    print(f"\n=== Reorganization Energy ===")
    print(f"  λ = C(0) / (2 k_B T) = {results['lambda_cm']:.2f} cm⁻¹ = {results['lambda_eV']:.4f} eV")
    
    # Store correlation parameters
    results['corr_ps'] = corr_ps
    
    return results


def plot_results(results: dict, output_dir: str = '.'):
    """Plot correlation functions and spectral densities."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    t_corr = results['t_corr']
    corr_ps = results['corr_ps']
    nu = results['nu']
    
    # Conversion factor: cm⁻² to eV²
    cm2_to_eV2 = cm_to_eV ** 2
    
    # Left: Correlation function in eV²
    ax = axes[0]
    ax.plot(t_corr / 1000, results['C'] * cm2_to_eV2, 'b-', linewidth=1.5)
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel(r'$C_{\Delta\epsilon}(t)$ (eV²)', fontsize=12)
    ax.set_title('Energy Difference Autocorrelation', fontsize=13, fontweight='bold')
    ax.set_xlim(0, corr_ps)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # Right: Spectral density
    ax = axes[1]
    ax.plot(nu, results['SD'], 'r-', linewidth=1.5)
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel(r'$J(\tilde{\nu})$ (cm⁻¹)', fontsize=12)
    ax.set_title('Spectral Density', fontsize=13, fontweight='bold')
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_dir / 'energy_diff_spectral.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n=== Output ===")
    print(f"  Figure saved: {output_file}")
    plt.close()
    
    # Plot spectral density only
    fig_sd, ax_sd = plt.subplots(figsize=(7, 5))
    ax_sd.plot(nu, results['SD'], 'r-', linewidth=1.5)
    ax_sd.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=15)
    ax_sd.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=15)
    ax_sd.tick_params(axis='both', labelsize=14)
    #ax_sd.set_title('Spectral Density', fontsize=13, fontweight='bold')
    ax_sd.set_xlim(0, 2000)

    ax_sd.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save spectral density figure in PNG and PDF
    sd_png = output_dir / 'spectral_density.png'
    sd_pdf = output_dir / 'spectral_density.pdf'
    plt.savefig(sd_png, dpi=300, bbox_inches='tight')
    plt.savefig(sd_pdf, bbox_inches='tight')
    print(f"  Spectral density figure saved: {sd_png}")
    print(f"  Spectral density figure saved: {sd_pdf}")
    plt.close()
    
    # Plot J(ω)/ω (spectral density divided by omega)
    fig_sd_over_omega, ax_sd_over_omega = plt.subplots(figsize=(7, 5))
    # Avoid division by zero at ω=0
    with np.errstate(divide='ignore', invalid='ignore'):
        SD_over_omega = np.where(nu > 0, results['SD'] / nu, 0.0)
    ax_sd_over_omega.plot(nu, SD_over_omega, 'b-', linewidth=1.5)
    ax_sd_over_omega.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=15)
    ax_sd_over_omega.set_ylabel(r'$J(\omega)/\omega$ (dimensionless)', fontsize=15)
    ax_sd_over_omega.tick_params(axis='both', labelsize=14)
    ax_sd_over_omega.set_xlim(0, 2000)
    ax_sd_over_omega.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save J(ω)/ω figure in PNG and PDF
    sd_over_omega_png = output_dir / 'spectral_density_over_omega.png'
    sd_over_omega_pdf = output_dir / 'spectral_density_over_omega.pdf'
    plt.savefig(sd_over_omega_png, dpi=300, bbox_inches='tight')
    plt.savefig(sd_over_omega_pdf, bbox_inches='tight')
    print(f"  J(ω)/ω figure saved: {sd_over_omega_png}")
    print(f"  J(ω)/ω figure saved: {sd_over_omega_pdf}")
    plt.close()
    
    # Save numerical data
    data_file = output_dir / 'correlation_function.dat'
    with open(data_file, 'w') as f:
        f.write("# Autocorrelation function of energy difference\n")
        f.write(f"# Temperature: {results['T']:.1f} K\n")
        f.write(f"# Reorganization energy: {results['lambda_cm']:.2f} cm⁻¹ = {results['lambda_eV']:.4f} eV\n")
        f.write("# Time(ps)    C(t)(cm⁻²)       C(t)(eV²)\n")
        for t, c in zip(t_corr, results['C']):
            f.write(f"{t/1000:10.4f}  {c:16.6e}  {c*cm2_to_eV2:16.6e}\n")
    print(f"  Correlation data saved: {data_file}")
    
    data_file = output_dir / 'spectral_density.dat'
    with open(data_file, 'w') as f:
        f.write("# Spectral density of energy difference\n")
        f.write(f"# Temperature: {results['T']:.1f} K\n")
        f.write(f"# Reorganization energy: {results['lambda_cm']:.2f} cm⁻¹ = {results['lambda_eV']:.4f} eV\n")
        f.write("# Wavenumber(cm⁻¹)    J(cm⁻¹)\n")
        for nu_val, j_val in zip(results['nu'], results['SD']):
            f.write(f"{nu_val:12.2f}  {j_val:16.6e}\n")
    print(f"  Spectral density saved: {data_file}")


def main():
    """Main analysis routine."""
    # Input parameters
    input_file = 'energy_diff_wopc.dat'
    dt = 4.0  # fs (time step between frames)
    T = 300.0  # K (temperature)
    
    # Analysis parameters
    segment_ps = 20.0  # ps (segment length for averaging)
    corr_ps = 10.0  # ps (correlation time to compute)
    
    print(f"=== Energy Difference Spectral Analysis ===")
    print(f"Input file: {input_file}")
    print(f"Time step: {dt} fs")
    print(f"Temperature: {T} K")
    
    # Read data
    data = read_energy_diff(input_file)
    print(f"\nLoaded {len(data['frame'])} frames")
    
    # Interpolate NaN values if present
    data = interpolate_nan(data)
    
    # Analyze
    results = analyze_energy_diff_fluctuations(
        data, dt, T, 
        segment_ps=segment_ps,
        corr_ps=corr_ps
    )
    
    # Plot
    plot_results(results, output_dir='.')
    
    print("\n=== Analysis complete ===")


if __name__ == '__main__':
    main()
