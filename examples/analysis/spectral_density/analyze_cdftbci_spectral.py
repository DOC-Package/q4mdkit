"""
Analyze J, dE1, dE2 fluctuations from CDFTB-CI calculations.

Read cdftbci_extracted.dat and compute:
1. Autocorrelation functions C(t) for J, dE1, dE2
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
from pathlib import Path

# =============================================================================
# Physical constants
# =============================================================================
kB_cm = 0.695034800  # cm⁻¹/K (Boltzmann constant)
c_cm_fs = 2.99792458e-5  # cm/fs (speed of light)
two_pi_c = 2 * np.pi * c_cm_fs  # rad·cm/fs

# Unit conversions
Ha_to_cm = 219474.63  # cm⁻¹/Ha
cm_to_eV = 1.23984e-4  # eV/cm⁻¹
meV_to_cm = 8.0655  # cm⁻¹/meV


# =============================================================================
# File I/O
# =============================================================================
def read_cdftbci_extracted(filepath: str) -> dict:
    """
    Read cdftbci_extracted.dat file.
    
    Expected format:
    # Frame  Time(fs)   J_lowdin(meV)      dE1(Ha)           dE2(Ha)
         0       0.00       -69.9948       0.2297844942       0.2389978084
        ...
    
    Returns
    -------
    data : dict
        Dictionary with keys: 'frame', 'time', 'J', 'dE1', 'dE2'
        J is converted to absolute values
    """
    frames = []
    times = []
    J_list = []
    dE1_list = []
    dE2_list = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                # Handle 'nan' values - convert to np.nan
                J_val = float(parts[2]) if parts[2].lower() != 'nan' else np.nan
                dE1_val = float(parts[3]) if parts[3].lower() != 'nan' else np.nan
                dE2_val = float(parts[4]) if parts[4].lower() != 'nan' else np.nan
                # J: use absolute value (ignore sign flips)
                J_list.append(np.abs(J_val))
                dE1_list.append(dE1_val)
                dE2_list.append(dE2_val)
    
    # Convert to numpy arrays
    J_arr = np.array(J_list)
    dE1_arr = np.array(dE1_list)
    dE2_arr = np.array(dE2_list)
    
    # Interpolate NaN values
    def interpolate_nans(arr):
        nan_mask = np.isnan(arr)
        n_nan = np.sum(nan_mask)
        if n_nan == 0:
            return arr, 0
        if np.all(nan_mask):
            return arr, n_nan  # All NaN, can't interpolate
        x = np.arange(len(arr))
        arr[nan_mask] = np.interp(x[nan_mask], x[~nan_mask], arr[~nan_mask])
        return arr, n_nan
    
    J_arr, n_nan_J = interpolate_nans(J_arr)
    dE1_arr, n_nan_dE1 = interpolate_nans(dE1_arr)
    dE2_arr, n_nan_dE2 = interpolate_nans(dE2_arr)
    
    if n_nan_J > 0 or n_nan_dE1 > 0 or n_nan_dE2 > 0:
        print(f"  Warning: Interpolated NaN values: J={n_nan_J}, dE1={n_nan_dE1}, dE2={n_nan_dE2}")
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'J': J_arr,      # meV (absolute value)
        'dE1': dE1_arr,  # Ha
        'dE2': dE2_arr,  # Ha
    }


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
def analyze_cdftbci_fluctuations(data: dict, dt: float, T: float, 
                                  segment_ps: float = 40.0,
                                  corr_ps_J: float = 30.0,
                                  corr_ps_E: float = 30.0) -> dict:
    """
    Analyze J, dE1, dE2 fluctuations.
    
    Parameters
    ----------
    data : dict
        Output from read_cdftbci_extracted()
    dt : float
        Time step (fs)
    T : float
        Temperature (K)
    segment_ps : float
        Segment length in ps
    corr_ps_J : float
        Correlation length for J (ps)
    corr_ps_E : float
        Correlation length for dE1, dE2 (ps)
        
    Returns
    -------
    results : dict
        Analysis results
    """
    # Convert units to cm⁻¹
    J_cm = data['J'] * meV_to_cm      # meV -> cm⁻¹
    dE1_cm = data['dE1'] * Ha_to_cm   # Ha -> cm⁻¹
    dE2_cm = data['dE2'] * Ha_to_cm   # Ha -> cm⁻¹
    
    # Compute fluctuations (subtract mean)
    dJ = J_cm - np.mean(J_cm)
    ddE1 = dE1_cm - np.mean(dE1_cm)
    ddE2 = dE2_cm - np.mean(dE2_cm)
    
    N = len(J_cm)
    total_time_ps = N * dt / 1000
    
    results = {
        'dt': dt,
        'T': T,
        'N_frames': N,
        'total_time_ps': total_time_ps,
        # Statistics (cm⁻¹)
        'J_mean_cm': np.mean(J_cm),
        'J_std_cm': np.std(dJ),
        'dE1_mean_cm': np.mean(dE1_cm),
        'dE1_std_cm': np.std(ddE1),
        'dE2_mean_cm': np.mean(dE2_cm),
        'dE2_std_cm': np.std(ddE2),
        # Statistics (meV, Ha)
        'J_mean_meV': np.mean(data['J']),
        'J_std_meV': np.std(data['J'] - np.mean(data['J'])),
        'dE1_mean_Ha': np.mean(data['dE1']),
        'dE1_std_Ha': np.std(data['dE1'] - np.mean(data['dE1'])),
        'dE2_mean_Ha': np.mean(data['dE2']),
        'dE2_std_Ha': np.std(data['dE2'] - np.mean(data['dE2'])),
    }
    
    # Determine analysis method based on data length
    segment_length = int(segment_ps * 1000 / dt)
    corr_length_J = int(corr_ps_J * 1000 / dt)
    corr_length_E = int(corr_ps_E * 1000 / dt)
    n_segments_expected = int(total_time_ps / segment_ps)
    
    if n_segments_expected >= 2:
        print(f"  Using segmented analysis:")
        print(f"    Segment length: {segment_ps} ps ({segment_length} frames)")
        print(f"    Correlation length (J): {corr_ps_J} ps ({corr_length_J} frames)")
        print(f"    Correlation length (E): {corr_ps_E} ps ({corr_length_E} frames)")
        print(f"    Expected segments: {n_segments_expected}")
        
        t_corr_J, C_J, n_seg = compute_correlation_segmented(dJ, dt, segment_length, corr_length_J)
        t_corr_E, C_dE1, _ = compute_correlation_segmented(ddE1, dt, segment_length, corr_length_E)
        _, C_dE2, _ = compute_correlation_segmented(ddE2, dt, segment_length, corr_length_E)
        results['n_segments'] = n_seg
    else:
        print(f"  Using full trajectory correlation...")
        t_corr_J, C_J = compute_correlation_fft(dJ, dt)
        t_corr_E, C_dE1 = compute_correlation_fft(ddE1, dt)
        _, C_dE2 = compute_correlation_fft(ddE2, dt)
        # Truncate to reasonable length
        max_corr_J = min(len(t_corr_J), int(corr_ps_J * 1000 / dt))
        max_corr_E = min(len(t_corr_E), int(corr_ps_E * 1000 / dt))
        t_corr_J = t_corr_J[:max_corr_J]
        C_J = C_J[:max_corr_J]
        t_corr_E = t_corr_E[:max_corr_E]
        C_dE1 = C_dE1[:max_corr_E]
        C_dE2 = C_dE2[:max_corr_E]
    
    results['t_corr_J'] = t_corr_J
    results['t_corr_E'] = t_corr_E
    results['C_J'] = C_J
    results['C_dE1'] = C_dE1
    results['C_dE2'] = C_dE2
    
    # Spectral density
    nu_max = 2000  # cm⁻¹
    nu_out = np.linspace(1, nu_max, 500)
    
    results['nu'] = nu_out
    results['SD_J'] = correlation_to_spectral_density(t_corr_J, C_J, nu_out, T)
    results['SD_dE1'] = correlation_to_spectral_density(t_corr_E, C_dE1, nu_out, T)
    results['SD_dE2'] = correlation_to_spectral_density(t_corr_E, C_dE2, nu_out, T)
    
    # FFT-based spectral density
    nu_fft_J, SD_J_fft = spectral_density_fft(t_corr_J, C_J, T)
    nu_fft_E, SD_dE1_fft = spectral_density_fft(t_corr_E, C_dE1, T)
    _, SD_dE2_fft = spectral_density_fft(t_corr_E, C_dE2, T)
    
    results['nu_fft_J'] = nu_fft_J
    results['nu_fft_E'] = nu_fft_E
    results['SD_J_fft'] = SD_J_fft
    results['SD_dE1_fft'] = SD_dE1_fft
    results['SD_dE2_fft'] = SD_dE2_fft
    
    # Reorganization energy from C(0) / (2 k_B T)
    kBT = kB_cm * T
    results['lambda_J'] = C_J[0] / (2 * kBT)
    results['lambda_dE1'] = C_dE1[0] / (2 * kBT)
    results['lambda_dE2'] = C_dE2[0] / (2 * kBT)
    
    # Store correlation parameters
    results['corr_ps_J'] = corr_ps_J
    results['corr_ps_E'] = corr_ps_E
    
    return results


def plot_results(results: dict, output_dir: str = '.'):
    """Plot correlation functions and spectral densities."""
    output_dir = Path(output_dir)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    t_corr_J = results['t_corr_J']
    t_corr_E = results['t_corr_E']
    corr_ps_J = results['corr_ps_J']
    corr_ps_E = results['corr_ps_E']
    nu = results['nu']
    
    # Conversion factor: cm⁻² to meV²
    cm2_to_meV2 = (cm_to_eV * 1000) ** 2
    
    # Row 1: Correlation functions
    # J correlation
    ax = axes[0, 0]
    ax.plot(t_corr_J / 1000, results['C_J'] * cm2_to_meV2, 'b-', linewidth=1)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('$C_J(t)$ (meV²)')
    ax.set_title('J Autocorrelation')
    ax.set_xlim(0, corr_ps_J)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # dE1 correlation
    ax = axes[0, 1]
    ax.plot(t_corr_E / 1000, results['C_dE1'] * cm2_to_meV2, 'r-', linewidth=1)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('$C_{\\delta E_1}(t)$ (meV²)')
    ax.set_title('δE₁ Autocorrelation')
    ax.set_xlim(0, corr_ps_E)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # dE2 correlation
    ax = axes[0, 2]
    ax.plot(t_corr_E / 1000, results['C_dE2'] * cm2_to_meV2, 'g-', linewidth=1)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('$C_{\\delta E_2}(t)$ (meV²)')
    ax.set_title('δE₂ Autocorrelation')
    ax.set_xlim(0, corr_ps_E)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # Row 2: Spectral densities (meV)
    cm_to_meV = cm_to_eV * 1000
    
    # J spectral density
    ax = axes[1, 0]
    ax.plot(nu, results['SD_J'] * cm_to_meV, 'b-', linewidth=1)
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('$J_J(\\tilde{\\nu})$ (meV)')
    ax.set_title('J Spectral Density')
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)
    
    # dE1 spectral density
    ax = axes[1, 1]
    ax.plot(nu, results['SD_dE1'] * cm_to_meV, 'r-', linewidth=1)
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('$J_{\\delta E_1}(\\tilde{\\nu})$ (meV)')
    ax.set_title('δE₁ Spectral Density')
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)
    
    # dE2 spectral density
    ax = axes[1, 2]
    ax.plot(nu, results['SD_dE2'] * cm_to_meV, 'g-', linewidth=1)
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('$J_{\\delta E_2}(\\tilde{\\nu})$ (meV)')
    ax.set_title('δE₂ Spectral Density')
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cdftbci_spectral_density.png', dpi=150)
    plt.close()
    
    # Additional plot: All spectral densities overlaid
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(nu, results['SD_J'] * cm_to_meV, 'b-', linewidth=1.5, label='J')
    ax.plot(nu, results['SD_dE1'] * cm_to_meV, 'r-', linewidth=1.5, label='δE₁')
    ax.plot(nu, results['SD_dE2'] * cm_to_meV, 'g-', linewidth=1.5, label='δE₂')
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Spectral Density (meV)', fontsize=12)
    ax.set_title('CDFTB-CI Spectral Densities', fontsize=14)
    ax.set_xlim(0, 2000)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'cdftbci_spectral_density_overlay.png', dpi=150)
    plt.close()
    
    print(f"Plots saved to {output_dir}")


def main():
    # =========================================================================
    # User parameters
    # =========================================================================
    T = 300  # K (temperature)
    dt = 4.0  # fs (time step)
    
    # Correlation function parameters
    corr_ps_J = 30.0    # Correlation length for J (ps)
    corr_ps_E = 8.0    # Correlation length for dE1, dE2 (ps)
    segment_ps = 200.0  # Segment length (ps)
    # =========================================================================
    
    script_dir = Path(__file__).parent
    input_file = script_dir / 'cdftbci_extracted.dat'
    
    print("=" * 70)
    print("CDFTB-CI Spectral Density Analysis")
    print("=" * 70)
    
    # Read data
    print(f"\nReading {input_file}...")
    data = read_cdftbci_extracted(input_file)
    
    N = len(data['frame'])
    total_time_ps = N * dt / 1000
    
    print(f"  Total frames: {N}")
    print(f"  Time step: {dt} fs")
    print(f"  Total time: {total_time_ps:.2f} ps")
    
    # Statistics
    print(f"\nStatistics:")
    print(f"  |J| (absolute value):")
    print(f"    mean = {np.mean(data['J']):.2f} meV")
    print(f"    std  = {np.std(data['J']):.2f} meV")
    print(f"  δE₁ (energy fluctuation):")
    print(f"    mean = {np.mean(data['dE1']):.6f} Ha = {np.mean(data['dE1']) * 27.2114:.4f} eV")
    print(f"    std  = {np.std(data['dE1']):.6e} Ha = {np.std(data['dE1']) * 27.2114 * 1000:.2f} meV")
    print(f"  δE₂ (energy fluctuation):")
    print(f"    mean = {np.mean(data['dE2']):.6f} Ha = {np.mean(data['dE2']) * 27.2114:.4f} eV")
    print(f"    std  = {np.std(data['dE2']):.6e} Ha = {np.std(data['dE2']) * 27.2114 * 1000:.2f} meV")
    
    # Analysis
    print("\n" + "-" * 70)
    print("Computing correlation functions and spectral densities...")
    print("-" * 70)
    
    results = analyze_cdftbci_fluctuations(data, dt, T, 
                                            segment_ps=segment_ps, 
                                            corr_ps_J=corr_ps_J,
                                            corr_ps_E=corr_ps_E)
    
    print(f"\nCorrelation function at t=0 (cm⁻²):")
    print(f"  C_J(0) = {results['C_J'][0]:.2f}")
    print(f"  C_δE1(0) = {results['C_dE1'][0]:.2f}")
    print(f"  C_δE2(0) = {results['C_dE2'][0]:.2f}")
    
    print(f"\nReorganization energy λ = C(0)/(2kBT) (cm⁻¹):")
    print(f"  λ_J   = {results['lambda_J']:.2f} cm⁻¹ = {results['lambda_J'] * cm_to_eV * 1000:.2f} meV")
    print(f"  λ_δE1 = {results['lambda_dE1']:.2f} cm⁻¹ = {results['lambda_dE1'] * cm_to_eV * 1000:.2f} meV")
    print(f"  λ_δE2 = {results['lambda_dE2']:.2f} cm⁻¹ = {results['lambda_dE2'] * cm_to_eV * 1000:.2f} meV")
    
    # Plot
    plot_results(results, script_dir)
    
    # Save numerical results
    output_file = script_dir / 'cdftbci_spectral_results.npz'
    np.savez(output_file,
             t_corr_J=results['t_corr_J'],
             t_corr_E=results['t_corr_E'],
             C_J=results['C_J'],
             C_dE1=results['C_dE1'],
             C_dE2=results['C_dE2'],
             nu=results['nu'],
             SD_J=results['SD_J'],
             SD_dE1=results['SD_dE1'],
             SD_dE2=results['SD_dE2'],
             nu_fft_J=results['nu_fft_J'],
             nu_fft_E=results['nu_fft_E'],
             SD_J_fft=results['SD_J_fft'],
             SD_dE1_fft=results['SD_dE1_fft'],
             SD_dE2_fft=results['SD_dE2_fft'],
             corr_ps_J=results['corr_ps_J'],
             corr_ps_E=results['corr_ps_E'])
    print(f"\nNumerical results saved to {output_file}")
    
    # Save spectral density in text format
    output_txt = script_dir / 'cdftbci_spectral_density.dat'
    with open(output_txt, 'w') as f:
        f.write("# CDFTB-CI Spectral Densities\n")
        f.write("# nu(cm-1)      SD_J(cm-1)      SD_dE1(cm-1)    SD_dE2(cm-1)\n")
        for i in range(len(results['nu'])):
            f.write(f"{results['nu'][i]:12.4f}  {results['SD_J'][i]:14.6f}  {results['SD_dE1'][i]:14.6f}  {results['SD_dE2'][i]:14.6f}\n")
    print(f"Text results saved to {output_txt}")
    
    # Save summary statistics (mean values and reorganization energies)
    summary_file = script_dir / 'cdftbci_summary.dat'
    dE_mean = 0.5 * (np.mean(data['dE1']) + np.mean(data['dE2']))  # Average site energy
    dE_diff_mean = np.mean(data['dE1'] - data['dE2'])  # Mean energy gap
    with open(summary_file, 'w') as f:
        f.write("# CDFTB-CI Analysis Summary\n")
        f.write("# ============================================================\n")
        f.write(f"# Total frames: {N}\n")
        f.write(f"# Time step: {dt} fs\n")
        f.write(f"# Total time: {total_time_ps:.2f} ps\n")
        f.write(f"# Temperature: {T} K\n")
        f.write("# ============================================================\n")
        f.write("#\n")
        f.write("# Mean values\n")
        f.write("# ------------------------------------------------------------\n")
        f.write(f"J_mean_meV          {np.mean(data['J']):12.4f}   # |J| (meV)\n")
        f.write(f"J_std_meV           {np.std(data['J']):12.4f}   # σ_J (meV)\n")
        f.write(f"dE1_mean_Ha         {np.mean(data['dE1']):12.8f}   # <δE₁> (Ha)\n")
        f.write(f"dE1_mean_eV         {np.mean(data['dE1']) * 27.2114:12.6f}   # <δE₁> (eV)\n")
        f.write(f"dE2_mean_Ha         {np.mean(data['dE2']):12.8f}   # <δE₂> (Ha)\n")
        f.write(f"dE2_mean_eV         {np.mean(data['dE2']) * 27.2114:12.6f}   # <δE₂> (eV)\n")
        f.write(f"dE_avg_mean_eV      {dE_mean * 27.2114:12.6f}   # (<δE₁>+<δE₂>)/2 (eV)\n")
        f.write(f"dE_diff_mean_meV    {dE_diff_mean * 27.2114 * 1000:12.4f}   # <δE₁-δE₂> (meV)\n")
        f.write("#\n")
        f.write("# Standard deviations (fluctuations)\n")
        f.write("# ------------------------------------------------------------\n")
        f.write(f"dE1_std_meV         {np.std(data['dE1']) * 27.2114 * 1000:12.4f}   # σ_δE₁ (meV)\n")
        f.write(f"dE2_std_meV         {np.std(data['dE2']) * 27.2114 * 1000:12.4f}   # σ_δE₂ (meV)\n")
        f.write("#\n")
        f.write("# Reorganization energies λ = C(0)/(2kBT)\n")
        f.write("# ------------------------------------------------------------\n")
        f.write(f"lambda_J_cm         {results['lambda_J']:12.4f}   # λ_J (cm⁻¹)\n")
        f.write(f"lambda_J_meV        {results['lambda_J'] * cm_to_eV * 1000:12.4f}   # λ_J (meV)\n")
        f.write(f"lambda_dE1_cm       {results['lambda_dE1']:12.4f}   # λ_δE₁ (cm⁻¹)\n")
        f.write(f"lambda_dE1_meV      {results['lambda_dE1'] * cm_to_eV * 1000:12.4f}   # λ_δE₁ (meV)\n")
        f.write(f"lambda_dE2_cm       {results['lambda_dE2']:12.4f}   # λ_δE₂ (cm⁻¹)\n")
        f.write(f"lambda_dE2_meV      {results['lambda_dE2'] * cm_to_eV * 1000:12.4f}   # λ_δE₂ (meV)\n")
    print(f"Summary saved to {summary_file}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
