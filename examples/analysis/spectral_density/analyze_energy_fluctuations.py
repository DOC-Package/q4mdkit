"""
Analyze energy fluctuations from CDFT calculations.

Read energies.dat and compute:
1. Energy fluctuation correlation function C(t)
2. Spectral density J(ω) from Fourier transform of C(t)

The relation between correlation function and spectral density:
    C_cl(t) = (2/π) k_B T ∫_0^∞ dν̃ (J(ν̃)/ν̃) cos(2πc ν̃ t)
    
Inverse relation:
    J(ν̃) = (π ν̃ / k_B T) ∫_{-∞}^∞ dt C_cl(t) exp(-i 2πc ν̃ t)
         = (2π ν̃ / k_B T) ∫_0^∞ dt C_cl(t) cos(2πc ν̃ t)
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

# Hartree to cm⁻¹
Ha_to_cm = 219474.63  # cm⁻¹/Ha

# cm⁻¹ to eV
cm_to_eV = 1.23984e-4  # eV/cm⁻¹


# =============================================================================
# File I/O
# =============================================================================
def read_energies(filepath: str) -> dict:
    """
    Read energies.dat file.
    
    Expected format:
    # Frame  Time[fs]  Energy_frag1[a.u.]  Energy_frag2[a.u.]  [Charge_frag1[e]  Charge_frag2[e]]
        0         0.000      -85.4832131177      -85.4900072463     [+0.999066     +0.999034]
        ...
    
    Returns
    -------
    data : dict
        Dictionary with keys: 'frame', 'time', 'E1', 'E2', 'Q1', 'Q2'
    """
    frames = []
    times = []
    E1_list = []
    E2_list = []
    Q1_list = []
    Q2_list = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                # Handle 'nan' values - convert to np.nan
                E1_val = float(parts[2]) if parts[2].lower() != 'nan' else np.nan
                E2_val = float(parts[3]) if parts[3].lower() != 'nan' else np.nan
                E1_list.append(E1_val)
                E2_list.append(E2_val)
                # Q1, Q2 are optional (may not be present in some files)
                if len(parts) >= 6:
                    Q1_list.append(float(parts[4]))
                    Q2_list.append(float(parts[5]))
                else:
                    Q1_list.append(0.0)
                    Q2_list.append(0.0)
    
    # Convert to numpy arrays
    E1 = np.array(E1_list)
    E2 = np.array(E2_list)
    
    # Interpolate NaN values
    n_nan_E1 = np.sum(np.isnan(E1))
    n_nan_E2 = np.sum(np.isnan(E2))
    
    if n_nan_E1 > 0 or n_nan_E2 > 0:
        print(f"  Warning: Found {n_nan_E1} NaN in E1, {n_nan_E2} NaN in E2. Interpolating...")
        
        # Linear interpolation for NaN values
        def interpolate_nans(arr):
            nan_mask = np.isnan(arr)
            if np.all(nan_mask):
                return arr  # All NaN, can't interpolate
            x = np.arange(len(arr))
            arr[nan_mask] = np.interp(x[nan_mask], x[~nan_mask], arr[~nan_mask])
            return arr
        
        E1 = interpolate_nans(E1)
        E2 = interpolate_nans(E2)
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'E1': E1,
        'E2': E2,
        'Q1': np.array(Q1_list),
        'Q2': np.array(Q2_list)
    }


# =============================================================================
# Correlation function computation
# =============================================================================
def compute_correlation_fft(delta_E: np.ndarray, dt: float) -> tuple:
    """
    Compute autocorrelation function using FFT (Wiener-Khinchin theorem).
    
    C(τ) = <δE(t) δE(t+τ)>
    
    Parameters
    ----------
    delta_E : np.ndarray
        Energy fluctuation time series (mean subtracted)
    dt : float
        Time step (fs)
        
    Returns
    -------
    t_corr : np.ndarray
        Correlation time array (fs)
    C : np.ndarray
        Correlation function (same units as δE²)
    """
    N = len(delta_E)
    
    # Zero-padding for linear (not circular) correlation
    delta_E_padded = np.concatenate([delta_E, np.zeros(N)])
    
    # FFT
    fft_delta_E = fft(delta_E_padded)
    
    # Power spectrum
    power_spectrum = np.abs(fft_delta_E)**2
    
    # Inverse FFT gives autocorrelation
    autocorr_full = np.real(ifft(power_spectrum))[:N]
    
    # Normalization: divide by number of pairs at each lag
    norm = np.arange(N, 0, -1)
    C = autocorr_full / norm
    
    t_corr = np.arange(N) * dt
    
    return t_corr, C


def compute_correlation_direct(delta_E: np.ndarray, dt: float, max_lag: int = None) -> tuple:
    """
    Compute autocorrelation function using direct summation.
    Slower but more transparent.
    
    C(τ_k) = (1/(N-k)) Σ_{i=0}^{N-k-1} δE(t_i) δE(t_{i+k})
    """
    N = len(delta_E)
    if max_lag is None:
        max_lag = N
    
    C = np.zeros(max_lag)
    for k in range(max_lag):
        C[k] = np.mean(delta_E[:N-k] * delta_E[k:])
    
    t_corr = np.arange(max_lag) * dt
    return t_corr, C


def compute_correlation_segmented(delta_E: np.ndarray, dt: float, 
                                   segment_length: int, corr_length: int) -> tuple:
    """
    Compute correlation function with segment averaging.
    
    Split the time series into segments and average the correlation
    functions. This reduces noise for long trajectories.
    
    Parameters
    ----------
    delta_E : np.ndarray
        Energy fluctuation time series
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
    N = len(delta_E)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    C_sum = np.zeros(corr_length)
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        segment = delta_E[start:end]
        
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
                                     use_window: bool = False) -> np.ndarray:
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
    # (factor of 2 for one-sided spectrum)
    S = 2 * np.real(C_fft) * dt
    
    # Spectral density: J(ν) = π ν S(ν) / (2 k_B T)
    # From S(ν) = (2/π) k_B T J(ν) / ν
    with np.errstate(divide='ignore', invalid='ignore'):
        J = np.pi * np.abs(nu) * S / (2 * kB_cm * T)
        J = np.nan_to_num(J, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Return only positive frequencies
    pos_mask = nu >= 0
    return nu[pos_mask], J[pos_mask]


# =============================================================================
# Analysis
# =============================================================================
def analyze_energy_fluctuations(data: dict, dt: float, T: float, 
                                 use_segments: bool = False,
                                 segment_ps: float = None,
                                 corr_ps: float = None) -> dict:
    """
    Analyze energy fluctuations from both fragments.
    
    Parameters
    ----------
    data : dict
        Output from read_energies()
    dt : float
        Time step (fs)
    T : float
        Temperature (K)
    use_segments : bool
        Whether to use segment averaging
    segment_ps : float
        Segment length in ps (for segmented analysis)
    corr_ps : float
        Correlation length in ps
        
    Returns
    -------
    results : dict
        Analysis results for each fragment and the difference
    """
    E1 = data['E1'] * Ha_to_cm  # Convert to cm⁻¹
    E2 = data['E2'] * Ha_to_cm
    
    # Fluctuations (subtract mean)
    dE1 = E1 - np.mean(E1)
    dE2 = E2 - np.mean(E2)
    
    N = len(E1)
    
    results = {
        'dt': dt,
        'T': T,
        'N_frames': N,
        'total_time_ps': N * dt / 1000,
    }
    
    # Statistics
    results['E1_mean'] = np.mean(E1)
    results['E2_mean'] = np.mean(E2)
    results['E1_std'] = np.std(dE1)
    results['E2_std'] = np.std(dE2)
    
    if use_segments and segment_ps is not None:
        segment_length = int(segment_ps * 1000 / dt)
        corr_length = int(corr_ps * 1000 / dt) if corr_ps else segment_length // 2
        
        # Segmented correlation
        t1, C1, n_seg = compute_correlation_segmented(dE1, dt, segment_length, corr_length)
        t2, C2, _ = compute_correlation_segmented(dE2, dt, segment_length, corr_length)
        
        results['n_segments'] = n_seg
        results['segment_length'] = segment_length
    else:
        # Full correlation using direct method
        print(f"  Using direct correlation method (max_lag = {len(dE1)} frames)")
        t1, C1 = compute_correlation_direct(dE1, dt, max_lag=None)
        t2, C2 = compute_correlation_direct(dE2, dt, max_lag=None)
    
    results['t_corr'] = t1
    results['C1'] = C1
    results['C2'] = C2
    
    # Spectral density
    nu_max = 2000  # cm⁻¹
    nu_out = np.linspace(1, nu_max, 500)
    
    results['nu'] = nu_out
    results['J1'] = correlation_to_spectral_density(t1, C1, nu_out, T)
    results['J2'] = correlation_to_spectral_density(t2, C2, nu_out, T)
    
    # Also compute FFT-based spectral density
    nu_fft, J1_fft = spectral_density_fft(t1, C1, T)
    _, J2_fft = spectral_density_fft(t2, C2, T)
    
    results['nu_fft'] = nu_fft
    results['J1_fft'] = J1_fft
    results['J2_fft'] = J2_fft
    
    # Reorganization energy from C(0)
    # Classical fluctuation-dissipation: C(0) = 2 k_B T λ
    # Therefore: λ = C(0) / (2 k_B T)
    kBT = kB_cm * T
    results['lambda1'] = C1[0] / (2 * kBT)
    results['lambda2'] = C2[0] / (2 * kBT)
    
    return results


def plot_results(results: dict, output_dir: str = '.'):
    """Plot correlation functions and spectral densities."""
    output_dir = Path(output_dir)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    t_corr = results['t_corr']
    nu = results['nu']
    
    # Conversion factor: cm⁻² to eV²
    cm2_to_eV2 = cm_to_eV ** 2
    
    # Correlation functions (convert cm⁻² to eV²)
    ax = axes[0]
    ax.plot(t_corr / 1000, results['C1'] * cm2_to_eV2, 'b-', label='Fragment 1', linewidth=1)
    ax.plot(t_corr / 1000, results['C2'] * cm2_to_eV2, 'r-', label='Fragment 2', linewidth=1)
    if 'C_avg' in results:
        ax.plot(t_corr / 1000, results['C_avg'] * cm2_to_eV2, 'k-', label='Average', linewidth=1.5)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('C(t) (eV²)')
    ax.set_title('Site Energy Correlation Functions')
    ax.legend()
    ax.set_xlim(0, min(t_corr[-1]/1000, 20))
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # Spectral density (convert cm⁻¹ to eV)
    ax = axes[1]
    ax.plot(nu, results['J1'] * cm_to_eV, 'b-', label='Fragment 1', linewidth=1, alpha=0.5)
    ax.plot(nu, results['J2'] * cm_to_eV, 'r-', label='Fragment 2', linewidth=1, alpha=0.5)
    if 'J_avg' in results:
        ax.plot(nu, results['J_avg'] * cm_to_eV, 'k-', label='Average', linewidth=1.5)
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('J(ν̃) (eV)')
    ax.set_title('Spectral Density')
    ax.legend()
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'energy_fluctuation_analysis.png', dpi=150)
    plt.close()
    
    print(f"Plots saved to {output_dir}")


def main():
    # =========================================================================
    # User parameters - modify these as needed
    # =========================================================================
    T = 300  # K (temperature)
    dt = 4.0  # fs (time step)
    
    # Time range selection (in ps)
    t_start_ps = 0.0    # Start time (ps) - set to skip initial equilibration
    duration_ps = None  # Duration (ps) - set to None to use all remaining data
    
    # Correlation function parameters
    corr_ps = 30.0      # Correlation length (ps)
    segment_ps = 40.0   # Segment length (ps), should be > corr_ps
    # =========================================================================
    
    # File path
    script_dir = Path(__file__).parent
    
    # Try corrfunc/energies.dat first (longer trajectory), then output/energies.dat
    energy_file = script_dir / 'energies.dat'
    if not energy_file.exists():
        energy_file = script_dir.parent / 'output' / 'energies.dat'
    
    print("=" * 70)
    print("Energy Fluctuation Analysis")
    print("=" * 70)
    
    # Read data
    print(f"\nReading {energy_file}...")
    data = read_energies(energy_file)
    
    N_total = len(data['frame'])
    total_time_full = N_total * dt / 1000  # ps
    
    print(f"  Total frames: {N_total}")
    print(f"  Time step: {dt} fs")
    print(f"  Total time: {total_time_full:.2f} ps")
    
    # Apply time range selection
    start_frame = int(t_start_ps * 1000 / dt)
    if duration_ps is not None:
        end_frame = start_frame + int(duration_ps * 1000 / dt)
    else:
        end_frame = N_total
    
    # Validate range
    if start_frame >= N_total:
        raise ValueError(f"t_start_ps ({t_start_ps} ps) is beyond data range ({total_time_full:.2f} ps)")
    if end_frame > N_total:
        end_frame = N_total
        print(f"  Warning: Requested duration exceeds data, using all remaining data")
    
    # Slice data
    data['frame'] = data['frame'][start_frame:end_frame]
    data['E1'] = data['E1'][start_frame:end_frame]
    data['E2'] = data['E2'][start_frame:end_frame]
    
    N = len(data['frame'])
    total_time = N * dt / 1000  # ps
    
    print(f"\n  Selected time range:")
    print(f"    Start: {t_start_ps:.2f} ps (frame {start_frame})")
    print(f"    End: {(start_frame + N) * dt / 1000:.2f} ps (frame {start_frame + N})")
    print(f"    Duration: {total_time:.2f} ps ({N} frames)")
    
    # Basic statistics
    E1_Ha = data['E1']
    E2_Ha = data['E2']
    print(f"\nEnergy statistics (Hartree):")
    print(f"  E1: mean = {np.mean(E1_Ha):.6f}, std = {np.std(E1_Ha):.6e}")
    print(f"  E2: mean = {np.mean(E2_Ha):.6f}, std = {np.std(E2_Ha):.6e}")
    print(f"  ΔE = E1-E2: mean = {np.mean(E1_Ha - E2_Ha):.6e}, std = {np.std(E1_Ha - E2_Ha):.6e}")
    
    # Convert to cm⁻¹
    E1_cm = E1_Ha * Ha_to_cm
    E2_cm = E2_Ha * Ha_to_cm
    dE1 = E1_cm - np.mean(E1_cm)
    dE2 = E2_cm - np.mean(E2_cm)
    
    print(f"\nEnergy fluctuation statistics (cm⁻¹):")
    print(f"  σ(E1) = {np.std(dE1):.2f} cm⁻¹")
    print(f"  σ(E2) = {np.std(dE2):.2f} cm⁻¹")
    print(f"  σ(ΔE) = {np.std(dE1 - dE2):.2f} cm⁻¹")
    
    # Theoretical relation: σ² = C(0) = 2 k_B T λ
    kBT = kB_cm * T
    print(f"\nTheoretical relation: C(0) = 2 k_B T λ")
    print(f"  k_B T = {kBT:.2f} cm⁻¹")
    print(f"  Estimated λ from variance:")
    print(f"    λ1 = σ²/(2 k_B T) = {np.var(dE1) / (2 * kBT):.2f} cm⁻¹")
    print(f"    λ2 = σ²/(2 k_B T) = {np.var(dE2) / (2 * kBT):.2f} cm⁻¹")
    
    # Full analysis with segment averaging
    # Use segments to get at least 10 ps correlation length
    print("\n" + "-" * 70)
    print("Computing correlation functions and spectral densities...")
    print("-" * 70)
    
    # Check if we have enough data for at least 2 segments
    total_time_ps = N * dt / 1000
    
    # Adjust segment_ps if needed
    if segment_ps < corr_ps:
        # If we can't fit 2 segments with corr_ps, reduce corr_ps
        corr_ps = segment_ps * 0.9
    
    n_segments_expected = int(total_time_ps / segment_ps)
    
    if n_segments_expected >= 2:
        print(f"  Using segmented analysis:")
        print(f"    Segment length: {segment_ps} ps")
        print(f"    Correlation length: {corr_ps} ps")
        print(f"    Expected segments: {n_segments_expected}")
        results = analyze_energy_fluctuations(data, dt, T, 
                                               use_segments=True,
                                               segment_ps=segment_ps,
                                               corr_ps=corr_ps)
        print(f"    Actual segments used: {results.get('n_segments', 'N/A')}")
    else:
        print(f"  Data too short for segmented analysis, using full trajectory")
        results = analyze_energy_fluctuations(data, dt, T, use_segments=False)
    
    print(f"\nCorrelation function at t=0:")
    print(f"  C1(0) = {results['C1'][0]:.2f} cm⁻²")
    print(f"  C2(0) = {results['C2'][0]:.2f} cm⁻²")
    
    print(f"\nReorganization energy from C(0)/(2kBT):")
    print(f"  λ1 = {results['lambda1']:.2f} cm⁻¹")
    print(f"  λ2 = {results['lambda2']:.2f} cm⁻¹")
    print(f"  λ_avg (fragment average) = {(results['lambda1'] + results['lambda2']) / 2:.2f} cm⁻¹")
    
    # Compute average spectral density (fragment 1 + 2) / 2
    results['J_avg'] = (results['J1'] + results['J2']) / 2
    results['C_avg'] = (results['C1'] + results['C2']) / 2
    results['lambda_avg'] = (results['lambda1'] + results['lambda2']) / 2
    
    # Plot
    plot_results(results, script_dir)
    
    # Save numerical results
    output_file = script_dir / 'correlation_results.npz'
    np.savez(output_file,
             t_corr=results['t_corr'],
             C1=results['C1'],
             C2=results['C2'],
             C_avg=results['C_avg'],
             nu=results['nu'],
             J1=results['J1'],
             J2=results['J2'],
             J_avg=results['J_avg'],
             nu_fft=results['nu_fft'],
             J1_fft=results['J1_fft'],
             J2_fft=results['J2_fft'])
    print(f"\nNumerical results saved to {output_file}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
