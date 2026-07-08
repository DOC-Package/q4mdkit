#!/usr/bin/env python3
"""
Plot coupling (J) spectral density.

Based on analyze_cdftbci_spectral.py and calc_spectral_density_diff.py
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
cm_to_eV = 1.23984e-4  # eV/cm⁻¹
cm_to_meV = cm_to_eV * 1000
meV_to_cm = 8.0655  # cm⁻¹/meV


# =============================================================================
# File I/O
# =============================================================================
def read_cdftbci_extracted(filepath: str) -> dict:
    """Read cdftbci_extracted.dat file."""
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
                J_val = float(parts[2]) if parts[2].lower() != 'nan' else np.nan
                dE1_val = float(parts[3]) if parts[3].lower() != 'nan' else np.nan
                dE2_val = float(parts[4]) if parts[4].lower() != 'nan' else np.nan
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
            return arr, n_nan
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
        J[i] = two_pi_c * nu / (kB_cm * T) * integral
    
    return J


# =============================================================================
# Main
# =============================================================================
def main():
    """Main analysis routine."""
    # Parameters
    T = 300.0  # K
    dt = 4.0   # fs
    segment_ps = 100.0
    corr_ps_J = 50.0  # Longer correlation for J (slower dynamics)
    nu_max = 2000  # cm⁻¹
    
    script_dir = Path(__file__).parent
    cdftbci_file = script_dir / 'cdftbci_extracted.dat'
    
    print("=" * 70)
    print("Coupling (J) Spectral Density Analysis")
    print("=" * 70)
    
    # Read CDFTB-CI data
    print(f"\nReading {cdftbci_file}...")
    data = read_cdftbci_extracted(cdftbci_file)
    
    N = len(data['frame'])
    total_time_ps = N * dt / 1000
    
    print(f"\n  Total frames: {N}")
    print(f"  Time step: {dt} fs")
    print(f"  Total time: {total_time_ps:.2f} ps")
    
    # Convert J to cm⁻¹
    J_cm = data['J'] * meV_to_cm  # meV -> cm⁻¹
    
    # Compute fluctuations
    delta_J = J_cm - np.mean(J_cm)
    
    print(f"\nStatistics:")
    print(f"  |J|: mean = {np.mean(data['J']):.2f} meV, std = {np.std(data['J']):.2f} meV")
    
    # Compute correlation functions
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps_J * 1000 / dt)
    
    print(f"\nComputing correlation functions...")
    print(f"  Segment length: {segment_ps} ps ({segment_length} frames)")
    print(f"  Correlation length: {corr_ps_J} ps ({corr_length} frames)")
    
    t_corr, C_J, n_seg = compute_correlation_segmented(delta_J, dt, segment_length, corr_length)
    
    print(f"  {n_seg} segments used")
    
    # Compute spectral densities
    nu_out = np.linspace(1, nu_max, 5000)
    
    print(f"\nComputing spectral densities...")
    SD_J = correlation_to_spectral_density(t_corr, C_J, nu_out, T)
    
    # Reorganization energies
    kBT = kB_cm * T
    lambda_J = C_J[0] / (2 * kBT)
    
    print(f"\nReorganization energies:")
    print(f"  λ_J = {lambda_J:.2f} cm⁻¹ = {lambda_J*cm_to_meV:.2f} meV")
    
    # ==========================================================================
    # Plot 1: Coupling spectral density
    # ==========================================================================
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(nu_out, SD_J * cm_to_meV, 'r-', linewidth=2, 
            label=f'$J$ (λ = {lambda_J*cm_to_meV:.2f} meV)')
    
    ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=14)
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=14)
    ax.set_title('Coupling Spectral Density', fontsize=16, fontweight='bold')
    ax.set_xlim(0, nu_max)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    plt.savefig(script_dir / 'coupling_spectral_density.png', dpi=300, bbox_inches='tight')
    plt.savefig(script_dir / 'coupling_spectral_density.pdf', bbox_inches='tight')
    print(f"\nFigures saved: coupling_spectral_density.png, .pdf")
    plt.close()
    
    # ==========================================================================
    # Plot 2: Correlation function and spectral density
    # ==========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Correlation function
    ax = axes[0]
    cm2_to_meV2 = cm_to_meV ** 2
    ax.plot(t_corr / 1000, C_J * cm2_to_meV2, 'r-', linewidth=2, label='$J$')
    ax.set_xlabel('Time (ps)', fontsize=14)
    ax.set_ylabel(r'Autocorrelation Function (meV$^2$)', fontsize=14)
    #ax.set_title('Autocorrelation Function', fontsize=14, fontweight='bold')
    ax.set_xlim(0, corr_ps_J)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    #ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=12)
    ax.text(0.95, 0.95, '(a)', transform=ax.transAxes, fontsize=14, fontweight='bold',
            verticalalignment='top', horizontalalignment='right')
    
    # Spectral density
    ax = axes[1]
    ax.plot(nu_out, SD_J * cm_to_meV, 'r-', linewidth=2, label='$J$')
    ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=16)
    #ax.set_title('Spectral Density', fontsize=16, fontweight='bold')
    ax.set_xlim(0, nu_max)
    ax.set_ylim(bottom=0)
    #ax.legend(fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=16)
    ax.text(0.95, 0.95, '(b)', transform=ax.transAxes, fontsize=16, fontweight='bold',
            verticalalignment='top', horizontalalignment='right')
    
    plt.tight_layout()
    plt.savefig(script_dir / 'coupling_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(script_dir / 'coupling_analysis.pdf', bbox_inches='tight')
    print(f"Figures saved: coupling_analysis.png, .pdf")
    plt.close()
    
    # Save data
    output_file = script_dir / 'coupling_spectral_density.dat'
    with open(output_file, 'w') as f:
        f.write("# Coupling Spectral Density\n")
        f.write(f"# Temperature: {T} K\n")
        f.write(f"# λ_J = {lambda_J*cm_to_meV:.4f} meV\n")
        f.write(f"# |J| mean = {np.mean(data['J']):.2f} meV\n")
        f.write(f"# |J| std = {np.std(data['J']):.2f} meV\n")
        f.write("# nu(cm-1)  SD_J(meV)\n")
        for i in range(len(nu_out)):
            f.write(f"{nu_out[i]:10.2f} {SD_J[i]*cm_to_meV:14.6e}\n")
    print(f"Data saved: {output_file}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
