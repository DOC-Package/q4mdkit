#!/usr/bin/env python3
"""
Plot local energy spectral densities on the same figure.

Based on analyze_cdftbci_spectral.py and calc_spectral_density_diff.py

Plots:
- Average local energy (E_1 + E_2) / 2 from CDFTB-CI
- Total energy difference from energy_diff.dat
- Electrostatic contribution from energy_diff_elstat.dat
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


# =============================================================================
# File I/O
# =============================================================================
def read_energy_diff(filepath: str) -> dict:
    """Read energy_diff*.dat file."""
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
                # Handle NaN values
                val_str = parts[2].lower()
                if val_str == 'nan' or val_str == '-nan':
                    delta_epsilon.append(np.nan)
                else:
                    delta_epsilon.append(float(parts[2]))  # dE(a.u.) column
    
    # Convert to numpy array
    delta_epsilon_arr = np.array(delta_epsilon)
    
    # Interpolate NaN values
    nan_mask = np.isnan(delta_epsilon_arr)
    n_nan = np.sum(nan_mask)
    if n_nan > 0:
        print(f"  Warning: Interpolating {n_nan} NaN values in {filepath}")
        x = np.arange(len(delta_epsilon_arr))
        valid_mask = ~nan_mask
        if np.sum(valid_mask) >= 2:
            delta_epsilon_arr[nan_mask] = np.interp(x[nan_mask], x[valid_mask], delta_epsilon_arr[valid_mask])
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'delta_epsilon': delta_epsilon_arr,  # Ha
    }


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
        'J': J_arr,
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
                                     use_window: bool = False) -> np.ndarray:
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
# Main
# =============================================================================
def main():
    """Main analysis routine."""
    # Parameters
    T = 300.0  # K
    dt = 4.0   # fs
    segment_ps = 100.0
    corr_ps = 8.0
    nu_max = 2000  # cm⁻¹
    
    script_dir = Path(__file__).parent
    cdftbci_file = script_dir / 'cdftbci_extracted.dat'
    energy_diff_file = script_dir / 'energy_diff.dat'
    energy_diff_elstat_file = script_dir / 'energy_diff_elstat.dat'
    
    print("=" * 70)
    print("Local Energy Spectral Density Analysis")
    print("=" * 70)
    
    # Read CDFTB-CI data
    print(f"\nReading {cdftbci_file}...")
    data = read_cdftbci_extracted(cdftbci_file)
    
    # Read energy_diff data
    print(f"Reading {energy_diff_file}...")
    data_total = read_energy_diff(energy_diff_file)
    
    print(f"Reading {energy_diff_elstat_file}...")
    data_elstat = read_energy_diff(energy_diff_elstat_file)
    
    # Match data lengths
    N = min(len(data['frame']), len(data_total['frame']), len(data_elstat['frame']))
    total_time_ps = N * dt / 1000
    
    print(f"\n  Total frames: {N}")
    print(f"  Time step: {dt} fs")
    print(f"  Total time: {total_time_ps:.2f} ps")
    
    # Convert to cm⁻¹
    dE1_cm = data['dE1'][:N] * Ha_to_cm
    dE2_cm = data['dE2'][:N] * Ha_to_cm
    
    # Energy differences from MM calculations
    eps_total_cm = data_total['delta_epsilon'][:N] * Ha_to_cm
    eps_elstat_cm = data_elstat['delta_epsilon'][:N] * Ha_to_cm
    
    # Compute fluctuations
    delta_dE1 = dE1_cm - np.mean(dE1_cm)
    delta_dE2 = dE2_cm - np.mean(dE2_cm)
    delta_total = eps_total_cm - np.mean(eps_total_cm)
    delta_elstat = eps_elstat_cm - np.mean(eps_elstat_cm)
    
    print(f"\nStatistics:")
    print(f"  E₁:      mean = {np.mean(dE1_cm)*cm_to_eV:.4f} eV, std = {np.std(delta_dE1)*cm_to_meV:.2f} meV")
    print(f"  E₂:      mean = {np.mean(dE2_cm)*cm_to_eV:.4f} eV, std = {np.std(delta_dE2)*cm_to_meV:.2f} meV")
    print(f"  Total:   mean = {np.mean(eps_total_cm)*cm_to_eV:.4f} eV, std = {np.std(delta_total)*cm_to_meV:.2f} meV")
    print(f"  Elstat:  mean = {np.mean(eps_elstat_cm)*cm_to_eV:.4f} eV, std = {np.std(delta_elstat)*cm_to_meV:.2f} meV")
    
    # Compute correlation functions
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    
    print(f"\nComputing correlation functions...")
    print(f"  Segment length: {segment_ps} ps ({segment_length} frames)")
    print(f"  Correlation length: {corr_ps} ps ({corr_length} frames)")
    
    t_corr, C_dE1, n_seg = compute_correlation_segmented(delta_dE1, dt, segment_length, corr_length)
    _, C_dE2, _ = compute_correlation_segmented(delta_dE2, dt, segment_length, corr_length)
    _, C_total, _ = compute_correlation_segmented(delta_total, dt, segment_length, corr_length)
    _, C_elstat, _ = compute_correlation_segmented(delta_elstat, dt, segment_length, corr_length)
    
    print(f"  {n_seg} segments used")
    
    # Compute spectral densities
    nu_out = np.linspace(1, nu_max, 500)
    
    print(f"\nComputing spectral densities...")
    SD_dE1 = correlation_to_spectral_density(t_corr, C_dE1, nu_out, T)
    SD_dE2 = correlation_to_spectral_density(t_corr, C_dE2, nu_out, T)
    SD_total = correlation_to_spectral_density(t_corr, C_total, nu_out, T)
    SD_elstat = correlation_to_spectral_density(t_corr, C_elstat, nu_out, T)
    
    # Reorganization energies
    kBT = kB_cm * T
    lambda_dE1 = C_dE1[0] / (2 * kBT)
    lambda_dE2 = C_dE2[0] / (2 * kBT)
    lambda_total = C_total[0] / (2 * kBT)
    lambda_elstat = C_elstat[0] / (2 * kBT)
    
    print(f"\nReorganization energies:")
    print(f"  λ_E₁     = {lambda_dE1:.2f} cm⁻¹ = {lambda_dE1*cm_to_meV:.2f} meV")
    print(f"  λ_E₂     = {lambda_dE2:.2f} cm⁻¹ = {lambda_dE2*cm_to_meV:.2f} meV")
    print(f"  λ_total  = {lambda_total:.2f} cm⁻¹ = {lambda_total*cm_to_meV:.2f} meV")
    print(f"  λ_elstat = {lambda_elstat:.2f} cm⁻¹ = {lambda_elstat*cm_to_meV:.2f} meV")
    
    # ==========================================================================
    # Plot: All spectral densities on the same figure
    # ==========================================================================
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot in meV with different line styles for visibility
    ax.plot(nu_out, SD_dE1 * cm_to_meV, 'r-', linewidth=2, 
            label=f'$E_I$ (λ = {lambda_dE1*cm_to_meV:.2f} meV)')
    ax.plot(nu_out, SD_dE2 * cm_to_meV, 'b--', linewidth=2, 
            label=f'$E_J$ (λ = {lambda_dE2*cm_to_meV:.2f} meV)')
    ax.plot(nu_out, SD_total * cm_to_meV, 'k:', linewidth=2.5, 
            label=f'Single total (λ = {lambda_total*cm_to_meV:.2f} meV)')
    ax.plot(nu_out, SD_elstat * cm_to_meV, 'g-.', linewidth=2, 
            label=f'Single electrostatic (λ = {lambda_elstat*cm_to_meV:.2f} meV)')
    
    ax.set_xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)', fontsize=14)
    ax.set_ylabel(r'$J(\tilde{\nu})$ (meV)', fontsize=14)
    ax.set_title('Local Energy Spectral Densities', fontsize=16, fontweight='bold')
    ax.set_xlim(0, nu_max)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    plt.savefig(script_dir / 'local_energy_spectral_densities.png', dpi=300, bbox_inches='tight')
    plt.savefig(script_dir / 'local_energy_spectral_densities.pdf', bbox_inches='tight')
    print(f"\nFigures saved: local_energy_spectral_densities.png, .pdf")
    plt.close()
    
    # ==========================================================================
    # Additional plot: With correlation functions
    # ==========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Correlation functions
    ax = axes[0]
    cm2_to_meV2 = cm_to_meV ** 2
    ax.plot(t_corr / 1000, C_dE1 * cm2_to_meV2, 'r-', linewidth=2, label='$E_I$')
    ax.plot(t_corr / 1000, C_dE2 * cm2_to_meV2, 'b--', linewidth=2, label='$E_J$')
    ax.plot(t_corr / 1000, C_total * cm2_to_meV2, 'k:', linewidth=2.5, label='Single total')
    ax.plot(t_corr / 1000, C_elstat * cm2_to_meV2, 'g-.', linewidth=2, label='Single electrostatic')
    ax.set_xlabel('Time (ps)', fontsize=16)
    ax.set_ylabel(r'Autocorrelation Function (meV$^2$)', fontsize=16)
    #ax.set_title('Autocorrelation Functions', fontsize=16, fontweight='bold')
    ax.set_xlim(0, corr_ps)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.legend(fontsize=16, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=16)
    ax.text(0.95, 0.95, '(a)', transform=ax.transAxes, fontsize=16, fontweight='bold',
            verticalalignment='top', horizontalalignment='right')
    
    # Spectral densities
    ax = axes[1]
    ax.plot(nu_out, SD_dE1 * cm_to_meV, 'r-', linewidth=2, label='$E_I$')
    ax.plot(nu_out, SD_dE2 * cm_to_meV, 'b--', linewidth=2, label='$E_J$')
    ax.plot(nu_out, SD_total * cm_to_meV, 'k:', linewidth=2.5, label='Single total')
    ax.plot(nu_out, SD_elstat * cm_to_meV, 'g-.', linewidth=2, label='Single electrostatic')
    ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=16)
    #ax.set_title('Spectral Densities', fontsize=14, fontweight='bold')
    ax.set_xlim(0, nu_max)
    ax.set_ylim(bottom=0)
    # Remove legend, add inset instead
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=16)
    ax.text(0.95, 0.95, '(b)', transform=ax.transAxes, fontsize=16, fontweight='bold',
            verticalalignment='top', horizontalalignment='right')
    
    # Inset for low frequency region (upper left)
    ax_inset = ax.inset_axes([0.09, 0.40, 0.55, 0.52])  # [x, y, width, height] in axes coordinates
    nu_max_inset = 400  # cm⁻¹
    ax_inset.plot(nu_out, SD_dE1 * cm_to_meV, 'r-', linewidth=1.5)
    ax_inset.plot(nu_out, SD_dE2 * cm_to_meV, 'b--', linewidth=1.5)
    ax_inset.plot(nu_out, SD_total * cm_to_meV, 'k:', linewidth=2)
    ax_inset.plot(nu_out, SD_elstat * cm_to_meV, 'g-.', linewidth=1.5)
    ax_inset.set_xlim(0, nu_max_inset)
    ax_inset.set_ylim(0, 180)
    ax_inset.tick_params(axis='both', labelsize=12)
    ax_inset.grid(True, alpha=0.3)
    #ax_inset.set_xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)', fontsize=12)
    #ax_inset.set_ylabel(r'$J$ (meV)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(script_dir / 'local_energy_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(script_dir / 'local_energy_analysis.pdf', bbox_inches='tight')
    print(f"Figures saved: local_energy_analysis.png, .pdf")
    plt.close()
    
    # Save data
    output_file = script_dir / 'local_energy_spectral_densities.dat'
    with open(output_file, 'w') as f:
        f.write("# Local Energy Spectral Densities\n")
        f.write(f"# Temperature: {T} K\n")
        f.write(f"# λ_EI = {lambda_dE1*cm_to_meV:.4f} meV\n")
        f.write(f"# λ_EJ = {lambda_dE2*cm_to_meV:.4f} meV\n")
        f.write(f"# λ_total = {lambda_total*cm_to_meV:.4f} meV\n")
        f.write(f"# λ_elstat = {lambda_elstat*cm_to_meV:.4f} meV\n")
        f.write("# nu(cm-1)  SD_EI(meV)  SD_EJ(meV)  SD_total(meV)  SD_elstat(meV)\n")
        for i in range(len(nu_out)):
            f.write(f"{nu_out[i]:10.2f} {SD_dE1[i]*cm_to_meV:14.6e} {SD_dE2[i]*cm_to_meV:14.6e} {SD_total[i]*cm_to_meV:14.6e} {SD_elstat[i]*cm_to_meV:14.6e}\n")
    print(f"Data saved: {output_file}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
