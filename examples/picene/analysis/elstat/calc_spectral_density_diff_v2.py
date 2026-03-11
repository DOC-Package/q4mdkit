#!/usr/bin/env python3
"""
Compute spectral densities from energy_diff.dat and energy_diff_elstat.dat
and calculate their difference and cross-term contribution.

The total spectral density J_tot can be decomposed as:
    J_tot(ω) = J_1(ω) + J_2(ω) + J_cross(ω)
    
where:
    J_1: spectral density from energy_diff (total)
    J_2: spectral density from energy_diff_elstat (external)
    J_cross: cross-term contribution
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
cm_to_meV = 0.123984  # meV/cm⁻¹


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
                delta_epsilon.append(float(parts[2]))  # dE(a.u.) column
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'delta_epsilon': np.array(delta_epsilon),  # Ha
    }


def interpolate_nan(data: dict) -> dict:
    """Interpolate NaN values in the delta_epsilon array."""
    delta_epsilon = data['delta_epsilon']
    time = data['time']
    
    nan_mask = np.isnan(delta_epsilon)
    n_nan = np.sum(nan_mask)
    
    if n_nan == 0:
        return data
    
    print(f"  Found {n_nan} NaN values, interpolating...")
    
    valid_mask = ~nan_mask
    if np.sum(valid_mask) < 2:
        raise ValueError("Not enough valid data points for interpolation")
    
    f_interp = interp1d(time[valid_mask], delta_epsilon[valid_mask], 
                        kind='linear', fill_value='extrapolate')
    
    delta_epsilon_interp = delta_epsilon.copy()
    delta_epsilon_interp[nan_mask] = f_interp(time[nan_mask])
    
    data['delta_epsilon'] = delta_epsilon_interp
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
    # Input parameters
    input_file1 = 'energy_diff.dat'  # Total energy difference
    input_file2 = 'energy_diff_elstat.dat'  # External contribution
    dt = 4.0  # fs
    T = 300.0  # K
    
    # Analysis parameters
    segment_ps = 20.0
    corr_ps = 10.0
    
    print(f"=== Spectral Density Difference Analysis ===")
    print(f"File 1 (total): {input_file1}")
    print(f"File 2 (elstat): {input_file2}")
    print(f"Time step: {dt} fs, Temperature: {T} K")
    
    # Read data
    data1 = read_energy_diff(input_file1)
    data2 = read_energy_diff(input_file2)
    
    # Interpolate NaN values
    data1 = interpolate_nan(data1)
    data2 = interpolate_nan(data2)
    
    # Match data lengths
    n_points = min(len(data1['delta_epsilon']), len(data2['delta_epsilon']))
    print(f"\nUsing {n_points} frames")
    
    # Convert to cm⁻¹ and compute fluctuations
    eps1_cm = data1['delta_epsilon'][:n_points] * Ha_to_cm
    eps2_cm = data2['delta_epsilon'][:n_points] * Ha_to_cm
    
    delta1 = eps1_cm - np.mean(eps1_cm)  # Total fluctuation
    delta2 = eps2_cm - np.mean(eps2_cm)  # Elstat fluctuation
    delta_diff = delta1 - delta2  # Non-elstat fluctuation
    
    print(f"\nStatistics (eV):")
    print(f"  Total:     mean = {np.mean(eps1_cm)*cm_to_eV:.4f}, std = {np.std(delta1)*cm_to_eV:.4f}")
    print(f"  Elstat:    mean = {np.mean(eps2_cm)*cm_to_eV:.4f}, std = {np.std(delta2)*cm_to_eV:.4f}")
    print(f"  Non-elstat: std = {np.std(delta_diff)*cm_to_eV:.4f}")
    
    # Compute correlation functions
    segment_length = int(segment_ps * 1000 / dt)
    corr_length = int(corr_ps * 1000 / dt)
    
    print(f"\nComputing correlation functions...")
    t_corr, C1, n_seg = compute_correlation_segmented(delta1, dt, segment_length, corr_length)
    _, C2, _ = compute_correlation_segmented(delta2, dt, segment_length, corr_length)
    _, C_diff, _ = compute_correlation_segmented(delta_diff, dt, segment_length, corr_length)
    # Cross-correlation 1: ⟨δ_non-elstat(t) × δ_elstat(t+τ)⟩
    _, C_cross1, _ = compute_cross_correlation_segmented(delta_diff, delta2, dt, segment_length, corr_length)
    # Cross-correlation 2: ⟨δ_elstat(t) × δ_non-elstat(t+τ)⟩
    _, C_cross2, _ = compute_cross_correlation_segmented(delta2, delta_diff, dt, segment_length, corr_length)
    
    print(f"  {n_seg} segments used")
    
    # Compute spectral densities
    nu_max = 2000
    nu_out = np.linspace(1, nu_max, 500)
    
    print(f"\nComputing spectral densities...")
    J1 = correlation_to_spectral_density(t_corr, C1, nu_out, T)  # Total
    J2 = correlation_to_spectral_density(t_corr, C2, nu_out, T)  # Elstat
    J_diff = correlation_to_spectral_density(t_corr, C_diff, nu_out, T)  # Non-elstat
    J_cross1 = correlation_to_spectral_density(t_corr, C_cross1, nu_out, T)  # Cross-term 1
    J_cross2 = correlation_to_spectral_density(t_corr, C_cross2, nu_out, T)  # Cross-term 2
    
    # Reorganization energies
    kBT = kB_cm * T
    lambda1 = C1[0] / (2 * kBT)
    lambda2 = C2[0] / (2 * kBT)
    lambda_diff = C_diff[0] / (2 * kBT)
    lambda_cross1 = C_cross1[0] / (2 * kBT)
    lambda_cross2 = C_cross2[0] / (2 * kBT)
    
    print(f"\nReorganization energies:")
    print(f"  λ_total = {lambda1:.2f} cm⁻¹ = {lambda1*cm_to_eV:.4f} eV")
    print(f"  λ_elstat = {lambda2:.2f} cm⁻¹ = {lambda2*cm_to_eV:.4f} eV")
    print(f"  λ_non-elstat = {lambda_diff:.2f} cm⁻¹ = {lambda_diff*cm_to_eV:.4f} eV")
    print(f"  λ_cross1 (non-elstat, elstat) = {lambda_cross1:.2f} cm⁻¹ = {lambda_cross1*cm_to_eV:.4f} eV")
    print(f"  λ_cross2 (elstat, non-elstat) = {lambda_cross2:.2f} cm⁻¹ = {lambda_cross2*cm_to_eV:.4f} eV")
    print(f"  λ_elstat + λ_non-elstat + λ_cross1 + λ_cross2 = {(lambda2 + lambda_diff + lambda_cross1 + lambda_cross2)*cm_to_eV:.4f} eV")
    
    # Plot spectral densities
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Spectral densities (convert to meV)
    ax.plot(nu_out, J1 * cm_to_meV, 'k-', linewidth=2, label='Total')
    ax.plot(nu_out, J2 * cm_to_meV, 'b-', linewidth=1.5, label='External')
    ax.plot(nu_out, J_diff * cm_to_meV, 'r-', linewidth=1.5, label='Internal')
    ax.plot(nu_out, J_cross1 * cm_to_meV, 'g--', linewidth=1.5, label='Cross1 (int,ext)')
    ax.plot(nu_out, J_cross2 * cm_to_meV, 'm:', linewidth=1.5, label='Cross2 (ext,int)')
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=16)
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=16)
    ax.set_xlim(0, nu_max)
    ax.tick_params(axis='both', labelsize=16)
    ax.legend(fontsize=16)
    ax.grid(True, alpha=0.3)
    #ax.set_title('Spectral Densities', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('spectral_density_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('spectral_density_comparison.pdf', bbox_inches='tight')
    print(f"\nFigures saved: spectral_density_comparison.png, .pdf")
    plt.close()
    
    # Save data
    J_diff_direct = J1 - J2
    with open('spectral_density_comparison.dat', 'w') as f:
        f.write("# Spectral density comparison\n")
        f.write(f"# Temperature: {T} K\n")
        f.write(f"# λ_total = {lambda1*cm_to_meV:.2f} meV\n")
        f.write(f"# λ_external = {lambda2*cm_to_meV:.2f} meV\n")
        f.write(f"# λ_internal = {lambda_diff*cm_to_meV:.2f} meV\n")
        f.write(f"# λ_cross1 (int,ext) = {lambda_cross1*cm_to_meV:.2f} meV\n")
        f.write(f"# λ_cross2 (ext,int) = {lambda_cross2*cm_to_meV:.2f} meV\n")
        f.write("# nu(cm-1)  J_total(meV)  J_external(meV)  J_internal(meV)  J_cross1(meV)  J_cross2(meV)  J_diff(direct)(meV)\n")
        for i in range(len(nu_out)):
            f.write(f"{nu_out[i]:10.2f} {J1[i]*cm_to_meV:14.6e} {J2[i]*cm_to_meV:14.6e} {J_diff[i]*cm_to_meV:14.6e} {J_cross1[i]*cm_to_meV:14.6e} {J_cross2[i]*cm_to_meV:14.6e} {J_diff_direct[i]*cm_to_meV:14.6e}\n")
    print(f"Data saved: spectral_density_comparison.dat")
    
    print("\n=== Analysis complete ===")


if __name__ == '__main__':
    main()
