#!/usr/bin/env python3
"""
Compute spectral densities from energy_diff.dat and energy_diff_elstat.dat
and calculate their difference and cross-term contribution.

The total spectral density J_tot can be decomposed as:
    J_tot(ω) = J_1(ω) + J_2(ω) + J_cross(ω)
    
where:
    J_1: spectral density from energy_diff (total)
    J_2: spectral density from energy_diff_elstat (electrostatic)
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


def compute_spectral_density_fft(delta: np.ndarray, dt: float, T: float,
                                  segment_length: int, nu_max: float = 2000) -> tuple:
    """
    Compute spectral density using FFT (Welch method).
    This guarantees non-negative spectral density for autocorrelation.
    
    J(ν̃) = (2πc ν̃ / k_B T) * S(ν̃)
    where S(ν̃) is the power spectrum (always non-negative).
    """
    N = len(delta)
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    # Hanning window
    window = np.hanning(segment_length)
    window_norm = np.sum(window**2) / segment_length  # Normalization factor
    
    # Frequency axis (convert to cm^-1)
    freq_fs = fftfreq(segment_length, d=dt)  # 1/fs
    freq_cm = freq_fs / c_cm_fs  # cm^-1
    
    # Only positive frequencies
    pos_mask = freq_cm > 0
    freq_pos = freq_cm[pos_mask]
    
    # Accumulate power spectrum
    power_sum = np.zeros(np.sum(pos_mask))
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        segment = delta[start:end]
        segment = segment - np.mean(segment)  # Remove mean
        
        # Apply window and FFT
        segment_windowed = segment * window
        fft_seg = fft(segment_windowed)
        
        # Power spectrum (one-sided)
        power = np.abs(fft_seg[pos_mask])**2
        power_sum += power
    
    # Average and normalize
    power_avg = power_sum / n_segments
    # Normalization: dt for time integral, 2 for one-sided spectrum, window_norm for window
    power_avg = power_avg * dt * 2 / (segment_length * window_norm)
    
    # Convert to spectral density: J(ν) = 2πc ν / (2 k_B T) * S(ν)
    kBT = kB_cm * T
    J = two_pi_c * freq_pos / (2.0 * kBT) * power_avg
    
    return freq_pos, J, n_segments


def compute_cross_spectral_density_fft(delta1: np.ndarray, delta2: np.ndarray, 
                                        dt: float, T: float,
                                        segment_length: int, nu_max: float = 2000) -> tuple:
    """
    Compute cross-spectral density using FFT.
    Note: Cross-spectrum can have negative real part.
    """
    N = min(len(delta1), len(delta2))
    n_segments = N // segment_length
    
    if n_segments < 1:
        raise ValueError(f"Segment length ({segment_length}) > data length ({N})")
    
    # Hanning window
    window = np.hanning(segment_length)
    window_norm = np.sum(window**2) / segment_length
    
    # Frequency axis
    freq_fs = fftfreq(segment_length, d=dt)
    freq_cm = freq_fs / c_cm_fs
    
    pos_mask = freq_cm > 0
    freq_pos = freq_cm[pos_mask]
    
    # Accumulate cross spectrum (real part)
    cross_sum = np.zeros(np.sum(pos_mask))
    
    for i in range(n_segments):
        start = i * segment_length
        end = start + segment_length
        seg1 = delta1[start:end] - np.mean(delta1[start:end])
        seg2 = delta2[start:end] - np.mean(delta2[start:end])
        
        seg1_windowed = seg1 * window
        seg2_windowed = seg2 * window
        
        fft_seg1 = fft(seg1_windowed)
        fft_seg2 = fft(seg2_windowed)
        
        # Cross spectrum: Re[conj(FFT1) * FFT2]
        cross = np.real(np.conj(fft_seg1[pos_mask]) * fft_seg2[pos_mask])
        cross_sum += cross
    
    cross_avg = cross_sum / n_segments
    cross_avg = cross_avg * dt * 2 / (segment_length * window_norm)
    
    kBT = kB_cm * T
    J = two_pi_c * freq_pos / (2.0 * kBT) * cross_avg
    
    return freq_pos, J, n_segments


# =============================================================================
# Main analysis
# =============================================================================
def main():
    """Main analysis routine."""
    # Input parameters
    input_file1 = 'energy_diff.dat'  # Total energy difference
    input_file2 = 'energy_diff_elstat.dat'  # Electrostatic contribution
    dt = 4.0  # fs
    T = 300.0  # K
    
    # Analysis parameters
    segment_ps = 100.0
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
    
    # Compute spectral densities using FFT (guarantees non-negative for auto-spectra)
    nu_max = 2000
    
    print(f"\nComputing spectral densities (FFT method)...")
    nu_out, J1, _ = compute_spectral_density_fft(delta1, dt, T, segment_length)  # Total
    _, J2, _ = compute_spectral_density_fft(delta2, dt, T, segment_length)  # Elstat
    _, J_diff, _ = compute_spectral_density_fft(delta_diff, dt, T, segment_length)  # Non-elstat
    _, J_cross1, _ = compute_cross_spectral_density_fft(delta_diff, delta2, dt, T, segment_length)  # Cross 1
    _, J_cross2, _ = compute_cross_spectral_density_fft(delta2, delta_diff, dt, T, segment_length)  # Cross 2
    
    # Limit to nu_max
    mask = nu_out <= nu_max
    nu_out = nu_out[mask]
    J1 = J1[mask]
    J2 = J2[mask]
    J_diff = J_diff[mask]
    J_cross1 = J_cross1[mask]
    J_cross2 = J_cross2[mask]
    
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
    
    # Plot correlation functions
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: Correlation functions
    ax = axes[0]
    t_corr_ps = t_corr / 1000  # fs -> ps
    ax.plot(t_corr_ps, C1, 'k-', linewidth=2, label='Total')
    ax.plot(t_corr_ps, C2, 'b-', linewidth=1.5, label='Electrostatic')
    ax.plot(t_corr_ps, C_diff, 'r-', linewidth=1.5, label='Non-electrostatic')
    ax.plot(t_corr_ps, C_cross1, 'g--', linewidth=1.5, label='Cross1 (non,elstat)')
    ax.plot(t_corr_ps, C_cross2, 'm:', linewidth=1.5, label='Cross2 (elstat,non)')
    ax.set_xlabel('Time (ps)', fontsize=14)
    ax.set_ylabel(r'$C(\tau)$ (cm$^{-2}$)', fontsize=14)
    ax.set_xlim(0, corr_ps)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Correlation Functions', fontsize=13, fontweight='bold')
    
    # Right: Zoomed correlation functions (first 5 ps)
    ax = axes[1]
    zoom_ps = 5.0
    zoom_idx = int(zoom_ps * 1000 / dt)
    ax.plot(t_corr_ps[:zoom_idx], C1[:zoom_idx], 'k-', linewidth=2, label='Total')
    ax.plot(t_corr_ps[:zoom_idx], C2[:zoom_idx], 'b-', linewidth=1.5, label='Electrostatic')
    ax.plot(t_corr_ps[:zoom_idx], C_diff[:zoom_idx], 'r-', linewidth=1.5, label='Non-electrostatic')
    ax.plot(t_corr_ps[:zoom_idx], C_cross1[:zoom_idx], 'g--', linewidth=1.5, label='Cross1')
    ax.plot(t_corr_ps[:zoom_idx], C_cross2[:zoom_idx], 'm:', linewidth=1.5, label='Cross2')
    ax.set_xlabel('Time (ps)', fontsize=14)
    ax.set_ylabel(r'$C(\tau)$ (cm$^{-2}$)', fontsize=14)
    ax.set_xlim(0, zoom_ps)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Correlation Functions (zoomed)', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('correlation_functions.png', dpi=300, bbox_inches='tight')
    plt.savefig('correlation_functions.pdf', bbox_inches='tight')
    print(f"\nFigures saved: correlation_functions.png, .pdf")
    plt.close()
    
    # Plot spectral densities
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: All spectral densities
    ax = axes[0]
    ax.plot(nu_out, J1, 'k-', linewidth=2, label='Total')
    ax.plot(nu_out, J2, 'b-', linewidth=1.5, label='Electrostatic')
    ax.plot(nu_out, J_diff, 'r-', linewidth=1.5, label='Non-electrostatic')
    ax.plot(nu_out, J_cross1, 'g--', linewidth=1.5, label='Cross1 (non,elstat)')
    ax.plot(nu_out, J_cross2, 'm:', linewidth=1.5, label='Cross2 (elstat,non)')
    ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=14)
    ax.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=14)
    ax.set_xlim(0, nu_max)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Spectral Densities', fontsize=13, fontweight='bold')
    
    # Right: Difference (J1 - J2 = J_diff + J_cross1 + J_cross2)
    ax = axes[1]
    J_diff_direct = J1 - J2
    J_reconstructed = J_diff + J_cross1 + J_cross2
    ax.plot(nu_out, J_diff_direct, 'k-', linewidth=2, label=r'$J_{total} - J_{elstat}$')
    ax.plot(nu_out, J_reconstructed, 'r--', linewidth=1.5, label=r'$J_{non} + J_{cross1} + J_{cross2}$')
    ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=14)
    ax.set_ylabel(r'$J(\omega)$ (cm⁻¹)', fontsize=14)
    ax.set_xlim(0, nu_max)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_title('Spectral Density Difference', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('spectral_density_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('spectral_density_comparison.pdf', bbox_inches='tight')
    print(f"\nFigures saved: spectral_density_comparison.png, .pdf")
    plt.close()
    
    # Save data
    with open('spectral_density_comparison.dat', 'w') as f:
        f.write("# Spectral density comparison\n")
        f.write(f"# Temperature: {T} K\n")
        f.write(f"# λ_total = {lambda1*cm_to_eV:.4f} eV\n")
        f.write(f"# λ_elstat = {lambda2*cm_to_eV:.4f} eV\n")
        f.write(f"# λ_non-elstat = {lambda_diff*cm_to_eV:.4f} eV\n")
        f.write(f"# λ_cross1 (non,elstat) = {lambda_cross1*cm_to_eV:.4f} eV\n")
        f.write(f"# λ_cross2 (elstat,non) = {lambda_cross2*cm_to_eV:.4f} eV\n")
        f.write("# nu(cm-1)  J_total  J_elstat  J_non-elstat  J_cross1  J_cross2  J_diff(direct)\n")
        for i in range(len(nu_out)):
            f.write(f"{nu_out[i]:10.2f} {J1[i]:14.6e} {J2[i]:14.6e} {J_diff[i]:14.6e} {J_cross1[i]:14.6e} {J_cross2[i]:14.6e} {J_diff_direct[i]:14.6e}\n")
    print(f"Data saved: spectral_density_comparison.dat")
    
    print("\n=== Analysis complete ===")


if __name__ == '__main__':
    main()
