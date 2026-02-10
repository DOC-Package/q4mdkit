#!/usr/bin/env python3
"""
Test the effect of windowing on peak sharpness.
Compare Hanning window vs no window.
"""
import numpy as np
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Constants
c_cm_fs = 2.99792458e-5  # cm/fs
dt = 4.0  # fs

print("="*70)
print("Testing Window Function Effect on Peak Sharpness")
print("="*70)

# Load correlation functions
print("\nLoading correlation functions...")
acf_matrix = np.load('acf_matrix.npy')
metadata = np.load('correlation_metadata.npz')

N_modes = metadata['N_modes']
max_lag_full = metadata['max_lag']
mode_freqs = metadata['mode_freqs']

print(f"  N_modes: {N_modes}")
print(f"  max_lag: {max_lag_full}")

# Frequency grid
freq_min, freq_max = 420, 1800
n_freqs = 1381
freq_grid = np.linspace(freq_min, freq_max, n_freqs)

def correlation_to_spectral_density_with_window(corr, dt, freq_grid, max_lag, use_window=True):
    """
    Convert correlation function to spectral density.
    """
    tau = np.arange(max_lag) * dt
    
    if use_window:
        window = np.hanning(max_lag)
        corr_windowed = corr[:max_lag] * window
    else:
        corr_windowed = corr[:max_lag]
    
    J = np.zeros(len(freq_grid))
    for i, nu in enumerate(freq_grid):
        omega = 2 * np.pi * c_cm_fs * nu
        integrand = corr_windowed * np.cos(omega * tau)
        integral = simpson(integrand, x=tau)
        prefactor = 2 * np.pi * c_cm_fs
        J[i] = prefactor * integral
    
    return J

def compute_fwhm(J, freq_grid):
    """Compute FWHM of peak"""
    max_val = np.max(J)
    if max_val <= 0:
        return 0.0, 0.0
    
    half_max = max_val / 2
    peak_idx = np.argmax(J)
    peak_freq = freq_grid[peak_idx]
    
    above_half = J > half_max
    if np.any(above_half):
        indices = np.where(above_half)[0]
        if len(indices) > 1:
            left_idx = indices[0]
            right_idx = indices[-1]
            fwhm = freq_grid[right_idx] - freq_grid[left_idx]
        else:
            fwhm = 0.0
    else:
        fwhm = 0.0
    
    return peak_freq, fwhm

# Test different scenarios
scenarios = [
    {"name": "max_lag=5000 (20ps), with Hanning", "max_lag": 5000, "use_window": True},
    {"name": "max_lag=5000 (20ps), NO window", "max_lag": 5000, "use_window": False},
    {"name": "max_lag=500 (2ps), with Hanning", "max_lag": 500, "use_window": True},
    {"name": "max_lag=500 (2ps), NO window", "max_lag": 500, "use_window": False},
    {"name": "max_lag=100 (0.4ps), with Hanning", "max_lag": 100, "use_window": True},
    {"name": "max_lag=100 (0.4ps), NO window", "max_lag": 100, "use_window": False},
]

# Analyze modes around 520 and 790 cm⁻¹
target_freqs = [520, 790]

for target_freq in target_freqs:
    # Find closest mode
    mode_idx = np.argmin(np.abs(mode_freqs - target_freq))
    actual_freq = mode_freqs[mode_idx]
    
    print(f"\n" + "="*70)
    print(f"Mode {mode_idx} (freq = {actual_freq:.1f} cm⁻¹)")
    print("="*70)
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    axes = axes.flatten()
    
    for i, scenario in enumerate(scenarios):
        max_lag = scenario["max_lag"]
        use_window = scenario["use_window"]
        
        # Compute J with/without window
        J_diag = correlation_to_spectral_density_with_window(
            acf_matrix[mode_idx, :], dt, freq_grid, max_lag, use_window
        )
        
        peak_freq, fwhm = compute_fwhm(J_diag, freq_grid)
        
        print(f"\n{scenario['name']}:")
        print(f"  Peak at: {peak_freq:.1f} cm⁻¹")
        print(f"  Peak height: {np.max(J_diag):.6e}")
        print(f"  FWHM: {fwhm:.2f} cm⁻¹")
        
        # Plot
        ax = axes[i]
        ax.plot(freq_grid, J_diag, 'b-', linewidth=1.5)
        ax.axvline(actual_freq, color='r', linestyle='--', alpha=0.5, 
                   label=f'Mode: {actual_freq:.1f} cm⁻¹')
        
        ax.set_xlabel('Frequency (cm⁻¹)', fontsize=10)
        ax.set_ylabel('J', fontsize=10)
        ax.set_title(f'{scenario["name"]}\nFWHM: {fwhm:.1f} cm⁻¹', 
                     fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # Zoom to peak region
        if peak_freq > 0:
            zoom_width = max(100, fwhm * 3)
            ax.set_xlim(peak_freq - zoom_width, peak_freq + zoom_width)
    
    plt.tight_layout()
    plt.savefig(f'window_effect_mode{mode_idx}_{int(actual_freq)}cm.png', 
                dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: window_effect_mode{mode_idx}_{int(actual_freq)}cm.png")

# Summary comparison
print("\n" + "="*70)
print("Summary: FWHM Comparison")
print("="*70)

summary_data = []
for target_freq in target_freqs:
    mode_idx = np.argmin(np.abs(mode_freqs - target_freq))
    actual_freq = mode_freqs[mode_idx]
    
    print(f"\nMode {mode_idx} ({actual_freq:.1f} cm⁻¹):")
    print(f"{'Scenario':<40} {'FWHM (cm⁻¹)':<15}")
    print("-" * 60)
    
    for scenario in scenarios:
        max_lag = scenario["max_lag"]
        use_window = scenario["use_window"]
        
        J_diag = correlation_to_spectral_density_with_window(
            acf_matrix[mode_idx, :], dt, freq_grid, max_lag, use_window
        )
        
        peak_freq, fwhm = compute_fwhm(J_diag, freq_grid)
        print(f"{scenario['name']:<40} {fwhm:>10.2f}")
        summary_data.append([mode_idx, actual_freq, scenario['name'], fwhm])

print("\n" + "="*70)
print("COMPLETE")
print("="*70)
