#!/usr/bin/env python3
"""Plot TD-DFTB absorption spectrum for Quinacridone."""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('absorption_spectrum.dat', comments='#')
energy = data[:, 0]
wavelength = data[:, 1]
intensity = data[:, 2]
std = data[:, 3]

# Find range with significant intensity
max_int = intensity.max()
mask = intensity > 0.001 * max_int
if mask.any():
    e_min = energy[mask].min()
    e_max = energy[mask].max()
    print(f"Energy range with signal: {e_min:.3f} - {e_max:.3f} eV")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Energy plot (0-6 eV range)
plot_mask = (energy >= 0.0) & (energy <= 6.0)
ax1.plot(energy[plot_mask], intensity[plot_mask], 'b-', lw=1.5, label='Average')
ax1.fill_between(energy[plot_mask], 
                 intensity[plot_mask] - std[plot_mask], 
                 intensity[plot_mask] + std[plot_mask],
                 alpha=0.3, color='blue', label='± Std. dev.')

ax1.set_xlabel('Energy (eV)', fontsize=12)
ax1.set_ylabel('Intensity (arb. units)', fontsize=12)
ax1.set_title('TD-DFTB Absorption Spectrum (Quinacridone)', fontsize=14)
ax1.legend()
ax1.set_xlim(0.0, 6.0)
ax1.grid(True, alpha=0.3)

# Wavelength plot (200-800 nm range)
plot_mask2 = (wavelength >= 200) & (wavelength <= 800) & np.isfinite(wavelength)
ax2.plot(wavelength[plot_mask2], intensity[plot_mask2], 'r-', lw=1.5, label='Average')
ax2.fill_between(wavelength[plot_mask2], 
                 intensity[plot_mask2] - std[plot_mask2], 
                 intensity[plot_mask2] + std[plot_mask2],
                 alpha=0.3, color='red', label='± Std. dev.')

ax2.set_xlabel('Wavelength (nm)', fontsize=12)
ax2.set_ylabel('Intensity (arb. units)', fontsize=12)
ax2.set_title('TD-DFTB Absorption Spectrum', fontsize=14)
ax2.legend()
ax2.set_xlim(200, 800)
ax2.invert_xaxis()  # Higher wavelength on left
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('absorption_spectrum.png', dpi=150)
print("Saved: absorption_spectrum.png")

# Print peak info
peak_idx = np.argmax(intensity)
print(f"\nPeak absorption:")
print(f"  Energy: {energy[peak_idx]:.3f} eV")
print(f"  Wavelength: {wavelength[peak_idx]:.1f} nm")
print(f"  Intensity: {intensity[peak_idx]:.6f}")

# Also check raw excitation energies
exc_data = np.loadtxt('excitations.dat', comments='#')
print(f"\nExcitation statistics:")
print(f"  Energy range: {exc_data[:,3].min():.4f} - {exc_data[:,3].max():.4f} eV")
# S1 (state 1) statistics
s1 = exc_data[exc_data[:,2] == 1]
print(f"  S1 energy: {s1[:,3].mean():.3f} ± {s1[:,3].std():.3f} eV")
