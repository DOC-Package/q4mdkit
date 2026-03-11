#!/usr/bin/env python3
"""
Plot spectral densities for selected molecules from energy_diff.dat.

Selected molecules for Pentacene: 126, 127, 130, 159, 167, 195, 197, 207, 250
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.integrate import simpson
from scipy.interpolate import interp1d

# =============================================================================
# Physical constants
# =============================================================================
kB_cm = 0.695034800  # cm⁻¹/K (Boltzmann constant)
c_cm_fs = 2.99792458e-5  # cm/fs (speed of light)
two_pi_c = 2 * np.pi * c_cm_fs  # rad·cm/fs

# Unit conversions
eV_to_cm = 8065.54  # cm⁻¹/eV
cm_to_eV = 1.23984e-4  # eV/cm⁻¹
cm_to_meV = 0.123984  # meV/cm⁻¹

# Analysis parameters
T = 300.0  # K
dt = 4.0  # fs
segment_ps = 20.0  # ps
corr_ps = 10.0  # ps
nu_max = 2000  # cm⁻¹

# Selected molecules and their colors (matching VMD ColorIDs)
selected_molecules = [126, 127, 130, 159, 167, 195, 197, 207, 250]
molecule_colors = {
    126: '#FF0000',   # ColorID 1: red
    127: '#FF7F00',   # ColorID 3: orange
    130: '#A020F0',   # ColorID 11: purple
    159: '#00FF00',   # ColorID 7: green
    167: '#CCCC00',   # ColorID 4: yellow (darker)
    195: '#00FFFF',   # ColorID 10: cyan
    197: '#FFC0CB',   # ColorID 9: pink
    207: '#0000FF',   # ColorID 0: blue
    250: '#8B4513',   # ColorID 2: brown
}


# =============================================================================
# File I/O
# =============================================================================
def read_energy_diff(filepath: str) -> dict:
    """Read energy_diff.dat file with multiple molecule columns."""
    frames = []
    times = []
    energies = []
    col_names = []
    mol_ids = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            if not line_stripped:
                continue
            if line_stripped.startswith('#'):
                # Parse column names from header
                if 'dE_mol' in line or 'E_mol' in line:
                    parts = line.replace('#', '').split()
                    for p in parts:
                        if p.startswith('dE_mol') or p.startswith('E_mol'):
                            col_names.append(p)
                            # Extract molecule ID
                            mol_id = int(p.replace('dE_mol', '').replace('E_mol', '').replace('(eV)', ''))
                            mol_ids.append(mol_id)
                continue
            
            parts = line_stripped.split()
            if len(parts) >= 3:
                frames.append(int(parts[0]))
                times.append(float(parts[1]))
                energy_vals = [float(c) for c in parts[2:]]
                energies.append(energy_vals)
    
    return {
        'frame': np.array(frames),
        'time': np.array(times),
        'energies': np.array(energies),  # eV
        'col_names': col_names,
        'mol_ids': mol_ids,
    }


def interpolate_nan(data: np.ndarray) -> np.ndarray:
    """Interpolate NaN values."""
    nan_mask = np.isnan(data)
    if np.sum(nan_mask) == 0:
        return data
    
    valid_mask = ~nan_mask
    if np.sum(valid_mask) < 2:
        return data
    
    indices = np.arange(len(data))
    f_interp = interp1d(indices[valid_mask], data[valid_mask], 
                        kind='linear', fill_value='extrapolate')
    
    data_interp = data.copy()
    data_interp[nan_mask] = f_interp(indices[nan_mask])
    return data_interp


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
# Main
# =============================================================================
print("="*70)
print("Spectral Density Plot - Selected Molecules")
print("="*70)

# Read data file
print("\nReading data files...")
data1 = read_energy_diff('energy_diff.dat')

print(f"  energy_diff.dat: {len(data1['mol_ids'])} molecules: {data1['mol_ids']}")

# Read energy_diff_elstat.dat (single column)
print("  Reading energy_diff_elstat.dat...")
elstat_data = np.loadtxt('energy_diff_elstat.dat', comments='#')
# Column 3 is dE(eV)
elstat_energy = elstat_data[:, 3]  # eV
print(f"    {len(elstat_energy)} frames")

# Setup
segment_length = int(segment_ps * 1000 / dt)
corr_length = int(corr_ps * 1000 / dt)
nu_out = np.linspace(1, nu_max, 500)
kBT = kB_cm * T

print(f"\nAnalysis parameters:")
print(f"  Segment length: {segment_length} frames ({segment_ps} ps)")
print(f"  Correlation length: {corr_length} frames ({corr_ps} ps)")
print(f"  Frequency range: 0-{nu_max} cm⁻¹")

# Compute spectral densities for selected molecules
results = {}

print(f"\nComputing spectral densities for selected molecules...")
for mol_id in selected_molecules:
    print(f"  Processing mol{mol_id}...", end=" ")
    
    # Find which file contains this molecule
    if mol_id in data1['mol_ids']:
        data = data1
        col_idx = data1['mol_ids'].index(mol_id)
    else:
        print(f"NOT FOUND")
        continue
    
    # Get energy data and interpolate NaN
    energy = data['energies'][:, col_idx]
    energy = interpolate_nan(energy)
    
    # Convert to cm⁻¹ and compute fluctuation
    energy_cm = energy * eV_to_cm
    delta = energy_cm - np.mean(energy_cm)
    
    # Compute correlation and spectral density
    t_corr, C, n_seg = compute_correlation_segmented(delta, dt, segment_length, corr_length)
    J = correlation_to_spectral_density(t_corr, C, nu_out, T)
    
    # Reorganization energy
    lambda_cm = C[0] / (2 * kBT)
    lambda_meV = lambda_cm * cm_to_meV
    
    results[mol_id] = {
        'J': J,
        'lambda_cm': lambda_cm,
        'lambda_meV': lambda_meV,
    }
    
    print(f"λ = {lambda_meV:.2f} meV")

# Compute spectral density for elstat data
print("  Processing elstat...", end=" ")
elstat_energy_interp = interpolate_nan(elstat_energy)
elstat_cm = elstat_energy_interp * eV_to_cm
elstat_delta = elstat_cm - np.mean(elstat_cm)
t_corr_elstat, C_elstat, n_seg_elstat = compute_correlation_segmented(elstat_delta, dt, segment_length, corr_length)
J_elstat = correlation_to_spectral_density(t_corr_elstat, C_elstat, nu_out, T)
lambda_elstat_cm = C_elstat[0] / (2 * kBT)
lambda_elstat_meV = lambda_elstat_cm * cm_to_meV
results['elstat'] = {
    'J': J_elstat,
    'lambda_cm': lambda_elstat_cm,
    'lambda_meV': lambda_elstat_meV,
}
print(f"λ = {lambda_elstat_meV:.2f} meV")

# =============================================================================
# Plotting
# =============================================================================
print("\nGenerating plots...")

# Plot 1: All spectral densities together
fig, ax = plt.subplots(figsize=(10, 6))

# Sort molecules by lambda in descending order
sorted_mols = sorted([m for m in selected_molecules if m in results], 
                     key=lambda m: results[m]['lambda_meV'], reverse=True)

for mol_id in sorted_mols:
    r = results[mol_id]
    ax.plot(nu_out, r['J'] * cm_to_meV, linewidth=1.5, color=molecule_colors[mol_id], 
            label=f'mol{mol_id} (λ={r["lambda_meV"]:.2f} meV)')

ax.set_xlabel(r'$\omega$ (cm⁻¹)', fontsize=14)
ax.set_ylabel(r'$J(\omega)$ (meV)', fontsize=14)
ax.set_xlim(0, nu_max)
ax.legend(fontsize=14, loc='upper right', bbox_to_anchor=(1.0, 1.0), framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_title('Spectral Densities - Selected Molecules', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig('spectral_density_selected.png', dpi=300, bbox_inches='tight')
plt.savefig('spectral_density_selected.pdf', bbox_inches='tight')
print("  Saved: spectral_density_selected.png, .pdf")
plt.close()

# Plot 2: Individual subplots (sorted by lambda)
n_mols = len([m for m in selected_molecules if m in results])
n_cols = min(4, n_mols)
n_rows = int(np.ceil(n_mols / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows))
if n_mols == 1:
    axes = np.array([[axes]])
elif n_rows == 1:
    axes = axes.reshape(1, -1)

idx = 0
for mol_id in sorted_mols:
    r = results[mol_id]
    ax = axes[idx // n_cols, idx % n_cols]
    
    ax.plot(nu_out, r['J'] * cm_to_meV, linewidth=1.5, color=molecule_colors[mol_id])
    ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=13)
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=13)
    ax.set_xlim(0, nu_max)
    #ax.set_title(f'mol{mol_id}\nλ = {r["lambda_meV"]:.2f} meV', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    idx += 1

# Hide empty subplots
for j in range(idx, n_rows * n_cols):
    axes[j // n_cols, j % n_cols].set_visible(False)

plt.tight_layout()
plt.savefig('spectral_density_selected_individual.png', dpi=300, bbox_inches='tight')
plt.savefig('spectral_density_selected_individual.pdf', bbox_inches='tight')
print("  Saved: spectral_density_selected_individual.png, .pdf")
plt.close()

# Save data (sorted by lambda)
with open('spectral_density_selected.dat', 'w') as f:
    f.write("# Spectral density analysis - Selected molecules (sorted by λ)\n")
    f.write("# Temperature: 300 K\n")
    f.write("# Reorganization energies:\n")
    for mol_id in sorted_mols:
        r = results[mol_id]
        f.write(f"#   mol{mol_id}: λ = {r['lambda_cm']:.2f} cm⁻¹ = {r['lambda_meV']:.2f} meV\n")
    
    # Header line
    header = "# nu(cm-1)"
    for mol_id in sorted_mols:
        header += f"  J_mol{mol_id}"
    f.write(header + "\n")
    
    for i in range(len(nu_out)):
        line = f"{nu_out[i]:10.2f}"
        for mol_id in sorted_mols:
            line += f" {results[mol_id]['J'][i]:14.6e}"
        f.write(line + "\n")

print("  Saved: spectral_density_selected.dat")

# =============================================================================
# Plot 3: Pair comparisons (with elstat as first panel)
# =============================================================================
pairs = [(128, 260), (166, 232), (169, 228), (200, 204)]

# Sort pairs by the maximum lambda in each pair (descending order)
pairs_sorted = sorted(pairs, key=lambda p: max(results.get(p[0], {'lambda_meV': 0})['lambda_meV'], 
                                                  results.get(p[1], {'lambda_meV': 0})['lambda_meV']), 
                       reverse=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# First panel: elstat (standalone) - top-left (row 0, col 0)
ax = axes[0]
if 'elstat' in results:
    r_elstat = results['elstat']
    J_meV = r_elstat['J'] * cm_to_meV
    ax.plot(nu_out, J_meV, linewidth=2, color='black',
            label=f'Total (λ={r_elstat["lambda_meV"]:.2f} meV)')
    ax.set_ylabel(r'Spectral Density (meV)', fontsize=14)
    ax.set_xlim(0, nu_max)
    ax.set_ylim(0, np.max(J_meV) * 1.25)
    ax.legend(fontsize=14, loc='upper right', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=14)
    #ax.set_title('Electrostatic (Total)', fontsize=12, fontweight='bold')

# Remaining panels: pair comparisons (sorted by lambda)
for idx, (mol_a, mol_b) in enumerate(pairs_sorted):
    ax = axes[idx + 1]  # Start from panel 1 (0 is elstat)
    
    if mol_a in results and mol_b in results:
        r_a = results[mol_a]
        r_b = results[mol_b]
        
        J_a_meV = r_a['J'] * cm_to_meV
        J_b_meV = r_b['J'] * cm_to_meV
        ax.plot(nu_out, J_a_meV, linewidth=1.5, color=molecule_colors[mol_a], linestyle='-',
                label=f'mol{mol_a} (λ={r_a["lambda_meV"]:.2f} meV)')
        ax.plot(nu_out, J_b_meV, linewidth=1.5, color=molecule_colors[mol_b], linestyle='--',
                label=f'mol{mol_b} (λ={r_b["lambda_meV"]:.2f} meV)')
        
        # Axis labels: xlabel only for bottom row (idx >= 2), ylabel only for left column (idx == 2)
        panel_idx = idx + 1  # actual panel index (0 is elstat)
        row = panel_idx // 3
        col = panel_idx % 3
        if row == 1:  # Bottom row
            ax.set_xlabel(r'Wavenumber (cm⁻¹)', fontsize=14)
        if col == 0:  # Left column
            ax.set_ylabel(r'Spectral Density (meV)', fontsize=14)
        ax.set_xlim(0, nu_max)
        ax.set_ylim(0, max(np.max(J_a_meV), np.max(J_b_meV)) * 1.25)
        ax.legend(fontsize=14, loc='upper right', framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='both', labelsize=14)
        #ax.set_title(f'mol{mol_a} vs mol{mol_b}', fontsize=12, fontweight='bold')
    else:
        ax.set_visible(False)

# Hide last unused panel (5 used: 1 elstat + 4 pairs)
axes[5].set_visible(False)

plt.tight_layout()
plt.savefig('spectral_density_pairs.png', dpi=300, bbox_inches='tight')
plt.savefig('spectral_density_pairs.pdf', bbox_inches='tight')
print("  Saved: spectral_density_pairs.png, .pdf")
plt.close()

# Summary (sorted by lambda)
print("\n" + "="*70)
print("Summary (sorted by λ)")
print("="*70)
print(f"{'Molecule':<12} {'λ (cm⁻¹)':>12} {'λ (meV)':>12}")
print("-" * 38)
for mol_id in sorted_mols:
    r = results[mol_id]
    print(f"mol{mol_id:<8} {r['lambda_cm']:>12.2f} {r['lambda_meV']:>12.2f}")

print("\n" + "="*70)
print("Complete!")
print("="*70)
