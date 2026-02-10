#!/usr/bin/env python3
"""
Visualization for quadratic form spectral fitting results.
Loads results from spectral_fit_cross.py and generates plots.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Unit conversion: cm⁻¹ to meV
CM_TO_MEV = 0.12398  # 1 cm⁻¹ = 0.12398 meV


def load_results():
    """Load spectral fitting results from files."""
    # Load spectral densities
    data = np.loadtxt('spectral_density_quadratic_form.dat')
    freq_grid = data[:, 0]
    J_total = data[:, 1]
    J_model = data[:, 2]
    J_diag = data[:, 3]
    J_offdiag = data[:, 4]
    residual = data[:, 5]
    
    # Load coefficients and extract metadata
    mode_indices = []
    mode_freqs = []
    g_opt = []
    R2 = None
    capture = None
    
    with open('g_k_quadratic_form.dat') as f:
        for line in f:
            if line.startswith('#'):
                if 'R²' in line or 'R^2' in line:
                    try:
                        R2 = float(line.split('=')[1].strip())
                    except:
                        pass
                elif 'Capture' in line:
                    try:
                        capture = float(line.split('=')[1].replace('%', '').strip())
                    except:
                        pass
            else:
                parts = line.split()
                if len(parts) >= 3:
                    mode_indices.append(int(parts[0]))
                    mode_freqs.append(float(parts[1]))
                    g_opt.append(float(parts[2]))
    
    mode_indices = np.array(mode_indices)
    mode_freqs = np.array(mode_freqs)
    g_opt = np.array(g_opt)
    
    # Calculate statistics if not found in header
    J_total_integral = np.trapezoid(J_total, freq_grid)
    J_model_integral = np.trapezoid(J_model, freq_grid)
    J_diag_integral = np.trapezoid(J_diag, freq_grid)
    J_offdiag_integral = np.trapezoid(J_offdiag, freq_grid)
    
    if R2 is None:
        ss_res = np.sum(residual**2)
        ss_tot = np.sum((J_total - np.mean(J_total))**2)
        R2 = 1 - ss_res / ss_tot
    
    if capture is None:
        capture = 100 * J_model_integral / J_total_integral
    
    return {
        'freq_grid': freq_grid,
        'J_total': J_total,
        'J_model': J_model,
        'J_diag': J_diag,
        'J_offdiag': J_offdiag,
        'residual': residual,
        'mode_indices': mode_indices,
        'mode_freqs': mode_freqs,
        'g_opt': g_opt,
        'R2': R2,
        'capture': capture,
        'J_total_integral': J_total_integral,
        'J_model_integral': J_model_integral,
        'J_diag_integral': J_diag_integral,
        'J_offdiag_integral': J_offdiag_integral,
    }


def plot_spectral_comparison(results=None, output='spectral_comparison.png', 
                             unit='meV', figsize=(10, 6)):
    """
    Plot spectral density comparison (Panel 1 only).
    
    Parameters
    ----------
    results : dict, optional
        Results from load_results(). If None, loads from files.
    output : str
        Output filename.
    unit : str
        Unit for y-axis: 'meV' or 'cm-1'.
    figsize : tuple
        Figure size.
    """
    if results is None:
        results = load_results()
    
    freq_grid = results['freq_grid']
    J_total = results['J_total']
    J_model = results['J_model']
    R2 = results['R2']
    capture = results['capture']
    
    # Unit conversion
    if unit == 'meV':
        # J has units of cm⁻¹ (energy), convert to meV
        J_total_plot = J_total * CM_TO_MEV
        J_model_plot = J_model * CM_TO_MEV
        ylabel = 'J(ω) (meV)'
    else:
        J_total_plot = J_total
        J_model_plot = J_model
        ylabel = 'J(ω) (cm⁻¹)'
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(freq_grid, J_total_plot, 'k-', lw=2, label='Total', alpha=0.8)
    ax.plot(freq_grid, J_model_plot, 'r-', lw=1.5, label=f'Model (Linear)', alpha=0.8)
    ax.fill_between(freq_grid, 0, J_model_plot, alpha=0.2, color='red')
    
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    #ax.set_title(f'Quadratic Form Fitting: J(ω) = g$^T$ J$_{{qq}}$(ω) g (Capture={capture:.1f}%)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=13, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(freq_grid[0], freq_grid[-1])
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


# ============================================================================
# Main script
# ============================================================================
if __name__ == '__main__':
    print("="*70)
    print("Visualization: Quadratic Form Spectral Fitting")
    print("="*70)

    # ============================================================================
    # Load results
    # ============================================================================
    print("\nLoading results...")
    results = load_results()
    
    freq_grid = results['freq_grid']
    J_total = results['J_total']
    J_model = results['J_model']
    J_diag = results['J_diag']
    J_offdiag = results['J_offdiag']
    residual = results['residual']
    mode_indices = results['mode_indices']
    mode_freqs = results['mode_freqs']
    g_opt = results['g_opt']
    R2 = results['R2']
    capture = results['capture']
    J_total_integral = results['J_total_integral']
    J_diag_integral = results['J_diag_integral']
    J_offdiag_integral = results['J_offdiag_integral']
    
    print(f"  Frequency grid: {len(freq_grid)} points ({freq_grid[0]:.0f}-{freq_grid[-1]:.0f} cm⁻¹)")
    print(f"  Modes: {len(g_opt)}")
    print(f"  R² = {R2:.4f}")
    print(f"  Capture = {capture:.1f}%")

    # ============================================================================
    # Panel 1 only: Spectral comparison (meV)
    # ============================================================================
    print("\nGenerating spectral comparison figure (meV)...")
    plot_spectral_comparison(results, output='spectral_comparison_meV.png', unit='meV')

    # ============================================================================
    # Main Figure: 4-panel overview
    # ============================================================================
    print("Generating main figure...")

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # Panel 1: Spectral density comparison
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8)
    ax1.plot(freq_grid, J_model, 'r-', lw=1.5, label=f'J_model (R²={R2:.4f})', alpha=0.8)
    ax1.fill_between(freq_grid, 0, J_model, alpha=0.2, color='red')
    ax1.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax1.set_title(f'Quadratic Form Fitting: J(ω) = g^T J_qq(ω) g (Capture={capture:.1f}%)', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper right')
    ax1.grid(True, alpha=0.3)

    # Panel 2: Residual
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(freq_grid, residual, 'b-', lw=1.5, label='Residual', alpha=0.7)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax2.set_ylabel('Residual (cm⁻¹²)', fontsize=11)
    ax2.set_title('Fitting Residual', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Diagonal vs Off-diagonal contributions
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(freq_grid, J_diag, 'g-', lw=1.5, 
             label=f'Diagonal ({100*J_diag_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax3.plot(freq_grid, J_offdiag, 'm-', lw=1.5, 
             label=f'Off-diagonal ({100*J_offdiag_integral/J_total_integral:.1f}%)', alpha=0.7)
    ax3.axhline(0, color='k', ls='--', lw=0.5)
    ax3.set_xlabel('Frequency (cm⁻¹)', fontsize=11)
    ax3.set_ylabel('J(ν) (cm⁻¹²)', fontsize=11)
    ax3.set_title('Contribution Decomposition', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)

    # Panel 4: Stacked contributions
    ax4 = fig.add_subplot(gs[2, :])
    ax4.plot(freq_grid, J_total, 'k-', lw=2, label='J_total (data)', alpha=0.8, zorder=3)
    ax4.fill_between(freq_grid, 0, J_diag, alpha=0.4, color='green', label='Diagonal', zorder=1)
    ax4.fill_between(freq_grid, J_diag, J_diag + J_offdiag, alpha=0.4, color='magenta', 
                     label='Off-diagonal', zorder=2)
    ax4.plot(freq_grid, J_model, 'r--', lw=1.5, label='Total model', alpha=0.8, zorder=3)
    ax4.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax4.set_ylabel('J(ν) (cm⁻¹²)', fontsize=12)
    ax4.set_title('Stacked Contributions', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=11)
    ax4.grid(True, alpha=0.3)

    plt.savefig('spectral_fit_quadratic_form.png', dpi=150, bbox_inches='tight')
    print(f"Saved: spectral_fit_quadratic_form.png")
    plt.close()

    # ============================================================================
    # Coefficients Figure
    # ============================================================================
    print("Generating coefficients figure...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: All coefficients
    ax = axes[0, 0]
    colors = ['green' if g > 0 else 'red' for g in g_opt]
    ax.bar(range(len(g_opt)), g_opt, color=colors, alpha=0.7)
    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.set_xlabel('Mode index (sorted)', fontsize=11)
    ax.set_ylabel('g_k', fontsize=11)
    ax.set_title('All Coupling Coefficients', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    # Panel 2: Coefficient vs frequency
    ax = axes[0, 1]
    ax.scatter(mode_freqs, g_opt, c=colors, s=50, alpha=0.7)
    ax.axhline(0, color='k', ls='--', lw=0.5)
    ax.set_xlabel('Mode Frequency (cm⁻¹)', fontsize=11)
    ax.set_ylabel('g_k', fontsize=11)
    ax.set_title('Coefficients vs Mode Frequency', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Panel 3: |g_k| distribution
    ax = axes[1, 0]
    ax.hist(np.abs(g_opt), bins=30, color='blue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('|g_k|', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Distribution of |g_k|', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Panel 4: Top modes
    ax = axes[1, 1]
    top_idx = np.argsort(np.abs(g_opt))[-15:][::-1]
    top_freqs = mode_freqs[top_idx]
    top_g = g_opt[top_idx]
    top_colors = ['green' if g > 0 else 'red' for g in top_g]
    bars = ax.barh(range(len(top_idx)), np.abs(top_g), color=top_colors, alpha=0.7)
    ax.set_yticks(range(len(top_idx)))
    ax.set_yticklabels([f'Mode {mode_indices[i]} ({mode_freqs[i]:.0f} cm⁻¹)' for i in top_idx])
    ax.set_xlabel('|g_k|', fontsize=11)
    ax.set_title('Top 15 Modes by |g_k|', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig('spectral_fit_quadratic_form_coefficients.png', dpi=150, bbox_inches='tight')
    print(f"Saved: spectral_fit_quadratic_form_coefficients.png")
    plt.close()

    # ============================================================================
    # Summary
    # ============================================================================
    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70)
    print(f"\nGenerated figures:")
    print(f"  - spectral_comparison_meV.png (spectral comparison in meV)")
    print(f"  - spectral_fit_quadratic_form.png (main overview)")
    print(f"  - spectral_fit_quadratic_form_coefficients.png (coefficient analysis)")
    print(f"\nFit statistics:")
    print(f"  R² = {R2:.4f}")
    print(f"  Capture = {capture:.1f}%")
    print(f"  Diagonal contribution: {100*J_diag_integral/J_total_integral:.1f}%")
    print(f"  Off-diagonal contribution: {100*J_offdiag_integral/J_total_integral:.1f}%")
