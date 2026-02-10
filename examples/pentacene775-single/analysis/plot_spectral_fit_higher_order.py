#!/usr/bin/env python3
"""
Visualization for higher-order spectral fitting results.
Loads results from spectral_fit_higher_order_selected.py and generates plots.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Unit conversion: cm⁻¹ to meV
CM_TO_MEV = 0.12398  # 1 cm⁻¹ = 0.12398 meV


def load_results():
    """Load higher-order spectral fitting results from files."""
    # Load spectral densities
    data = np.loadtxt('spectral_density_higher_order_selected.dat')
    freq_grid = data[:, 0]
    J_total = data[:, 1]
    J_model_2nd = data[:, 2]
    J_3rd = data[:, 3]
    J_4th = data[:, 4]
    J_model_total = data[:, 5]
    residual_2nd = data[:, 6]
    residual_final = data[:, 7]
    
    # Load d_kl coefficients
    mode_pairs = []
    pair_info = {}
    d_kl = []
    R2_residual = None
    R2_total = None
    selected_modes = []
    
    with open('d_kl_higher_order_selected.dat') as f:
        for line in f:
            if line.startswith('#'):
                if 'R² (residual)' in line or 'R^2 (residual)' in line:
                    try:
                        R2_residual = float(line.split('=')[1].strip())
                    except:
                        pass
                elif 'R² (total)' in line or 'R^2 (total)' in line:
                    try:
                        R2_total = float(line.split('=')[1].strip())
                    except:
                        pass
                elif 'Selected modes:' in line:
                    try:
                        modes_str = line.split(':')[1].strip()
                        modes_str = modes_str.replace('[', '').replace(']', '')
                        selected_modes = [int(m.strip()) for m in modes_str.split(',')]
                    except:
                        pass
            else:
                parts = line.split()
                if len(parts) >= 5:
                    mode_k = int(parts[0])
                    mode_l = int(parts[1])
                    freq_k = float(parts[2])
                    freq_l = float(parts[3])
                    d_val = float(parts[4])
                    target = int(parts[5]) if len(parts) > 5 else 0
                    match_type = parts[6] if len(parts) > 6 else 'none'
                    
                    idx = len(mode_pairs)
                    mode_pairs.append((mode_k, mode_l))
                    pair_info[idx] = {
                        'mode_i': mode_k, 'mode_j': mode_l,
                        'freq_i': freq_k, 'freq_j': freq_l,
                        'target': target, 'type': match_type
                    }
                    d_kl.append(d_val)
    
    d_kl = np.array(d_kl)
    
    # Calculate statistics if not found in header
    J_total_integral = np.trapezoid(J_total, freq_grid)
    J_model_total_integral = np.trapezoid(J_model_total, freq_grid)
    J_2nd_integral = np.trapezoid(J_model_2nd, freq_grid)
    J_3rd_integral = np.trapezoid(J_3rd, freq_grid)
    J_4th_integral = np.trapezoid(J_4th, freq_grid)
    
    if R2_total is None:
        ss_res = np.sum((J_total - J_model_total)**2)
        ss_tot = np.sum((J_total - np.mean(J_total))**2)
        R2_total = 1 - ss_res / ss_tot
    
    if R2_residual is None:
        ss_res = np.sum(residual_final**2)
        ss_tot = np.sum(residual_2nd**2)
        R2_residual = 1 - ss_res / ss_tot
    
    return {
        'freq_grid': freq_grid,
        'J_total': J_total,
        'J_model_2nd': J_model_2nd,
        'J_3rd': J_3rd,
        'J_4th': J_4th,
        'J_model_total': J_model_total,
        'residual_2nd': residual_2nd,
        'residual_final': residual_final,
        'mode_pairs': mode_pairs,
        'pair_info': pair_info,
        'd_kl': d_kl,
        'R2_residual': R2_residual,
        'R2_total': R2_total,
        'selected_modes': selected_modes,
        'J_total_integral': J_total_integral,
        'J_model_total_integral': J_model_total_integral,
        'J_2nd_integral': J_2nd_integral,
        'J_3rd_integral': J_3rd_integral,
        'J_4th_integral': J_4th_integral,
    }


def plot_total_comparison(results=None, output='spectral_higher_total.png',
                          unit='meV', figsize=(10, 6)):
    """
    Plot total spectral density comparison.
    
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
    J_model_total = results['J_model_total']
    R2_total = results['R2_total']
    selected_modes = results['selected_modes']
    
    # Unit conversion
    if unit == 'meV':
        J_total_plot = J_total * CM_TO_MEV
        J_model_plot = J_model_total * CM_TO_MEV
        ylabel = r'$J(\omega)$ (meV)'
    else:
        J_total_plot = J_total
        J_model_plot = J_model_total
        ylabel = r'$J(\omega)$ (cm$^{-1}$)'
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(freq_grid, J_total_plot, 'k-', lw=2, label=r'$J_{\mathrm{total}}$ (data)', alpha=0.8)
    ax.plot(freq_grid, J_model_plot, 'r-', lw=1.5, label=rf'$J_{{\mathrm{{model}}}}$ ($R^2$={R2_total:.4f})', alpha=0.8)
    ax.fill_between(freq_grid, 0, J_model_plot, alpha=0.2, color='red')
    
    ax.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(r'Total Fit: 2nd + 3rd + 4th Order', fontsize=13, fontweight='bold')
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(freq_grid[0], freq_grid[-1])
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


def plot_residual_comparison(results=None, output='spectral_higher_residual.png',
                             unit='meV', figsize=(10, 6), target_freqs=None):
    """
    Plot residual comparison (before and after higher-order fit).
    
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
    target_freqs : list, optional
        Target frequencies to mark with vertical lines.
    """
    if results is None:
        results = load_results()
    
    if target_freqs is None:
        target_freqs = [524, 800, 1092, 1231]
    
    freq_grid = results['freq_grid']
    residual_2nd = results['residual_2nd']
    residual_final = results['residual_final']
    R2_residual = results['R2_residual']
    
    # Unit conversion
    if unit == 'meV':
        residual_2nd_plot = residual_2nd * CM_TO_MEV
        residual_final_plot = residual_final * CM_TO_MEV
        ylabel = r'Residual (meV)'
    else:
        residual_2nd_plot = residual_2nd
        residual_final_plot = residual_final
        ylabel = r'Residual (cm$^{-1}$)'
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(freq_grid, residual_2nd_plot, 'b-', lw=1.5, label='2nd order residual', alpha=0.7)
    ax.plot(freq_grid, residual_final_plot, 'g-', lw=1.5, label='After 3rd+4th fit', alpha=0.7)
    ax.axhline(0, color='k', ls='--', lw=0.5)
    
    for target in target_freqs:
        ax.axvline(target, color='r', ls=':', alpha=0.5)
    
    ax.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(rf'Residual Comparison ($R^2$={R2_residual:.4f})', fontsize=13, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(freq_grid[0], freq_grid[-1])
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


def plot_contributions(results=None, output='spectral_higher_contributions.png',
                       unit='meV', figsize=(12, 8), target_freqs=None):
    """
    Plot 3rd and 4th order contributions separately.
    
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
    target_freqs : list, optional
        Target frequencies to mark with vertical lines.
    """
    if results is None:
        results = load_results()
    
    if target_freqs is None:
        target_freqs = [524, 800, 1092, 1231]
    
    freq_grid = results['freq_grid']
    J_3rd = results['J_3rd']
    J_4th = results['J_4th']
    J_total_integral = results['J_total_integral']
    J_3rd_integral = results['J_3rd_integral']
    J_4th_integral = results['J_4th_integral']
    
    # Unit conversion
    if unit == 'meV':
        J_3rd_plot = J_3rd * CM_TO_MEV
        J_4th_plot = J_4th * CM_TO_MEV
        ylabel = r'$J(\omega)$ (meV)'
    else:
        J_3rd_plot = J_3rd
        J_4th_plot = J_4th
        ylabel = r'$J(\omega)$ (cm$^{-1}$)'
    
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    
    # 3rd order
    ax1 = axes[0]
    ax1.plot(freq_grid, J_3rd_plot, 'm-', lw=1.5, alpha=0.8)
    ax1.axhline(0, color='k', ls='--', lw=0.5)
    for target in target_freqs:
        ax1.axvline(target, color='r', ls=':', alpha=0.5)
    ax1.set_ylabel(ylabel, fontsize=13)
    pct_3rd = 100 * J_3rd_integral / J_total_integral
    ax1.set_title(rf'3rd Order Contribution ({pct_3rd:.1f}\%)', fontsize=13)
    ax1.grid(True, alpha=0.3)
    
    # 4th order
    ax2 = axes[1]
    ax2.plot(freq_grid, J_4th_plot, 'c-', lw=1.5, alpha=0.8)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    for target in target_freqs:
        ax2.axvline(target, color='r', ls=':', alpha=0.5)
    ax2.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=13)
    ax2.set_ylabel(ylabel, fontsize=13)
    pct_4th = 100 * J_4th_integral / J_total_integral
    ax2.set_title(rf'4th Order Contribution ({pct_4th:.1f}\%)', fontsize=13)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


def plot_stacked(results=None, output='spectral_higher_stacked.png',
                 unit='meV', figsize=(10, 6)):
    """
    Plot stacked contributions (2nd, 3rd, 4th order).
    
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
    J_model_2nd = results['J_model_2nd']
    J_3rd = results['J_3rd']
    J_model_total = results['J_model_total']
    
    # Unit conversion
    if unit == 'meV':
        J_total_plot = J_total * CM_TO_MEV
        J_2nd_plot = J_model_2nd * CM_TO_MEV
        J_3rd_plot = J_3rd * CM_TO_MEV
        J_total_model_plot = J_model_total * CM_TO_MEV
        ylabel = r'$J(\omega)$ (meV)'
    else:
        J_total_plot = J_total
        J_2nd_plot = J_model_2nd
        J_3rd_plot = J_3rd
        J_total_model_plot = J_model_total
        ylabel = r'$J(\omega)$ (cm$^{-1}$)'
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(freq_grid, J_total_plot, 'k-', lw=2, label=r'$J_{\mathrm{total}}$ (data)', alpha=0.8)
    ax.fill_between(freq_grid, 0, J_2nd_plot, alpha=0.4, color='blue', label='2nd order')
    ax.fill_between(freq_grid, J_2nd_plot, J_2nd_plot + J_3rd_plot,
                    alpha=0.4, color='magenta', label='3rd order')
    ax.fill_between(freq_grid, J_2nd_plot + J_3rd_plot, J_total_model_plot,
                    alpha=0.4, color='cyan', label='4th order')
    
    ax.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title('Stacked Contributions', fontsize=13, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(freq_grid[0], freq_grid[-1])
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


def plot_d_kl_coefficients(results=None, output='spectral_higher_d_kl.png',
                           figsize=(14, 6)):
    """
    Plot d_kl coefficient bar chart.
    
    Parameters
    ----------
    results : dict, optional
        Results from load_results(). If None, loads from files.
    output : str
        Output filename.
    figsize : tuple
        Figure size.
    """
    if results is None:
        results = load_results()
    
    mode_pairs = results['mode_pairs']
    pair_info = results['pair_info']
    d_kl = results['d_kl']
    N_pairs = len(mode_pairs)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    x_labels = [rf"{pair_info[i]['mode_i']}$\times${pair_info[i]['mode_j']}" for i in range(N_pairs)]
    colors = ['green' if d > 0 else 'red' for d in d_kl]
    ax.bar(range(N_pairs), d_kl, color=colors, alpha=0.7)
    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.set_xticks(range(N_pairs))
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_xlabel('Mode pair', fontsize=13)
    ax.set_ylabel(r'$d_{kl}$', fontsize=13)
    ax.set_title(rf'$d_{{kl}}$ Coefficients ({N_pairs} pairs)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


def plot_all_panels(results=None, output='spectral_fit_higher_order_selected.png',
                    unit='meV', figsize=(14, 12), target_freqs=None):
    """
    Generate the full 6-panel figure (original layout).
    
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
    target_freqs : list, optional
        Target frequencies to mark with vertical lines.
    """
    if results is None:
        results = load_results()
    
    if target_freqs is None:
        target_freqs = [524, 800, 1092, 1231]
    
    freq_grid = results['freq_grid']
    J_total = results['J_total']
    J_model_2nd = results['J_model_2nd']
    J_3rd = results['J_3rd']
    J_4th = results['J_4th']
    J_model_total = results['J_model_total']
    residual_2nd = results['residual_2nd']
    residual_final = results['residual_final']
    mode_pairs = results['mode_pairs']
    pair_info = results['pair_info']
    d_kl = results['d_kl']
    R2_residual = results['R2_residual']
    R2_total = results['R2_total']
    selected_modes = results['selected_modes']
    J_total_integral = results['J_total_integral']
    J_3rd_integral = results['J_3rd_integral']
    J_4th_integral = results['J_4th_integral']
    N_pairs = len(mode_pairs)
    
    # Unit conversion
    if unit == 'meV':
        scale = CM_TO_MEV
        ylabel = r'$J(\omega)$ (meV)'
        ylabel_res = 'Residual (meV)'
    else:
        scale = 1.0
        ylabel = r'$J(\omega)$ (cm$^{-1}$)'
        ylabel_res = r'Residual (cm$^{-1}$)'
    
    fig, axes = plt.subplots(3, 2, figsize=figsize)
    
    # Panel 1: Original vs Total model
    ax1 = axes[0, 0]
    ax1.plot(freq_grid, J_total * scale, 'k-', lw=2, label=r'$J_{\mathrm{total}}$ (data)', alpha=0.8)
    ax1.plot(freq_grid, J_model_total * scale, 'r-', lw=1.5, 
             label=rf'$J_{{\mathrm{{model}}}}$ ($R^2$={R2_total:.4f})', alpha=0.8)
    ax1.set_ylabel(ylabel)
    ax1.set_title(rf'Total Fit: 2nd + 3rd + 4th Order')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Residual comparison
    ax2 = axes[0, 1]
    ax2.plot(freq_grid, residual_2nd * scale, 'b-', lw=1.5, label='2nd order residual', alpha=0.7)
    ax2.plot(freq_grid, residual_final * scale, 'g-', lw=1.5, label='After 3rd+4th fit', alpha=0.7)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    for target in target_freqs:
        ax2.axvline(target, color='r', ls=':', alpha=0.5)
    ax2.set_ylabel(ylabel_res)
    ax2.set_title(rf'Residual Comparison ($R^2$={R2_residual:.4f})')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: 3rd order contribution
    ax3 = axes[1, 0]
    ax3.plot(freq_grid, J_3rd * scale, 'm-', lw=1.5, alpha=0.8)
    ax3.axhline(0, color='k', ls='--', lw=0.5)
    for target in target_freqs:
        ax3.axvline(target, color='r', ls=':', alpha=0.5)
    ax3.set_ylabel(ylabel)
    pct_3rd = 100 * J_3rd_integral / J_total_integral
    ax3.set_title(rf'3rd Order Contribution ({pct_3rd:.1f}\%)')
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: 4th order contribution
    ax4 = axes[1, 1]
    ax4.plot(freq_grid, J_4th * scale, 'c-', lw=1.5, alpha=0.8)
    ax4.axhline(0, color='k', ls='--', lw=0.5)
    for target in target_freqs:
        ax4.axvline(target, color='r', ls=':', alpha=0.5)
    ax4.set_ylabel(ylabel)
    pct_4th = 100 * J_4th_integral / J_total_integral
    ax4.set_title(rf'4th Order Contribution ({pct_4th:.1f}\%)')
    ax4.grid(True, alpha=0.3)
    
    # Panel 5: Stacked contributions
    ax5 = axes[2, 0]
    ax5.plot(freq_grid, J_total * scale, 'k-', lw=2, label=r'$J_{\mathrm{total}}$ (data)', alpha=0.8)
    ax5.fill_between(freq_grid, 0, J_model_2nd * scale, alpha=0.4, color='blue', label='2nd order')
    ax5.fill_between(freq_grid, J_model_2nd * scale, (J_model_2nd + J_3rd) * scale,
                     alpha=0.4, color='magenta', label='3rd order')
    ax5.fill_between(freq_grid, (J_model_2nd + J_3rd) * scale, J_model_total * scale,
                     alpha=0.4, color='cyan', label='4th order')
    ax5.set_xlabel(r'Frequency (cm$^{-1}$)')
    ax5.set_ylabel(ylabel)
    ax5.set_title('Stacked Contributions')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Panel 6: d_kl bar chart
    ax6 = axes[2, 1]
    x_labels = [rf"{pair_info[i]['mode_i']}$\times${pair_info[i]['mode_j']}" for i in range(N_pairs)]
    colors = ['green' if d > 0 else 'red' for d in d_kl]
    ax6.bar(range(N_pairs), d_kl, color=colors, alpha=0.7)
    ax6.axhline(0, color='k', ls='-', lw=0.5)
    ax6.set_xticks(range(N_pairs))
    ax6.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax6.set_xlabel('Mode pair')
    ax6.set_ylabel(r'$d_{kl}$')
    ax6.set_title(rf'$d_{{kl}}$ Coefficients ({N_pairs} pairs)')
    ax6.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved: {output}")
    plt.close()
    
    return fig


# ============================================================================
# Main script
# ============================================================================
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Plot higher-order spectral fitting results')
    parser.add_argument('--unit', choices=['meV', 'cm-1'], default='meV',
                        help='Unit for y-axis (default: meV)')
    parser.add_argument('--all', action='store_true',
                        help='Generate all individual plots')
    parser.add_argument('--output', type=str, default=None,
                        help='Output filename for combined plot')
    args = parser.parse_args()
    
    print("Loading results...")
    results = load_results()
    
    print(f"\nResults summary:")
    print(f"  R² (total): {results['R2_total']:.4f}")
    print(f"  R² (residual): {results['R2_residual']:.4f}")
    print(f"  Mode pairs: {len(results['mode_pairs'])}")
    print(f"  Selected modes: {results['selected_modes']}")
    
    # Generate combined plot
    output = args.output or 'spectral_fit_higher_order_selected.png'
    plot_all_panels(results, output=output, unit=args.unit)
    
    if args.all:
        print("\nGenerating individual plots...")
        plot_total_comparison(results, unit=args.unit)
        plot_residual_comparison(results, unit=args.unit)
        plot_contributions(results, unit=args.unit)
        plot_stacked(results, unit=args.unit)
        plot_d_kl_coefficients(results)
    
    print("\nDone.")
