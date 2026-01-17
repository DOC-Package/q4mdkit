#!/usr/bin/env python3
"""
Plot time evolution of selected normal mode coordinates.

Usage:
    python plot_mode_time_evolution.py [--modes 7 8 9] [--dt 0.5] [--combined]
"""

import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from q4mdkit.analysis.normal_mode_analysis import (
    plot_mode_time_evolution,
    plot_mode_time_evolution_combined,
)


def main():
    parser = argparse.ArgumentParser(description='Plot mode time evolution')
    parser.add_argument('--input', '-i', default='mode_coords.npz',
                       help='Input mode coordinates file (default: mode_coords.npz)')
    parser.add_argument('--modes', '-m', type=int, nargs='+', default=None,
                       help='Mode indices to plot (1-indexed). Default: first 6 modes')
    parser.add_argument('--dt', type=float, default=1.0,
                       help='Time step between frames (default: 1.0)')
    parser.add_argument('--time-unit', '-u', default='frame',
                       help='Time unit label (default: frame)')
    parser.add_argument('--output', '-o', default='mode_time_evolution',
                       help='Output file prefix (default: mode_time_evolution)')
    parser.add_argument('--combined', '-c', action='store_true',
                       help='Also create combined plot with all modes')
    parser.add_argument('--normalize', '-n', action='store_true',
                       help='Normalize modes in combined plot')
    parser.add_argument('--no-mean', action='store_true',
                       help='Do not show mean line in subplots')
    parser.add_argument('--share-y', action='store_true',
                       help='Share y-axis across subplots')
    
    args = parser.parse_args()
    
    # Load data
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    print(f"Loading data from {args.input}")
    data = np.load(args.input)
    mode_coords = data['mode_coords']
    mode_indices = data['mode_indices'] + 1  # Convert to 1-indexed
    
    n_frames, n_modes = mode_coords.shape
    print(f"  Frames: {n_frames}")
    print(f"  Modes: {n_modes} (indices {mode_indices[0]} to {mode_indices[-1]})")
    
    # Select modes
    if args.modes is None:
        # Default: first 6 vibrational modes
        selected_modes = list(mode_indices[:6])
    else:
        selected_modes = args.modes
    
    print(f"\nSelected modes: {selected_modes}")
    
    # Plot individual subplots
    output_file = f"{args.output}.png"
    print(f"\nCreating time evolution plot...")
    
    fig = plot_mode_time_evolution(
        mode_coords,
        selected_modes,
        mode_indices=mode_indices,
        dt=args.dt,
        time_unit=args.time_unit,
        output_file=output_file,
        show_mean=not args.no_mean,
        share_y=args.share_y
    )
    
    # Combined plot
    if args.combined:
        combined_file = f"{args.output}_combined.png"
        print(f"\nCreating combined plot...")
        
        fig2 = plot_mode_time_evolution_combined(
            mode_coords,
            selected_modes,
            mode_indices=mode_indices,
            dt=args.dt,
            time_unit=args.time_unit,
            output_file=combined_file,
            normalize=args.normalize
        )
    
    print("\nDone!")


if __name__ == "__main__":
    main()
