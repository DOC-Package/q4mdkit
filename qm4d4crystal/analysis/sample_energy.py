"""
Sample qm_energy.dat file to extract energy lines corresponding to trajectory frames.

This module extracts energy lines that match trajectory frame sampling intervals.

Usage:
    python sample_energy.py qm_energy.dat --energy-interval 1 --sample-interval 20
    python sample_energy.py qm_energy.dat --start 100 --n-frames 1000
"""

import argparse
from pathlib import Path
from typing import Optional


def sample_energy_file(
    energy_path: str,
    output_path: Optional[str] = None,
    energy_interval: float = 1.0,
    sample_interval: float = 1.0,
    start_frame: int = 0,
    n_frames: Optional[int] = None,
) -> str:
    """
    Sample qm_energy.dat file, extracting lines that correspond to trajectory frames.
    
    Parameters
    ----------
    energy_path : str
        Path to qm_energy.dat file.
    output_path : str, optional
        Output file path. If None, creates '<input>_sampled.dat'.
    energy_interval : float
        Time interval between energy output steps in qm_energy.dat (in fs).
    sample_interval : float
        Time interval between DCD trajectory frames (in fs).
    start_frame : int
        First trajectory frame to include (0-indexed).
    n_frames : int, optional
        Total number of trajectory frames to process.
        If None, process all frames from start_frame to end.
    
    Returns
    -------
    str
        Path to the output file.
    
    Notes
    -----
    For trajectory frame i, the corresponding energy line index is:
        energy_idx = i * (sample_interval / energy_interval)
    """
    energy_path = Path(energy_path)
    
    if output_path is None:
        output_path = energy_path.parent / f"{energy_path.stem}_sampled.dat"
    else:
        output_path = Path(output_path)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Sampling energy file: {energy_path}")
    print(f"  Energy output interval: {energy_interval} fs")
    print(f"  Sample (DCD) interval: {sample_interval} fs")
    print(f"  Start frame: {start_frame}")
    
    # Read energy file
    with open(energy_path, 'r') as f:
        lines = f.readlines()
    
    # Separate header and data
    header_lines = []
    data_lines = []
    for line in lines:
        if line.startswith('#'):
            header_lines.append(line)
        else:
            data_lines.append(line)
    
    n_energy_steps = len(data_lines)
    print(f"  Total energy steps in file: {n_energy_steps}")
    
    # Calculate step ratio: how many energy lines per trajectory frame
    step_ratio = sample_interval / energy_interval
    if step_ratio != int(step_ratio):
        print(f"  Warning: sample_interval ({sample_interval}) is not a multiple of "
              f"energy_interval ({energy_interval}). Using rounded indices.")
    step_ratio = int(step_ratio)
    
    print(f"  Step ratio: {step_ratio} (energy lines per trajectory frame)")
    
    # Calculate the energy line index for start_frame
    start_energy_idx = start_frame * step_ratio
    if start_energy_idx >= n_energy_steps:
        raise ValueError(f"start_frame {start_frame} exceeds available data "
                        f"(max energy idx: {n_energy_steps - 1})")
    
    # Calculate available trajectory frames from start_frame
    max_available_frames = (n_energy_steps - start_energy_idx - 1) // step_ratio + 1
    
    if n_frames is None:
        n_frames = max_available_frames
    else:
        n_frames = min(n_frames, max_available_frames)
    
    print(f"  Processing {n_frames} frames starting from frame {start_frame}")
    
    # Build list of energy indices corresponding to trajectory frames
    sampled_indices = []
    for i in range(n_frames):
        frame_idx = start_frame + i
        energy_idx = frame_idx * step_ratio
        if energy_idx < n_energy_steps:
            sampled_indices.append(energy_idx)
    
    n_sampled = len(sampled_indices)
    print(f"  Sampled {n_sampled} energy lines")
    
    if n_sampled == 0:
        print("  Error: No frames to process!")
        return str(output_path)
    
    # Write output file
    with open(output_path, 'w') as f:
        # Write header
        for h in header_lines:
            f.write(h)
        # Write sampled data lines
        for idx in sampled_indices:
            f.write(data_lines[idx])
    
    print(f"  -> Saved: {output_path}")
    return str(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Sample qm_energy.dat to match trajectory frames"
    )
    parser.add_argument(
        "energy_file",
        help="Input qm_energy.dat file"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output file path (default: <input>_sampled.dat)"
    )
    parser.add_argument(
        "--energy-interval",
        type=float,
        default=1.0,
        help="Time interval between energy output steps in fs (default: 1)"
    )
    parser.add_argument(
        "--sample-interval",
        type=float,
        default=1.0,
        help="Time interval between DCD trajectory frames in fs (default: 1)"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="First trajectory frame to include (0-indexed, default: 0)"
    )
    parser.add_argument(
        "--n-frames",
        type=int,
        default=None,
        help="Number of frames to process (default: all from start)"
    )
    
    args = parser.parse_args()
    
    sample_energy_file(
        args.energy_file,
        output_path=args.output,
        energy_interval=args.energy_interval,
        sample_interval=args.sample_interval,
        start_frame=args.start,
        n_frames=args.n_frames,
    )


if __name__ == "__main__":
    main()
