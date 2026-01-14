"""
Sample qm_energy.dat file to extract energy lines corresponding to trajectory frames.

This module extracts energy lines that match trajectory frame sampling intervals.

Usage:
    python sample_energy.py qm_energy.dat --dt-fs 1 --dt-sample-fs 20
    python sample_energy.py qm_energy.dat --start 100 --n-frames 1000
"""

import argparse
from pathlib import Path
from typing import Optional


def sample_energy_file(
    energy_path: str,
    output_path: Optional[str] = None,
    dt_fs: float = 1.0,
    dt_sample_fs: float = 1.0,
    start_frame: int = 0,
    n_frames: Optional[int] = None,
    t0_fs: float = 0.0,
) -> str:
    """
    Sample qm_energy.dat file, extracting lines that correspond to trajectory frames.
    
    Parameters
    ----------
    energy_path : str
        Path to qm_energy.dat file.
    output_path : str, optional
        Output file path. If None, creates '<input>_sampled.dat'.
    dt_fs : float
        Time interval between energy output steps in qm_energy.dat (in fs).
    dt_sample_fs : float
        Time interval between DCD trajectory frames (in fs).
    start_frame : int
        First trajectory frame to include (0-indexed).
    n_frames : int, optional
        Total number of trajectory frames to process.
        If None, process all frames from start_frame to end.
    t0_fs : float
        Initial time in fs for first output frame (default: 0.0).
    
    Returns
    -------
    str
        Path to the output file.
    
    Notes
    -----
    For trajectory frame i, the corresponding energy line index is:
        energy_idx = i * (dt_sample_fs / dt_fs)
    """
    energy_path = Path(energy_path)
    
    if output_path is None:
        output_path = energy_path.parent / f"{energy_path.stem}_sampled.dat"
    else:
        output_path = Path(output_path)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Sampling energy file: {energy_path}")
    print(f"  Energy output interval: {dt_fs} fs")
    print(f"  Sample (DCD) interval: {dt_sample_fs} fs")
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
    step_ratio = dt_sample_fs / dt_fs
    if step_ratio != int(step_ratio):
        raise ValueError(f"dt_sample_fs ({dt_sample_fs}) must be a multiple of "
                        f"dt_fs ({dt_fs}). Got ratio: {step_ratio}")
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
        f.write("# call_index  time[fs]  energy[a.u.]\n")
        # Write sampled data lines with call_index, time and total energy
        for i, idx in enumerate(sampled_indices):
            time_fs = t0_fs + i * dt_sample_fs
            # Extract call_index and total energy (columns: index, electronic, repulsive, total)
            parts = data_lines[idx].strip().split()
            call_index = parts[0]
            total_energy = parts[3]  # Total energy is the 4th column
            f.write(f"{call_index:>8}  {time_fs:12.3f}  {total_energy}\n")
    
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
        "--dt-fs",
        type=float,
        default=1.0,
        help="Time interval between energy output steps in fs (default: 1)"
    )
    parser.add_argument(
        "--dt-sample-fs",
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
    parser.add_argument(
        "--t0",
        type=float,
        default=0.0,
        help="Initial time in fs for first output frame (default: 0.0)"
    )
    
    args = parser.parse_args()
    
    sample_energy_file(
        args.energy_file,
        output_path=args.output,
        dt_fs=args.dt_fs,
        dt_sample_fs=args.dt_sample_fs,
        start_frame=args.start,
        n_frames=args.n_frames,
        t0_fs=args.t0,
    )


if __name__ == "__main__":
    main()
