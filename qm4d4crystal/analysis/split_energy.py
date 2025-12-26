"""
Split qm_energy.dat file into multiple parts.

This module extracts energy lines that correspond to trajectory frames,
with support for specifying start position.

Usage:
    python split_energy.py qm_energy.dat -n 10 --energy-interval 1 --sample-interval 20
    python split_energy.py qm_energy.dat -n 10 --start 100
"""

import argparse
from pathlib import Path
from typing import Optional, List


def split_energy_file(
    energy_path: str,
    n_splits: int,
    output_dir: str,
    energy_interval: float = 1.0,
    sample_interval: float = 1.0,
    start_frame: int = 0,
    n_frames: Optional[int] = None,
) -> None:
    """
    Split qm_energy.dat file into multiple parts, extracting lines that
    correspond to trajectory frames.
    
    Parameters
    ----------
    energy_path : str
        Path to qm_energy.dat file.
    n_splits : int
        Number of parts to split into.
    output_dir : str
        Output directory.
    energy_interval : float
        Time interval between energy output steps in qm_energy.dat (in fs).
        This is the MD timestep.
    sample_interval : float
        Time interval between DCD trajectory frames (in fs).
        e.g., if trajectory is saved every 20 fs, set to 20.
    start_frame : int
        First trajectory frame to include (0-indexed).
        Used to skip initial equilibration.
    n_frames : int, optional
        Total number of trajectory frames to process.
        If None, process all frames from start_frame to end.
    
    Notes
    -----
    For trajectory frame i, the corresponding energy line index is:
        energy_idx = i * (sample_interval / energy_interval)
    
    Example: if sample_interval=20 and energy_interval=1,
        frame 0 -> energy line 0
        frame 1 -> energy line 20
        frame 2 -> energy line 40
    """
    energy_path = Path(energy_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Splitting energy file: {energy_path}")
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
    # Frame (start_frame + i) -> energy line (start_frame + i) * step_ratio
    traj_energy_indices = []
    for i in range(n_frames):
        frame_idx = start_frame + i
        energy_idx = frame_idx * step_ratio
        if energy_idx < n_energy_steps:
            traj_energy_indices.append(energy_idx)
    
    n_traj_frames = len(traj_energy_indices)
    print(f"  Trajectory frames with energy data: {n_traj_frames}")
    
    if n_traj_frames == 0:
        print("  Error: No frames to process!")
        return
    
    # Calculate split ranges (in trajectory frame units, relative to processed data)
    frames_per_split = n_traj_frames // n_splits
    remainder = n_traj_frames % n_splits
    
    split_ranges = []
    start = 0
    for i in range(n_splits):
        n_this = frames_per_split + (1 if i < remainder else 0)
        split_ranges.append((start, start + n_this))
        start += n_this
    
    # Base name for output files
    base_name = energy_path.stem
    
    # Process each split
    for i, (split_start, split_end) in enumerate(split_ranges):
        output_file = output_dir / f"{base_name}_part{i+1:02d}.dat"
        
        # Get energy indices for this split's trajectory frames
        split_energy_indices = traj_energy_indices[split_start:split_end]
        n_steps_this = len(split_energy_indices)
        
        # Calculate actual frame numbers
        actual_start_frame = start_frame + split_start
        actual_end_frame = start_frame + split_end - 1
        
        print(f"  Part {i+1}/{n_splits}: frames {actual_start_frame}-{actual_end_frame}, "
              f"energy lines {split_energy_indices[0]}-{split_energy_indices[-1]} "
              f"({n_steps_this} lines)")
        
        # Write output file
        with open(output_file, 'w') as f:
            # Write header
            for h in header_lines:
                f.write(h)
            # Write data lines for this split (only lines corresponding to traj frames)
            for idx in split_energy_indices:
                f.write(data_lines[idx])
        
        print(f"    -> Saved: {output_file.name}")
    
    print(f"  Done! {n_splits} energy files saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Split qm_energy.dat file into multiple parts"
    )
    parser.add_argument(
        "energy_file",
        help="Input qm_energy.dat file"
    )
    parser.add_argument(
        "-n", "--n-splits",
        type=int,
        required=True,
        help="Number of parts to split into"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output directory (default: same as input)"
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
    
    output_dir = args.output
    if output_dir is None:
        output_dir = str(Path(args.energy_file).parent)
    
    split_energy_file(
        args.energy_file,
        args.n_splits,
        output_dir,
        energy_interval=args.energy_interval,
        sample_interval=args.sample_interval,
        start_frame=args.start,
        n_frames=args.n_frames,
    )


if __name__ == "__main__":
    main()
