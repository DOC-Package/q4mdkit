"""
Trajectory splitting module.

Splits DCD trajectory files into multiple parts for parallel processing.
Supports skipping initial equilibration frames and limiting total frames.
"""

import os
import argparse
import mdtraj as md
from typing import Optional


def get_n_frames(trajectory_path: str) -> int:
    """Get the number of frames in a trajectory file.
    
    Parameters
    ----------
    trajectory_path : str
        Path to the trajectory file (DCD format).
    
    Returns
    -------
    int
        Number of frames in the trajectory.
    """
    traj = md.load(trajectory_path, top=None)
    return traj.n_frames


def split_trajectory(
    trajectory_path: str,
    topology_path: str,
    n_splits: int,
    output_dir: str,
    start_frame: int = 0,
    n_frames: Optional[int] = None
) -> list[str]:
    """Split a trajectory file into multiple parts.
    
    Parameters
    ----------
    trajectory_path : str
        Path to the trajectory file (DCD format).
    topology_path : str
        Path to the topology file (PDB or GRO format).
    n_splits : int
        Number of parts to split the trajectory into.
    output_dir : str
        Directory to save the split trajectory files.
    start_frame : int, optional
        Starting frame index (0-indexed). Default is 0.
        Use this to skip initial equilibration frames.
    n_frames : int, optional
        Number of frames to process from start_frame.
        If None, process all remaining frames.
    
    Returns
    -------
    list[str]
        List of paths to the split trajectory files.
    
    Example
    -------
    # Skip first 100 frames and process next 1000 frames, split into 10 parts
    >>> split_trajectory("prod.dcd", "prod.pdb", 10, "split", start_frame=100, n_frames=1000)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Load full trajectory
    traj = md.load(trajectory_path, top=topology_path)
    total_frames = traj.n_frames
    
    # Determine frame range
    if start_frame >= total_frames:
        raise ValueError(f"start_frame ({start_frame}) >= total frames ({total_frames})")
    
    end_frame = total_frames
    if n_frames is not None:
        end_frame = min(start_frame + n_frames, total_frames)
    
    actual_frames = end_frame - start_frame
    print(f"Trajectory: {total_frames} total frames")
    print(f"Processing: frames {start_frame} to {end_frame-1} ({actual_frames} frames)")
    
    # Calculate frames per split
    frames_per_split = actual_frames // n_splits
    remainder = actual_frames % n_splits
    
    output_paths = []
    current_frame = start_frame
    
    for i in range(n_splits):
        # Distribute remainder frames across first splits
        split_frames = frames_per_split + (1 if i < remainder else 0)
        
        if split_frames == 0:
            continue
        
        # Extract frames for this split
        split_traj = traj[current_frame:current_frame + split_frames]
        
        # Save split trajectory
        output_path = os.path.join(output_dir, f"traj_part{i+1:02d}.dcd")
        split_traj.save_dcd(output_path)
        output_paths.append(output_path)
        
        print(f"Part {i+1}: frames {current_frame}-{current_frame + split_frames - 1} "
              f"({split_frames} frames) -> {output_path}")
        
        current_frame += split_frames
    
    print(f"\nSplit into {len(output_paths)} files")
    return output_paths


def main():
    """CLI for trajectory splitting."""
    parser = argparse.ArgumentParser(
        description="Split trajectory files into multiple parts"
    )
    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Path to trajectory file (DCD)"
    )
    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Path to topology file (PDB/GRO)"
    )
    parser.add_argument(
        "-n", "--n-splits",
        type=int,
        required=True,
        help="Number of splits"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="split",
        help="Output directory (default: split)"
    )
    parser.add_argument(
        "--start-frame",
        type=int,
        default=0,
        help="Starting frame index (0-indexed, default: 0)"
    )
    parser.add_argument(
        "--n-frames",
        type=int,
        default=None,
        help="Number of frames to process (default: all)"
    )
    
    args = parser.parse_args()
    
    split_trajectory(
        trajectory_path=args.trajectory,
        topology_path=args.topology,
        n_splits=args.n_splits,
        output_dir=args.output_dir,
        start_frame=args.start_frame,
        n_frames=args.n_frames
    )


if __name__ == "__main__":
    main()
