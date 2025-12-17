"""
Split a DCD trajectory file into multiple parts.

Usage:
    python split_dcd.py input.dcd topology.gro -n 5
    python split_dcd.py input.dcd topology.gro -n 5 -o output_dir
"""

import argparse
import mdtraj as md
from pathlib import Path
import numpy as np


def get_n_frames(traj_path: str, top_path: str) -> int:
    """Get total number of frames without loading entire trajectory."""
    # Load just first frame to get structure info
    traj = md.load_frame(str(traj_path), 0, top=str(top_path))
    
    # Count frames by iterating (memory efficient)
    n_frames = 0
    for chunk in md.iterload(str(traj_path), top=str(top_path), chunk=100):
        n_frames += chunk.n_frames
    return n_frames


def split_trajectory(traj_path: str, top_path: str, n_splits: int, output_dir: str = None):
    """
    Split a DCD trajectory into n_splits equal parts.
    Memory-efficient: processes trajectory in chunks.
    
    Parameters:
    -----------
    traj_path : str
        Path to input DCD trajectory file
    top_path : str
        Path to topology file (gro, pdb, etc.)
    n_splits : int
        Number of parts to split into
    output_dir : str, optional
        Output directory (default: same as input)
    """
    traj_path = Path(traj_path)
    top_path = Path(top_path)
    
    if output_dir is None:
        output_dir = traj_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading trajectory: {traj_path}")
    print(f"Topology: {top_path}")
    
    # Get total frames (memory efficient)
    print("Counting frames...")
    n_frames = get_n_frames(str(traj_path), str(top_path))
    print(f"Total frames: {n_frames}")
    print(f"Splitting into {n_splits} parts")
    
    if n_splits > n_frames:
        raise ValueError(f"Cannot split {n_frames} frames into {n_splits} parts")
    
    # Calculate frame ranges for each split
    frames_per_split = n_frames // n_splits
    remainder = n_frames % n_splits
    
    split_ranges = []
    start = 0
    for i in range(n_splits):
        n_this = frames_per_split + (1 if i < remainder else 0)
        split_ranges.append((start, start + n_this))
        start += n_this
    
    print(f"Frames per split: ~{frames_per_split}")
    print("-" * 50)
    
    # Base name for output files
    base_name = traj_path.stem
    
    # Process each split
    for i, (start_frame, end_frame) in enumerate(split_ranges):
        output_file = output_dir / f"{base_name}_part{i+1:02d}.dcd"
        n_frames_this = end_frame - start_frame
        
        print(f"Part {i+1}/{n_splits}: frames {start_frame}-{end_frame-1} ({n_frames_this} frames)")
        
        # Collect frames for this split
        frames_collected = []
        current_frame = 0
        
        for chunk in md.iterload(str(traj_path), top=str(top_path), chunk=100):
            chunk_start = current_frame
            chunk_end = current_frame + chunk.n_frames
            
            # Check if this chunk overlaps with our target range
            if chunk_end > start_frame and chunk_start < end_frame:
                # Calculate local indices within this chunk
                local_start = max(0, start_frame - chunk_start)
                local_end = min(chunk.n_frames, end_frame - chunk_start)
                
                frames_collected.append(chunk[local_start:local_end])
            
            current_frame = chunk_end
            
            # Stop if we've passed our target range
            if current_frame >= end_frame:
                break
        
        # Join and save
        if frames_collected:
            traj_split = md.join(frames_collected)
            traj_split.save_dcd(str(output_file))
            print(f"  -> Saved: {output_file.name}")
        
    print("-" * 50)
    print(f"Done! {n_splits} DCD files saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Split a DCD trajectory file into multiple parts"
    )
    parser.add_argument(
        "trajectory",
        help="Input DCD trajectory file"
    )
    parser.add_argument(
        "topology",
        help="Topology file (gro, pdb, etc.)"
    )
    parser.add_argument(
        "-n", "--n_splits",
        type=int,
        required=True,
        help="Number of parts to split into"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output directory (default: same as input)"
    )
    
    args = parser.parse_args()
    
    split_trajectory(
        args.trajectory,
        args.topology,
        args.n_splits,
        args.output
    )


if __name__ == "__main__":
    main()
