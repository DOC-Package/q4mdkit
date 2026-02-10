from q4mdkit.analysis import convert_dcd_to_xyz

# Input files
trajectory = "mol1.dcd"
topology = "mol1.pdb"
    
# Output XYZ file
output_xyz = "mol1.xyz"
    
# Convert DCD to XYZ (extract first frame by default)
convert_dcd_to_xyz(
    trajectory_file=trajectory,
    topology_file=topology,
    output_file=output_xyz,
    frame=0,  # First frame (can use -1 for last frame, etc.)
    verbose=True
)
