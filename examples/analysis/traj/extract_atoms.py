from q4mdkit.analysis.extract_atoms import extract_atoms

# Input files
trajectory = "nve.dcd"
topology = "pentacene.pdb"
    
# Extract molecule 1 (qmatoms1)
extract_atoms(
    trajectory_path=trajectory,
    topology_path=topology,
    atom_indices="qmatoms1",
    output_path="mol1.dcd",
    output_topology="mol1.pdb"
)
