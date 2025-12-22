# Visualize one layer of pentacene crystal
# Load structure
mol new pentacene.pdb

# Check z-coordinate range
set all [atomselect top "all"]
set minmax [measure minmax $all]
puts "Full structure z-range: [lindex $minmax 0 2] to [lindex $minmax 1 2] Angstrom"

# Extract one layer (bottom layer: z < 15)
# Adjust the z-range as needed
set layer [atomselect top "z > -3 and z < 15"]
puts "Number of atoms in layer: [$layer num]"

# Delete default representation
mol delrep 0 top

# Add representation for the layer
mol representation CPK 1.0 0.3 12 12
mol color Name
mol selection "z > -3 and z < 12"
mol material Opaque
mol addrep top

# Optional: Add bonds
mol representation DynamicBonds 1.6 0.1 12
mol color Name
mol selection "z > -3 and z < 12"
mol material Opaque
mol addrep top

# Center view on the layer
display resetview

puts "Layer visualization complete"
puts "To extract different layers, modify the z-range in the selection"
puts "Middle layer example: z > 12 and z < 28"
puts "Top layer example: z > 25 and z < 42"
