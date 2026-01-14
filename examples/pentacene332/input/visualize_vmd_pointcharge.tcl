# VMD visualization script - QM + MM point charges
# QM region: CPK representation
# MM region: Small spheres (point charge representation)
# Usage: vmd -e /home/takahashi/python/qm4d4crystal/examples/pentacene775/input/visualize_vmd_pointcharge.tcl

# Load structure
mol new /home/takahashi/python/qm4d4crystal/examples/pentacene775/input/pentacene_qmmm.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# MM region: Small points (point charge representation)
mol representation Points 3.000000
mol color ColorID 8  ; # white/gray color for point charges
mol selection {beta < 50}
mol material Opaque
mol addrep top

# Alternative: Very small VDW spheres for point charges
mol representation VDW 0.100000 12.000000
mol color ColorID 2  ; # gray
mol selection {beta < 50}
mol material Transparent
mol addrep top

# QM region: CPK balls (normal atomic view)
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol color Name
mol selection {beta > 50}
mol material Opaque
mol addrep top

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on QM region
set qm_sel [atomselect top "beta > 50"]
if {[$qm_sel num] > 0} {
    set center [measure center $qm_sel]
    molinfo top set center [list $center]
}

puts "QM region: CPK representation (colored by element)"
puts "MM region: Point charges (small gray spheres)"
puts ""
puts "Useful commands:"
puts "  mol showrep top 0 off    - Hide MM points"
puts "  mol showrep top 1 off    - Hide MM small spheres"
puts "  mol showrep top 2 off    - Hide QM region"
puts ""
puts "To change point charge size:"
puts "  mol modstyle 0 top Points 5.0     - Larger points"
puts "  mol modstyle 1 top VDW 0.05 12    - Smaller spheres"
