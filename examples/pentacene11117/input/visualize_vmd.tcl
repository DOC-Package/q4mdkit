# VMD visualization script for QM/MM system
# Usage: vmd -e /home/takahashi/python/qm4dcrystal/examples/pentacene71111/input/visualize_vmd.tcl

# Load structure
mol new /home/takahashi/python/qm4dcrystal/examples/pentacene71111/input/pentacene_qmmm.pdb type pdb waitfor all

# Color by B-factor (QM/MM region)
mol delrep 0 top
mol representation Lines 1.000000
mol color Beta
mol selection {beta < 50}
mol material Opaque
mol addrep top

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
set center [measure center $qm_sel]
set box_size 30
molinfo top set center [list $center]

puts "QM region: CPK representation (colored by element)"
puts "MM region: Lines (colored by B-factor)"
puts ""
puts "Useful commands:"
puts "  mol showrep top 0 off    - Hide MM region"
puts "  mol showrep top 1 off    - Hide QM region"
