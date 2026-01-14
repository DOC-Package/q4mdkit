# VMD visualization script - All QM-like view
# Shows entire system with CPK representation (QM-style)
# Usage: vmd -e /home/takahashi/python/qm4d4crystal/examples/pentacene775/input/visualize_vmd_allqm.tcl

# Load structure
mol new /home/takahashi/python/qm4d4crystal/examples/pentacene775/input/pentacene_qmmm.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# All atoms: CPK balls (QM-like style)
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol color Name
mol selection {all}
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

puts "All atoms shown with CPK representation (QM-like style)"
puts "Colored by element name"
