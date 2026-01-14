# VMD script to visualize pentacene with 4 bond types
# Run with: vmd -e vmd_bonds_unitcell.tcl

# Load structure
mol new pentacene_bonds_viz.xyz type xyz waitfor all

# Set molecule representation - all atoms as licorice
mol delrep 0 top
mol representation Licorice 0.12 12.0 12.0
mol color Element
mol selection {all}
mol material Opaque
mol addrep top

# Set background
color Display Background white
display projection Orthographic
axes location Off

# Draw molecular centers for central unit cell
graphics top color cyan
graphics top sphere { 3.531 4.513 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 3.981 6.442 9.263 } radius 0.5 resolution 20

# Type 1: A-B herringbone bond (4.89 A) - RED
graphics top color red
graphics top material Opaque
graphics top cylinder { 3.531 4.513 4.794 } { 3.981 6.442 9.263 } radius 0.2 resolution 20 filled yes

# Type 2: A-A a-axis bond (6.28 A) - BLUE
graphics top color blue
graphics top cylinder { 3.531 4.513 4.794 } { -2.744 4.513 4.794 } radius 0.2 resolution 20 filled yes

# Type 3: A-B diagonal bond (7.29 A) - GREEN
graphics top color green
graphics top cylinder { 3.531 4.513 4.794 } { 3.245 -1.236 9.263 } radius 0.2 resolution 20 filled yes

# Type 4: A-A b-axis bond (7.71 A) - ORANGE
graphics top color orange
graphics top cylinder { 3.531 4.513 4.794 } { 2.795 -3.166 4.794 } radius 0.2 resolution 20 filled yes

# Draw unit cell box
graphics top color gray
graphics top line { 0.000 0.000 0.000 } { 6.275 0.000 0.000 } width 2 style solid
graphics top line { 0.000 0.000 0.000 } { 0.736 7.679 0.000 } width 2 style solid
graphics top line { 6.275 0.000 0.000 } { 7.011 7.679 0.000 } width 2 style solid
graphics top line { 0.736 7.679 0.000 } { 7.011 7.679 0.000 } width 2 style solid
graphics top line { 0.501 3.277 14.057 } { 6.777 3.277 14.057 } width 2 style solid
graphics top line { 0.501 3.277 14.057 } { 1.237 10.955 14.057 } width 2 style solid
graphics top line { 6.777 3.277 14.057 } { 7.513 10.955 14.057 } width 2 style solid
graphics top line { 1.237 10.955 14.057 } { 7.513 10.955 14.057 } width 2 style solid
graphics top line { 0.000 0.000 0.000 } { 0.501 3.277 14.057 } width 2 style solid
graphics top line { 6.275 0.000 0.000 } { 6.777 3.277 14.057 } width 2 style solid
graphics top line { 0.736 7.679 0.000 } { 1.237 10.955 14.057 } width 2 style solid
graphics top line { 7.011 7.679 0.000 } { 7.513 10.955 14.057 } width 2 style solid

# Add legend
graphics top color black
graphics top text { -8 -2 0 } "Bond Types:" size 1.2
graphics top color red
graphics top text { -8 -4 0 } "Red: A-B herringbone (4.89 A)" size 1.0
graphics top color blue
graphics top text { -8 -6 0 } "Blue: A-A a-axis (6.28 A)" size 1.0
graphics top color green
graphics top text { -8 -8 0 } "Green: A-B diagonal (7.29 A)" size 1.0
graphics top color orange
graphics top text { -8 -10 0 } "Orange: A-A b-axis (7.71 A)" size 1.0

# Set initial view (top-down along c-axis to see ab-plane)
display resetview
rotate x by -90
scale by 1.5

puts "Visualization loaded successfully!"
puts "Bond types (one example each):"
puts "  Red:    A-B herringbone (4.89 A)"
puts "  Blue:   A-A a-axis (6.28 A)"
puts "  Green:  A-B diagonal (7.29 A)"
puts "  Orange: A-A b-axis (7.71 A)"
