# VMD script to visualize central pairs of each bond type in pentacene
# Only selected molecules are shown (others excluded)
# Run with: vmd -e vmd_bond_types.tcl

# Load only selected molecules
mol new vmd_bond_types_selected.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Set background
color Display Background white
display projection Orthographic
axes location Off
display depthcue off


# Type 1: herringbone (4.73 Å)
# Molecules 156 and 207
mol representation Licorice 0.15 12.0 12.0
mol color ColorID 1
mol selection {resid 156 207}
mol material Opaque
mol addrep top

# Type 2: second nearest (5.20 Å)
# Molecules 126 and 156
mol representation Licorice 0.15 12.0 12.0
mol color ColorID 10
mol selection {resid 126 156}
mol material Opaque
mol addrep top

# Type 3: a-axis (6.28 Å)
# Molecules 156 and 250
mol representation Licorice 0.15 12.0 12.0
mol color ColorID 0
mol selection {resid 156 250}
mol material Opaque
mol addrep top

# Type 4: b-axis (7.71 Å)
# Molecules 151 and 156
mol representation Licorice 0.15 12.0 12.0
mol color ColorID 3
mol selection {resid 151 156}
mol material Opaque
mol addrep top

# Draw bond vectors between molecule centers

# Type 1: herringbone
graphics top color red
graphics top cylinder { 24.478 37.850 34.138 } { 27.247 34.010 34.138 } radius 0.25 resolution 20 filled yes
graphics top text { 25.862 37.430 34.138 } "T1: 4.73A" size 1.2

# Type 2: second nearest
graphics top color cyan
graphics top cylinder { 20.972 34.010 34.138 } { 24.478 37.850 34.138 } radius 0.25 resolution 20 filled yes
graphics top text { 22.725 37.430 34.138 } "T2: 5.20A" size 1.2

# Type 3: a-axis
graphics top color blue
graphics top cylinder { 24.478 37.850 34.138 } { 30.753 37.850 34.138 } radius 0.25 resolution 20 filled yes
graphics top text { 27.615 39.350 34.138 } "T3: 6.28A" size 1.2

# Type 4: b-axis
graphics top color orange
graphics top cylinder { 23.742 30.171 34.138 } { 24.478 37.850 34.138 } radius 0.25 resolution 20 filled yes
graphics top text { 24.110 35.510 34.138 } "T4: 7.71A" size 1.2

# Draw molecule center markers
graphics top color red
graphics top sphere { 24.478 37.850 34.138 } radius 0.25 resolution 20
graphics top sphere { 27.247 34.010 34.138 } radius 0.25 resolution 20
graphics top color cyan
graphics top sphere { 20.972 34.010 34.138 } radius 0.25 resolution 20
graphics top sphere { 24.478 37.850 34.138 } radius 0.25 resolution 20
graphics top color blue
graphics top sphere { 24.478 37.850 34.138 } radius 0.25 resolution 20
graphics top sphere { 30.753 37.850 34.138 } radius 0.25 resolution 20
graphics top color orange
graphics top sphere { 23.742 30.171 34.138 } radius 0.25 resolution 20
graphics top sphere { 24.478 37.850 34.138 } radius 0.25 resolution 20

# Bounding box aligned with crystal axes
graphics top color gray
graphics top line { 18.350 26.884 29.921 } { 31.527 26.884 29.921 } width 1 style dashed
graphics top line { 18.350 26.884 29.921 } { 19.528 39.170 29.921 } width 1 style dashed
graphics top line { 31.527 26.884 29.921 } { 32.706 39.170 29.921 } width 1 style dashed
graphics top line { 19.528 39.170 29.921 } { 32.706 39.170 29.921 } width 1 style dashed
graphics top line { 18.651 28.851 38.355 } { 31.828 28.851 38.355 } width 1 style dashed
graphics top line { 18.651 28.851 38.355 } { 19.829 41.136 38.355 } width 1 style dashed
graphics top line { 31.828 28.851 38.355 } { 33.007 41.136 38.355 } width 1 style dashed
graphics top line { 19.829 41.136 38.355 } { 33.007 41.136 38.355 } width 1 style dashed
graphics top line { 18.350 26.884 29.921 } { 18.651 28.851 38.355 } width 1 style dashed
graphics top line { 31.527 26.884 29.921 } { 31.828 28.851 38.355 } width 1 style dashed
graphics top line { 19.528 39.170 29.921 } { 19.829 41.136 38.355 } width 1 style dashed
graphics top line { 32.706 39.170 29.921 } { 33.007 41.136 38.355 } width 1 style dashed

# Unit cell vectors (a, b, c)
graphics top color black
graphics top text { 15.4 23.8 29.9 } "O" size 1.0
graphics top color green
graphics top cylinder { 16.947 25.349 29.921 } { 20.085 25.349 29.921 } radius 0.12 resolution 20 filled yes
graphics top cone { 20.085 25.349 29.921 } { 20.685 25.349 29.921 } radius 0.3 resolution 20
graphics top text { 21.185 25.349 29.921 } "a" size 1.2
graphics top cylinder { 16.947 25.349 29.921 } { 17.316 29.188 29.921 } radius 0.12 resolution 20 filled yes
graphics top cone { 17.316 29.188 29.921 } { 17.373 29.785 29.921 } radius 0.3 resolution 20
graphics top text { 17.373 30.285 29.921 } "b" size 1.2
graphics top color blue
graphics top cylinder { 16.947 25.349 29.921 } { 17.198 26.987 36.949 } radius 0.12 resolution 20 filled yes
graphics top cone { 17.198 26.987 36.949 } { 17.219 27.123 37.533 } radius 0.3 resolution 20
graphics top text { 17.219 27.123 38.033 } "c" size 1.2

# Legend
graphics top color black
graphics top text { 35.1 20.9 34.1 } "Bond Types:" size 1.2
graphics top color red
graphics top text { 35.1 18.4 34.1 } "Type 1: herringbone (4.73 A)" size 1.0
graphics top color cyan
graphics top text { 35.1 15.9 34.1 } "Type 2: second nearest (5.20 A)" size 1.0
graphics top color blue
graphics top text { 35.1 13.4 34.1 } "Type 3: a-axis (6.28 A)" size 1.0
graphics top color orange
graphics top text { 35.1 10.9 34.1 } "Type 4: b-axis (7.71 A)" size 1.0

# Set view centered on selected molecules
display resetview
rotate x by -60
rotate y by 10
scale by 1.2

puts "Visualization loaded!"
puts "Only selected molecules are shown."
puts "Bond types:"
puts "  Type 1 (herringbone): Mol 156-207, 4.73 A"
puts "  Type 2 (second nearest): Mol 126-156, 5.20 A"
puts "  Type 3 (a-axis): Mol 156-250, 6.28 A"
puts "  Type 4 (b-axis): Mol 151-156, 7.71 A"
