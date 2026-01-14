# VMD script to visualize pentacene crystal structure
# Similar to crystallographic figure with unit cell, axes, and transfer vectors
# Run with: vmd -e vmd_crystal_view.tcl

# Load structure
mol new pentacene_221.xyz type xyz waitfor all

# Set molecule representation - gray tubes
mol delrep 0 top
mol representation Licorice 0.15 20.0 20.0
mol color ColorID 2
mol selection {all}
mol material Opaque
mol addrep top

# Set background
color Display Background white
display projection Orthographic
axes location Off
display depthcue off

# Draw unit cell box (black lines)
graphics top color black
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

# Draw coordinate axes
graphics top color green
graphics top cone { 0.000 0.000 0.000 } { 3.765 0.000 0.000 } radius 0.3 resolution 20
graphics top cylinder { 0.000 0.000 0.000 } { 3.765 0.000 0.000 } radius 0.1 resolution 20 filled yes
graphics top text { 4.393 0.000 0.000 } "a" size 1.5
graphics top cone { 0.000 0.000 0.000 } { 0.442 4.607 0.000 } radius 0.3 resolution 20
graphics top cylinder { 0.000 0.000 0.000 } { 0.442 4.607 0.000 } radius 0.1 resolution 20 filled yes
graphics top text { 0.515 5.375 0.000 } "b" size 1.5
graphics top color blue
graphics top cone { 0.000 0.000 0.000 } { 0.301 1.966 8.434 } radius 0.3 resolution 20
graphics top cylinder { 0.000 0.000 0.000 } { 0.301 1.966 8.434 } radius 0.1 resolution 20 filled yes
graphics top text { 0.351 2.294 9.840 } "c" size 1.5
graphics top color black
graphics top text { -0.548 -0.498 0.000 } "O" size 1.5

# Draw transfer vectors t1 and t2 (red)
graphics top color red
graphics top cylinder { 3.531 4.513 4.794 } { 3.891 6.057 8.369 } radius 0.15 resolution 20 filled yes
graphics top cone { 3.891 6.057 8.369 } { 3.965 6.372 9.100 } radius 0.4 resolution 20
graphics top text { 4.030 6.649 9.740 } "t1" size 1.5
graphics top cylinder { 3.531 4.513 4.794 } { -0.234 4.513 4.794 } radius 0.15 resolution 20 filled yes
graphics top cone { -0.234 4.513 4.794 } { -1.034 4.513 4.794 } radius 0.4 resolution 20
graphics top text { -1.734 4.513 4.794 } "t2" size 1.5

# Set view - similar to crystallographic view
display resetview
rotate x by -70
rotate y by 20
rotate z by -30
scale by 1.2

# Center on the structure
display resize 800 600

puts "Crystal structure visualization loaded!"
puts "Axes: a (green), b (green), c (blue)"
puts "Transfer vectors: t1, t2 (red)"
