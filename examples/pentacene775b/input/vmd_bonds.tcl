# VMD script to visualize pentacene crystal with 4 bond types
# Run with: vmd -e vmd_bonds.tcl
# Or in VMD: source vmd_bonds.tcl

# Load pentacene structure
mol new pentacene_331_vmd.pdb type pdb waitfor all

# Set molecule representation
mol delrep 0 top
mol representation Licorice 0.15 12.0 12.0
mol color Element
mol selection {all}
mol material Transparent
mol addrep top

# Set background
color Display Background white
display projection Orthographic
axes location Off

# Define colors for bond types
# 1=red, 0=blue, 7=green, 3=orange
set bond_colors [list 1 0 7 3]
set bond_names [list "Type1_AB_herringbone_4.89A" "Type2_AA_a-axis_6.28A" "Type3_AB_diagonal_7.29A" "Type4_AA_b-axis_7.71A"]

# Molecular center coordinates
set mol_centers {
  { 3.5313 4.5130 4.7941 }
  { 3.9814 6.4424 9.2627 }
  { 4.2674 12.1916 4.7941 }
  { 4.7175 14.1210 9.2627 }
  { 5.0035 19.8702 4.7941 }
  { 5.4537 21.7996 9.2627 }
  { 9.8066 4.5130 4.7941 }
  { 10.2567 6.4424 9.2627 }
  { 10.5427 12.1916 4.7941 }
  { 10.9928 14.1210 9.2627 }
  { 11.2788 19.8702 4.7941 }
  { 11.7290 21.7996 9.2627 }
  { 16.0819 4.5130 4.7941 }
  { 16.5320 6.4424 9.2627 }
  { 16.8180 12.1916 4.7941 }
  { 17.2681 14.1210 9.2627 }
  { 17.5541 19.8702 4.7941 }
  { 18.0043 21.7996 9.2627 }
}

set mol_types {
  A
  B
  A
  B
  A
  B
  A
  B
  A
  B
  A
  B
  A
  B
  A
  B
  A
  B
}

# Type 1: A-B herringbone bonds (4.89 A) - RED
graphics top color red
graphics top material Opaque
graphics top cylinder { 3.531 4.513 4.794 } { 3.981 6.442 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.267 12.192 4.794 } { 4.718 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 5.003 19.870 4.794 } { 5.454 21.800 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 9.807 4.513 4.794 } { 10.257 6.442 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.543 12.192 4.794 } { 10.993 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 11.279 19.870 4.794 } { 11.729 21.800 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.082 4.513 4.794 } { 16.532 6.442 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.818 12.192 4.794 } { 17.268 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 17.554 19.870 4.794 } { 18.004 21.800 9.263 } radius 0.2 resolution 20 filled yes

# Type 2: A-A/B-B a-axis bonds (6.28 A) - BLUE
graphics top color blue
graphics top cylinder { 3.531 4.513 4.794 } { 9.807 4.513 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 3.981 6.442 9.263 } { 10.257 6.442 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.267 12.192 4.794 } { 10.543 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.718 14.121 9.263 } { 10.993 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 5.003 19.870 4.794 } { 11.279 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 5.454 21.800 9.263 } { 11.729 21.800 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 9.807 4.513 4.794 } { 16.082 4.513 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.257 6.442 9.263 } { 16.532 6.442 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.543 12.192 4.794 } { 16.818 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.993 14.121 9.263 } { 17.268 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 11.279 19.870 4.794 } { 17.554 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 11.729 21.800 9.263 } { 18.004 21.800 9.263 } radius 0.2 resolution 20 filled yes

# Type 3: A-B diagonal bonds (7.29 A) - GREEN
graphics top color green
graphics top cylinder { 3.981 6.442 9.263 } { 4.267 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.718 14.121 9.263 } { 5.003 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.257 6.442 9.263 } { 10.543 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.993 14.121 9.263 } { 11.279 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.532 6.442 9.263 } { 16.818 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 17.268 14.121 9.263 } { 17.554 19.870 4.794 } radius 0.2 resolution 20 filled yes

# Type 4: A-A/B-B b-axis bonds (7.71 A) - ORANGE
graphics top color orange
graphics top cylinder { 3.531 4.513 4.794 } { 4.267 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 3.981 6.442 9.263 } { 4.718 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.267 12.192 4.794 } { 5.003 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 4.718 14.121 9.263 } { 5.454 21.800 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 9.807 4.513 4.794 } { 10.543 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.257 6.442 9.263 } { 10.993 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.543 12.192 4.794 } { 11.279 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 10.993 14.121 9.263 } { 11.729 21.800 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.082 4.513 4.794 } { 16.818 12.192 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.532 6.442 9.263 } { 17.268 14.121 9.263 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 16.818 12.192 4.794 } { 17.554 19.870 4.794 } radius 0.2 resolution 20 filled yes
graphics top cylinder { 17.268 14.121 9.263 } { 18.004 21.800 9.263 } radius 0.2 resolution 20 filled yes

# Draw molecular centers
# Type A = cyan spheres, Type B = pink spheres
graphics top color cyan
graphics top sphere { 3.531 4.513 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 3.981 6.442 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 4.267 12.192 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 4.718 14.121 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 5.003 19.870 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 5.454 21.800 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 9.807 4.513 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 10.257 6.442 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 10.543 12.192 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 10.993 14.121 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 11.279 19.870 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 11.729 21.800 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 16.082 4.513 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 16.532 6.442 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 16.818 12.192 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 17.268 14.121 9.263 } radius 0.5 resolution 20
graphics top color cyan
graphics top sphere { 17.554 19.870 4.794 } radius 0.5 resolution 20
graphics top color pink
graphics top sphere { 18.004 21.800 9.263 } radius 0.5 resolution 20

# Draw supercell box
graphics top color black
graphics top line { 0.000 0.000 0.000 } { 18.826 0.000 0.000 } width 3 style solid
graphics top line { 0.000 0.000 0.000 } { 2.208 23.036 0.000 } width 3 style solid
graphics top line { 18.826 0.000 0.000 } { 21.034 23.036 0.000 } width 3 style solid
graphics top line { 2.208 23.036 0.000 } { 21.034 23.036 0.000 } width 3 style solid
graphics top line { 0.501 3.277 14.057 } { 19.327 3.277 14.057 } width 3 style solid
graphics top line { 0.501 3.277 14.057 } { 2.710 26.313 14.057 } width 3 style solid
graphics top line { 19.327 3.277 14.057 } { 21.536 26.313 14.057 } width 3 style solid
graphics top line { 2.710 26.313 14.057 } { 21.536 26.313 14.057 } width 3 style solid
graphics top line { 0.000 0.000 0.000 } { 0.501 3.277 14.057 } width 3 style solid
graphics top line { 18.826 0.000 0.000 } { 19.327 3.277 14.057 } width 3 style solid
graphics top line { 2.208 23.036 0.000 } { 2.710 26.313 14.057 } width 3 style solid
graphics top line { 21.034 23.036 0.000 } { 21.536 26.313 14.057 } width 3 style solid

# Add legend
graphics top color black
graphics top text { -3 -3 -3 } "Bond Types:" size 1.2
graphics top color red
graphics top text { -3 -5 -3 } "Red: A-B herringbone (4.89 A)" size 1.0
graphics top color blue
graphics top text { -3 -7 -3 } "Blue: A-A/B-B a-axis (6.28 A)" size 1.0
graphics top color green
graphics top text { -3 -9 -3 } "Green: A-B diagonal (7.29 A)" size 1.0
graphics top color orange
graphics top text { -3 -11 -3 } "Orange: A-A/B-B b-axis (7.71 A)" size 1.0

# Set initial view (top-down along c-axis)
display resetview
rotate x by -90
scale by 1.5

puts "Visualization loaded successfully!"
puts "Bond types:"
puts "  Red:    A-B herringbone (4.89 A)"
puts "  Blue:   A-A/B-B a-axis (6.28 A)"
puts "  Green:  A-B diagonal (7.29 A)"
puts "  Orange: A-A/B-B b-axis (7.71 A)"
