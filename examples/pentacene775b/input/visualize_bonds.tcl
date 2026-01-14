# VMD script to visualize pentacene crystal bonds
# Load molecule structures
mol new pentacene_331_viz.pdb type pdb
mol new mol_centers.xyz type xyz

# Style for pentacene molecules
mol modstyle 0 0 Licorice 0.1 12.0 12.0
mol modcolor 0 0 Element
mol modmaterial 0 0 Transparent

# Style for molecular centers
mol modstyle 0 1 VDW 0.5 12.0
mol modcolor 0 1 Name

# Define colors for bond types
# Type 1 (herringbone A-B): Red
# Type 2 (a-axis A-A/B-B): Blue  
# Type 3 (diagonal A-B): Green
# Type 4 (b-axis A-A/B-B): Orange

# Molecular center coordinates
set mol_centers {
  {3.5313 4.5130 4.7941}
  {3.9814 6.4424 9.2627}
  {4.2674 12.1916 4.7941}
  {4.7175 14.1210 9.2627}
  {5.0035 19.8702 4.7941}
  {5.4537 21.7996 9.2627}
  {9.8066 4.5130 4.7941}
  {10.2567 6.4424 9.2627}
  {10.5427 12.1916 4.7941}
  {10.9928 14.1210 9.2627}
  {11.2788 19.8702 4.7941}
  {11.7290 21.7996 9.2627}
  {16.0819 4.5130 4.7941}
  {16.5320 6.4424 9.2627}
  {16.8180 12.1916 4.7941}
  {17.2681 14.1210 9.2627}
  {17.5541 19.8702 4.7941}
  {18.0043 21.7996 9.2627}
}

# Draw Type 1 bonds (herringbone A-B, red)
draw color red
draw material Opaque
draw cylinder {3.531 4.513 4.794} {3.981 6.442 9.263} radius 0.15
draw cylinder {4.267 12.192 4.794} {4.718 14.121 9.263} radius 0.15
draw cylinder {5.003 19.870 4.794} {5.454 21.800 9.263} radius 0.15
draw cylinder {9.807 4.513 4.794} {10.257 6.442 9.263} radius 0.15
draw cylinder {10.543 12.192 4.794} {10.993 14.121 9.263} radius 0.15
draw cylinder {11.279 19.870 4.794} {11.729 21.800 9.263} radius 0.15
draw cylinder {16.082 4.513 4.794} {16.532 6.442 9.263} radius 0.15
draw cylinder {16.818 12.192 4.794} {17.268 14.121 9.263} radius 0.15
draw cylinder {17.554 19.870 4.794} {18.004 21.800 9.263} radius 0.15

# Draw Type 2 bonds (a-axis A-A/B-B, blue)
draw color blue
draw cylinder {3.531 4.513 4.794} {9.807 4.513 4.794} radius 0.15
draw cylinder {3.981 6.442 9.263} {10.257 6.442 9.263} radius 0.15
draw cylinder {4.267 12.192 4.794} {10.543 12.192 4.794} radius 0.15
draw cylinder {4.718 14.121 9.263} {10.993 14.121 9.263} radius 0.15
draw cylinder {5.003 19.870 4.794} {11.279 19.870 4.794} radius 0.15
draw cylinder {5.454 21.800 9.263} {11.729 21.800 9.263} radius 0.15
draw cylinder {9.807 4.513 4.794} {16.082 4.513 4.794} radius 0.15
draw cylinder {10.257 6.442 9.263} {16.532 6.442 9.263} radius 0.15
draw cylinder {10.543 12.192 4.794} {16.818 12.192 4.794} radius 0.15
draw cylinder {10.993 14.121 9.263} {17.268 14.121 9.263} radius 0.15
draw cylinder {11.279 19.870 4.794} {17.554 19.870 4.794} radius 0.15
draw cylinder {11.729 21.800 9.263} {18.004 21.800 9.263} radius 0.15

# Draw Type 3 bonds (diagonal A-B, green)
draw color green
draw cylinder {3.981 6.442 9.263} {4.267 12.192 4.794} radius 0.15
draw cylinder {4.718 14.121 9.263} {5.003 19.870 4.794} radius 0.15
draw cylinder {10.257 6.442 9.263} {10.543 12.192 4.794} radius 0.15
draw cylinder {10.993 14.121 9.263} {11.279 19.870 4.794} radius 0.15
draw cylinder {16.532 6.442 9.263} {16.818 12.192 4.794} radius 0.15
draw cylinder {17.268 14.121 9.263} {17.554 19.870 4.794} radius 0.15

# Draw Type 4 bonds (b-axis A-A/B-B, orange)
draw color orange
draw cylinder {3.531 4.513 4.794} {4.267 12.192 4.794} radius 0.15
draw cylinder {3.981 6.442 9.263} {4.718 14.121 9.263} radius 0.15
draw cylinder {4.267 12.192 4.794} {5.003 19.870 4.794} radius 0.15
draw cylinder {4.718 14.121 9.263} {5.454 21.800 9.263} radius 0.15
draw cylinder {9.807 4.513 4.794} {10.543 12.192 4.794} radius 0.15
draw cylinder {10.257 6.442 9.263} {10.993 14.121 9.263} radius 0.15
draw cylinder {10.543 12.192 4.794} {11.279 19.870 4.794} radius 0.15
draw cylinder {10.993 14.121 9.263} {11.729 21.800 9.263} radius 0.15
draw cylinder {16.082 4.513 4.794} {16.818 12.192 4.794} radius 0.15
draw cylinder {16.532 6.442 9.263} {17.268 14.121 9.263} radius 0.15
draw cylinder {16.818 12.192 4.794} {17.554 19.870 4.794} radius 0.15
draw cylinder {17.268 14.121 9.263} {18.004 21.800 9.263} radius 0.15

# Draw unit cell box
draw color white
draw line {0.000 0.000 0.000} {18.826 0.000 0.000} width 2
draw line {0.000 0.000 0.000} {2.208 23.036 0.000} width 2
draw line {18.826 0.000 0.000} {21.034 23.036 0.000} width 2
draw line {2.208 23.036 0.000} {21.034 23.036 0.000} width 2
draw line {0.501 3.277 14.057} {19.327 3.277 14.057} width 2
draw line {0.501 3.277 14.057} {2.710 26.313 14.057} width 2
draw line {19.327 3.277 14.057} {21.536 26.313 14.057} width 2
draw line {2.710 26.313 14.057} {21.536 26.313 14.057} width 2
draw line {0.000 0.000 0.000} {0.501 3.277 14.057} width 2
draw line {18.826 0.000 0.000} {19.327 3.277 14.057} width 2
draw line {2.208 23.036 0.000} {2.710 26.313 14.057} width 2
draw line {21.034 23.036 0.000} {21.536 26.313 14.057} width 2

# Add labels
draw color white
draw text {-2 -2 -2} "Red: A-B herringbone (4.89 A)" size 1
draw text {-2 -3 -2} "Blue: A-A/B-B a-axis (6.28 A)" size 1
draw text {-2 -4 -2} "Green: A-B diagonal (7.29 A)" size 1
draw text {-2 -5 -2} "Orange: A-A/B-B b-axis (7.71 A)" size 1

# Set view
display projection Orthographic
axes location Off
