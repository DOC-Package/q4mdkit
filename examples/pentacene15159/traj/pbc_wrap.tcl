# VMD script for wrapping molecules into PBC box
# Usage: In VMD, source this script: source pbc_wrap.tcl

# Load the PBC package
package require pbc

# Wrap all atoms into the primary unit cell
# -compound res: wrap whole residues (molecules) together
# -all: process all frames
pbc wrap -compound res -all

# Draw the unit cell box
pbc box -color blue -width 2

# Alternative: if you want to center the box around the center of mass
# pbc wrap -compound res -center com -centersel "all" -all

puts "PBC wrapping completed. Box is now displayed."
