`# VMD visualization script for QM molecule and adjacent layer below
# Molecules: QM(76), 150, 86, 112, 114 (layer z ≈ 25)
# Sorted by minimum atom-atom distance
# Usage: vmd -e visualize_qm_below.tcl

# Load structure
mol new /home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/ph-btbt-c10.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Molecule 150 (min_dist=2.58Å): Licorice (red, ColorID 1)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 1
mol selection {index 9600 9601 9602 9603 9604 9605 9606 9607 9608 9609 9610 9611 9612 9613 9614 9615 9616 9617 9618 9619 9620 9621 9622 9623 9624 9625 9626 9627 9628 9629 9630 9631 9632 9633 9634 9635 9636 9637 9638 9639 9640 9641 9642 9643 9644 9645 9646 9647 9648 9649 9650 9651 9652 9653 9654 9655 9656 9657 9658 9659 9660 9661 9662 9663}
mol material Opaque
mol addrep top

# Molecule 86 (min_dist=2.72Å): Licorice (orange, ColorID 3)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 3
mol selection {index 5504 5505 5506 5507 5508 5509 5510 5511 5512 5513 5514 5515 5516 5517 5518 5519 5520 5521 5522 5523 5524 5525 5526 5527 5528 5529 5530 5531 5532 5533 5534 5535 5536 5537 5538 5539 5540 5541 5542 5543 5544 5545 5546 5547 5548 5549 5550 5551 5552 5553 5554 5555 5556 5557 5558 5559 5560 5561 5562 5563 5564 5565 5566 5567}
mol material Opaque
mol addrep top

# Molecule 112 (min_dist=3.04Å): Licorice (green, ColorID 7)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 7
mol selection {index 7168 7169 7170 7171 7172 7173 7174 7175 7176 7177 7178 7179 7180 7181 7182 7183 7184 7185 7186 7187 7188 7189 7190 7191 7192 7193 7194 7195 7196 7197 7198 7199 7200 7201 7202 7203 7204 7205 7206 7207 7208 7209 7210 7211 7212 7213 7214 7215 7216 7217 7218 7219 7220 7221 7222 7223 7224 7225 7226 7227 7228 7229 7230 7231}
mol material Opaque
mol addrep top

# Molecule 114 (min_dist=3.05Å): Licorice (cyan, ColorID 10)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 10
mol selection {index 7296 7297 7298 7299 7300 7301 7302 7303 7304 7305 7306 7307 7308 7309 7310 7311 7312 7313 7314 7315 7316 7317 7318 7319 7320 7321 7322 7323 7324 7325 7326 7327 7328 7329 7330 7331 7332 7333 7334 7335 7336 7337 7338 7339 7340 7341 7342 7343 7344 7345 7346 7347 7348 7349 7350 7351 7352 7353 7354 7355 7356 7357 7358 7359}
mol material Opaque
mol addrep top

# QM molecule 76: CPK representation (element colors)
mol representation CPK 1.0 0.3 12.0 12.0
mol color Name
mol selection {index 4864 4865 4866 4867 4868 4869 4870 4871 4872 4873 4874 4875 4876 4877 4878 4879 4880 4881 4882 4883 4884 4885 4886 4887 4888 4889 4890 4891 4892 4893 4894 4895 4896 4897 4898 4899 4900 4901 4902 4903 4904 4905 4906 4907 4908 4909 4910 4911 4912 4913 4914 4915 4916 4917 4918 4919 4920 4921 4922 4923 4924 4925 4926 4927}
mol material Opaque
mol addrep top

# Change carbon color to gray for QM region
color Name C gray

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on all selected molecules
set sel [atomselect top "index 4864 to 4927 5504 to 5567 7168 to 7231 7296 to 7359 9600 to 9663"]
set center [measure center $sel]
molinfo top set center [list $center]

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_below.tcl"}} {
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {rotate_matrix} {[molinfo top get rotate_matrix]}"
    puts $fp "molinfo top set {center_matrix} {[molinfo top get center_matrix]}"
    puts $fp "molinfo top set {scale_matrix} {[molinfo top get scale_matrix]}"
    puts $fp "molinfo top set {global_matrix} {[molinfo top get global_matrix]}"
    close $fp
    puts "Viewpoint saved to: $filename"
}

proc load_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_below.tcl"}} {
    if {[file exists $filename]} {
        source $filename
        puts "Viewpoint loaded from: $filename"
    } else {
        puts "Viewpoint file not found: $filename"
    }
}

# Auto-load viewpoint if exists
if {[file exists "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_below.tcl"]} {
    load_viewpoint
}

puts ""
puts "=== QM + Adjacent Layer Below (z ≈ 25) ==="
puts "  Sorted by minimum atom-atom distance"
puts "  150: red (min_dist=2.58Å)"
puts "  86: orange (min_dist=2.72Å)"
puts "  112: green (min_dist=3.04Å)"
puts "  114: cyan (min_dist=3.05Å)"
puts "  76: gray (QM molecule)"
puts ""
puts "Toggle representations:"
puts "  mol showrep top 0 off/on  - Molecule 150 (red)"
puts "  mol showrep top 1 off/on  - Molecule 86 (orange)"
puts "  mol showrep top 2 off/on  - Molecule 112 (green)"
puts "  mol showrep top 3 off/on  - Molecule 114 (cyan)"
puts "  mol showrep top 4 off/on  - Molecule 76 (QM)"
puts ""
puts "Commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
