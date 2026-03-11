# VMD visualization script for selected molecules
# Molecules: QM(156), 89, 126, 127, 130, 159, 167, 195, 197, 207, 250
# Usage: vmd -e visualize_selected.tcl

# Load structure
mol new /home/takahashi/python/q4mdkit/examples/pentacene775-single/input/pentacene.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Molecule 89: Licorice (red, ColorID 1)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 1
mol selection {index 3168 3169 3170 3171 3172 3173 3174 3175 3176 3177 3178 3179 3180 3181 3182 3183 3184 3185 3186 3187 3188 3189 3190 3191 3192 3193 3194 3195 3196 3197 3198 3199 3200 3201 3202 3203}
mol material Opaque
mol addrep top

# Molecule 126: Licorice (orange, ColorID 3)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 3
mol selection {index 4500 4501 4502 4503 4504 4505 4506 4507 4508 4509 4510 4511 4512 4513 4514 4515 4516 4517 4518 4519 4520 4521 4522 4523 4524 4525 4526 4527 4528 4529 4530 4531 4532 4533 4534 4535}
mol material Opaque
mol addrep top

# Molecule 127: Licorice (purple, ColorID 11)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 11
mol selection {index 4536 4537 4538 4539 4540 4541 4542 4543 4544 4545 4546 4547 4548 4549 4550 4551 4552 4553 4554 4555 4556 4557 4558 4559 4560 4561 4562 4563 4564 4565 4566 4567 4568 4569 4570 4571}
mol material Opaque
mol addrep top

# Molecule 130: Licorice (green, ColorID 7)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 7
mol selection {index 4644 4645 4646 4647 4648 4649 4650 4651 4652 4653 4654 4655 4656 4657 4658 4659 4660 4661 4662 4663 4664 4665 4666 4667 4668 4669 4670 4671 4672 4673 4674 4675 4676 4677 4678 4679}
mol material Opaque
mol addrep top

# Molecule 159: Licorice (cyan, ColorID 10)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 10
mol selection {index 5688 5689 5690 5691 5692 5693 5694 5695 5696 5697 5698 5699 5700 5701 5702 5703 5704 5705 5706 5707 5708 5709 5710 5711 5712 5713 5714 5715 5716 5717 5718 5719 5720 5721 5722 5723}
mol material Opaque
mol addrep top

# Molecule 167: Licorice (blue, ColorID 0)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 0
mol selection {index 5976 5977 5978 5979 5980 5981 5982 5983 5984 5985 5986 5987 5988 5989 5990 5991 5992 5993 5994 5995 5996 5997 5998 5999 6000 6001 6002 6003 6004 6005 6006 6007 6008 6009 6010 6011}
mol material Opaque
mol addrep top

# Molecule 195: Licorice (yellow, ColorID 4)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 4
mol selection {index 6984 6985 6986 6987 6988 6989 6990 6991 6992 6993 6994 6995 6996 6997 6998 6999 7000 7001 7002 7003 7004 7005 7006 7007 7008 7009 7010 7011 7012 7013 7014 7015 7016 7017 7018 7019}
mol material Opaque
mol addrep top

# Molecule 197: Licorice (pink, ColorID 9)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 9
mol selection {index 7056 7057 7058 7059 7060 7061 7062 7063 7064 7065 7066 7067 7068 7069 7070 7071 7072 7073 7074 7075 7076 7077 7078 7079 7080 7081 7082 7083 7084 7085 7086 7087 7088 7089 7090 7091}
mol material Opaque
mol addrep top

# Molecule 207: Licorice (silver, ColorID 6)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 6
mol selection {index 7416 7417 7418 7419 7420 7421 7422 7423 7424 7425 7426 7427 7428 7429 7430 7431 7432 7433 7434 7435 7436 7437 7438 7439 7440 7441 7442 7443 7444 7445 7446 7447 7448 7449 7450 7451}
mol material Opaque
mol addrep top

# Molecule 250: Licorice (ochre, ColorID 14)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 14
mol selection {index 8964 8965 8966 8967 8968 8969 8970 8971 8972 8973 8974 8975 8976 8977 8978 8979 8980 8981 8982 8983 8984 8985 8986 8987 8988 8989 8990 8991 8992 8993 8994 8995 8996 8997 8998 8999}
mol material Opaque
mol addrep top

# QM molecule 156: CPK representation (element colors to show C-H bonds)
mol representation CPK 1.0 0.3 12.0 12.0
mol color Name
mol selection {index 5580 5581 5582 5583 5584 5585 5586 5587 5588 5589 5590 5591 5592 5593 5594 5595 5596 5597 5598 5599 5600 5601 5602 5603 5604 5605 5606 5607 5608 5609 5610 5611 5612 5613 5614 5615}
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
set sel [atomselect top "index 3168 to 3203 4500 to 4535 4536 to 4571 4644 to 4679 5580 to 5615 5688 to 5723 5976 to 6011 6984 to 7019 7056 to 7091 7416 to 7451 8964 to 8999"]
set center [measure center $sel]
molinfo top set center [list $center]

# ============================================================
# Add molecule labels at center of each molecule
# ============================================================
proc add_mol_label {mol_id first_atom last_atom color_id} {
    set sel [atomselect top "index $first_atom to $last_atom"]
    set center [measure center $sel]
    $sel delete
    graphics top color $color_id
    graphics top text $center "$mol_id" size 0.8 thickness 1.5
}

# Labels with matching colors (disabled)
# add_mol_label 89 3168 3203 1     ;# red
# add_mol_label 126 4500 4535 3    ;# orange
# add_mol_label 127 4536 4571 11   ;# purple
# add_mol_label 130 4644 4679 7    ;# green
# add_mol_label 159 5688 5723 10   ;# cyan
# add_mol_label 167 5976 6011 0    ;# blue
# add_mol_label 195 6984 7019 4    ;# yellow
# add_mol_label 197 7056 7091 9    ;# pink
# add_mol_label 207 7416 7451 6    ;# silver
# add_mol_label 250 8964 8999 14   ;# ochre
# add_mol_label 156 5580 5615 2    ;# gray (QM)

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/pentacene775-single/input/viewpoint_selected.tcl"}} {
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {rotate_matrix} {[molinfo top get rotate_matrix]}"
    puts $fp "molinfo top set {center_matrix} {[molinfo top get center_matrix]}"
    puts $fp "molinfo top set {scale_matrix} {[molinfo top get scale_matrix]}"
    puts $fp "molinfo top set {global_matrix} {[molinfo top get global_matrix]}"
    close $fp
    puts "Viewpoint saved to: $filename"
}

proc load_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/pentacene775-single/input/viewpoint_selected.tcl"}} {
    if {[file exists $filename]} {
        source $filename
        puts "Viewpoint loaded from: $filename"
    } else {
        puts "Viewpoint file not found: $filename"
    }
}

# Auto-load viewpoint if exists
if {[file exists "/home/takahashi/python/q4mdkit/examples/pentacene775-single/input/viewpoint_selected.tcl"]} {
    load_viewpoint
}

puts ""
puts "=== Selected Molecules ==="
puts "  89: red"
puts "  126: orange"
puts "  127: purple"
puts "  130: yellow"
puts "  159: tan"
puts "  167: blue"
puts "  195: ochre"
puts "  197: cyan"
puts "  207: pink"
puts "  250: lime"
puts "  156: gray (QM molecule)"
puts ""
puts "Toggle representations:"
puts "  mol showrep top 0 off/on  - Molecule 89 (red)"
puts "  mol showrep top 1 off/on  - Molecule 126 (orange)"
puts "  mol showrep top 2 off/on  - Molecule 127 (purple)"
puts "  mol showrep top 3 off/on  - Molecule 130 (yellow)"
puts "  mol showrep top 4 off/on  - Molecule 159 (tan)"
puts "  mol showrep top 5 off/on  - Molecule 167 (blue)"
puts "  mol showrep top 6 off/on  - Molecule 195 (ochre)"
puts "  mol showrep top 7 off/on  - Molecule 197 (cyan)"
puts "  mol showrep top 8 off/on  - Molecule 207 (pink)"
puts "  mol showrep top 9 off/on  - Molecule 250 (lime)"
puts "  mol showrep top 10 off/on  - Molecule 156 (QM)"
puts ""
puts "Commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
puts "  graphics top delete all - Hide labels"
