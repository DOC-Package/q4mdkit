# VMD visualization script for selected molecules
# Molecules: QM(77), 57, 106, 58, 17, 137, 112, 114, 86
# Usage: vmd -e visualize_selected.tcl

# Load structure
mol new /home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/ph-btbt-c10.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Molecule 57: Licorice (red, ColorID 1)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 1
mol selection {index 3584 3585 3586 3587 3588 3589 3590 3591 3592 3593 3594 3595 3596 3597 3598 3599 3600 3601 3602 3603 3604 3605 3606 3607 3608 3609 3610 3611 3612 3613 3614 3615 3616 3617 3618 3619 3620 3621 3622 3623 3624 3625 3626 3627 3628 3629 3630 3631 3632 3633 3634 3635 3636 3637 3638 3639 3640 3641 3642 3643 3644 3645 3646 3647}
mol material Opaque
mol addrep top

# Molecule 106: Licorice (blue, ColorID 0)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 0
mol selection {index 6720 6721 6722 6723 6724 6725 6726 6727 6728 6729 6730 6731 6732 6733 6734 6735 6736 6737 6738 6739 6740 6741 6742 6743 6744 6745 6746 6747 6748 6749 6750 6751 6752 6753 6754 6755 6756 6757 6758 6759 6760 6761 6762 6763 6764 6765 6766 6767 6768 6769 6770 6771 6772 6773 6774 6775 6776 6777 6778 6779 6780 6781 6782 6783}
mol material Opaque
mol addrep top

# Molecule 58: Licorice (orange, ColorID 3)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 3
mol selection {index 3648 3649 3650 3651 3652 3653 3654 3655 3656 3657 3658 3659 3660 3661 3662 3663 3664 3665 3666 3667 3668 3669 3670 3671 3672 3673 3674 3675 3676 3677 3678 3679 3680 3681 3682 3683 3684 3685 3686 3687 3688 3689 3690 3691 3692 3693 3694 3695 3696 3697 3698 3699 3700 3701 3702 3703 3704 3705 3706 3707 3708 3709 3710 3711}
mol material Opaque
mol addrep top

# Molecule 17: Licorice (purple, ColorID 11)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 11
mol selection {index 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 1035 1036 1037 1038 1039 1040 1041 1042 1043 1044 1045 1046 1047 1048 1049 1050 1051 1052 1053 1054 1055 1056 1057 1058 1059 1060 1061 1062 1063 1064 1065 1066 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 1077 1078 1079 1080 1081 1082 1083 1084 1085 1086 1087}
mol material Opaque
mol addrep top

# Molecule 137: Licorice (green, ColorID 7)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 7
mol selection {index 8704 8705 8706 8707 8708 8709 8710 8711 8712 8713 8714 8715 8716 8717 8718 8719 8720 8721 8722 8723 8724 8725 8726 8727 8728 8729 8730 8731 8732 8733 8734 8735 8736 8737 8738 8739 8740 8741 8742 8743 8744 8745 8746 8747 8748 8749 8750 8751 8752 8753 8754 8755 8756 8757 8758 8759 8760 8761 8762 8763 8764 8765 8766 8767}
mol material Opaque
mol addrep top

# Molecule 112: Licorice (yellow, ColorID 4)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 4
mol selection {index 7168 7169 7170 7171 7172 7173 7174 7175 7176 7177 7178 7179 7180 7181 7182 7183 7184 7185 7186 7187 7188 7189 7190 7191 7192 7193 7194 7195 7196 7197 7198 7199 7200 7201 7202 7203 7204 7205 7206 7207 7208 7209 7210 7211 7212 7213 7214 7215 7216 7217 7218 7219 7220 7221 7222 7223 7224 7225 7226 7227 7228 7229 7230 7231}
mol material Opaque
mol addrep top

# Molecule 114: Licorice (cyan, ColorID 10)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 10
mol selection {index 7296 7297 7298 7299 7300 7301 7302 7303 7304 7305 7306 7307 7308 7309 7310 7311 7312 7313 7314 7315 7316 7317 7318 7319 7320 7321 7322 7323 7324 7325 7326 7327 7328 7329 7330 7331 7332 7333 7334 7335 7336 7337 7338 7339 7340 7341 7342 7343 7344 7345 7346 7347 7348 7349 7350 7351 7352 7353 7354 7355 7356 7357 7358 7359}
mol material Opaque
mol addrep top

# Molecule 86: Licorice (pink, ColorID 9)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 9
mol selection {index 5504 5505 5506 5507 5508 5509 5510 5511 5512 5513 5514 5515 5516 5517 5518 5519 5520 5521 5522 5523 5524 5525 5526 5527 5528 5529 5530 5531 5532 5533 5534 5535 5536 5537 5538 5539 5540 5541 5542 5543 5544 5545 5546 5547 5548 5549 5550 5551 5552 5553 5554 5555 5556 5557 5558 5559 5560 5561 5562 5563 5564 5565 5566 5567}
mol material Opaque
mol addrep top

# QM molecule 77: CPK representation (element colors to show C-H bonds)
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
set sel [atomselect top "index 1024 to 1087 3584 to 3647 3648 to 3711 4864 to 4927 5504 to 5567 6720 to 6783 7168 to 7231 7296 to 7359 8704 to 8767"]
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
# add_mol_label 57 3584 3647 1     ;# red
# add_mol_label 106 6720 6783 0    ;# blue
# add_mol_label 58 3648 3711 3     ;# orange
# add_mol_label 17 1024 1087 11    ;# purple
# add_mol_label 137 8704 8767 7    ;# green
# add_mol_label 112 7168 7231 4    ;# yellow
# add_mol_label 114 7296 7359 10   ;# cyan
# add_mol_label 86 5504 5567 9     ;# pink
# add_mol_label 77 4864 4927 2     ;# gray (QM)

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_selected.tcl"}} {
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {rotate_matrix} {[molinfo top get rotate_matrix]}"
    puts $fp "molinfo top set {center_matrix} {[molinfo top get center_matrix]}"
    puts $fp "molinfo top set {scale_matrix} {[molinfo top get scale_matrix]}"
    puts $fp "molinfo top set {global_matrix} {[molinfo top get global_matrix]}"
    close $fp
    puts "Viewpoint saved to: $filename"
}

proc load_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_selected.tcl"}} {
    if {[file exists $filename]} {
        source $filename
        puts "Viewpoint loaded from: $filename"
    } else {
        puts "Viewpoint file not found: $filename"
    }
}

# Auto-load viewpoint if exists
if {[file exists "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_selected.tcl"]} {
    load_viewpoint
}

puts ""
puts "=== Selected Molecules ==="
puts "  57: red"
puts "  106: blue"
puts "  58: orange"
puts "  17: purple"
puts "  137: green"
puts "  112: yellow"
puts "  114: cyan"
puts "  86: pink"
puts "  77: gray (QM molecule)"
puts ""
puts "Toggle representations:"
puts "  mol showrep top 0 off/on  - Molecule 57 (red)"
puts "  mol showrep top 1 off/on  - Molecule 106 (blue)"
puts "  mol showrep top 2 off/on  - Molecule 58 (orange)"
puts "  mol showrep top 3 off/on  - Molecule 17 (purple)"
puts "  mol showrep top 4 off/on  - Molecule 137 (green)"
puts "  mol showrep top 5 off/on  - Molecule 112 (yellow)"
puts "  mol showrep top 6 off/on  - Molecule 114 (cyan)"
puts "  mol showrep top 7 off/on  - Molecule 86 (pink)"
puts "  mol showrep top 8 off/on  - Molecule 77 (QM)"
puts ""
puts "Commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
puts "  graphics top delete all - Hide labels"
