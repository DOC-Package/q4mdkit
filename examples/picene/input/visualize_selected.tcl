# VMD visualization script for selected molecules
# Molecules: QM(202), 128, 166, 169, 200, 204, 228, 232, 260
# Usage: vmd -e visualize_selected.tcl

# Load structure
mol new /home/takahashi/python/q4mdkit/examples/picene/input/picene.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Molecule 128: Licorice (red, ColorID 1)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 1
mol selection {index 4572 4573 4574 4575 4576 4577 4578 4579 4580 4581 4582 4583 4584 4585 4586 4587 4588 4589 4590 4591 4592 4593 4594 4595 4596 4597 4598 4599 4600 4601 4602 4603 4604 4605 4606 4607}
mol material Opaque
mol addrep top

# Molecule 166: Licorice (orange, ColorID 3)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 3
mol selection {index 5940 5941 5942 5943 5944 5945 5946 5947 5948 5949 5950 5951 5952 5953 5954 5955 5956 5957 5958 5959 5960 5961 5962 5963 5964 5965 5966 5967 5968 5969 5970 5971 5972 5973 5974 5975}
mol material Opaque
mol addrep top

# Molecule 169: Licorice (purple, ColorID 11)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 11
mol selection {index 6048 6049 6050 6051 6052 6053 6054 6055 6056 6057 6058 6059 6060 6061 6062 6063 6064 6065 6066 6067 6068 6069 6070 6071 6072 6073 6074 6075 6076 6077 6078 6079 6080 6081 6082 6083}
mol material Opaque
mol addrep top

# Molecule 200: Licorice (green, ColorID 7)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 7
mol selection {index 7164 7165 7166 7167 7168 7169 7170 7171 7172 7173 7174 7175 7176 7177 7178 7179 7180 7181 7182 7183 7184 7185 7186 7187 7188 7189 7190 7191 7192 7193 7194 7195 7196 7197 7198 7199}
mol material Opaque
mol addrep top

# Molecule 204: Licorice (yellow, ColorID 4)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 4
mol selection {index 7308 7309 7310 7311 7312 7313 7314 7315 7316 7317 7318 7319 7320 7321 7322 7323 7324 7325 7326 7327 7328 7329 7330 7331 7332 7333 7334 7335 7336 7337 7338 7339 7340 7341 7342 7343}
mol material Opaque
mol addrep top

# Molecule 228: Licorice (cyan, ColorID 10)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 10
mol selection {index 8172 8173 8174 8175 8176 8177 8178 8179 8180 8181 8182 8183 8184 8185 8186 8187 8188 8189 8190 8191 8192 8193 8194 8195 8196 8197 8198 8199 8200 8201 8202 8203 8204 8205 8206 8207}
mol material Opaque
mol addrep top

# Molecule 232: Licorice (pink, ColorID 9)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 9
mol selection {index 8316 8317 8318 8319 8320 8321 8322 8323 8324 8325 8326 8327 8328 8329 8330 8331 8332 8333 8334 8335 8336 8337 8338 8339 8340 8341 8342 8343 8344 8345 8346 8347 8348 8349 8350 8351}
mol material Opaque
mol addrep top

# Molecule 260: Licorice (blue, ColorID 0)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 0
mol selection {index 9324 9325 9326 9327 9328 9329 9330 9331 9332 9333 9334 9335 9336 9337 9338 9339 9340 9341 9342 9343 9344 9345 9346 9347 9348 9349 9350 9351 9352 9353 9354 9355 9356 9357 9358 9359}
mol material Opaque
mol addrep top

# QM molecule 202: CPK representation (element colors to show C-H bonds)
mol representation CPK 1.0 0.3 12.0 12.0
mol color Name
mol selection {index 7236 7237 7238 7239 7240 7241 7242 7243 7244 7245 7246 7247 7248 7249 7250 7251 7252 7253 7254 7255 7256 7257 7258 7259 7260 7261 7262 7263 7264 7265 7266 7267 7268 7269 7270 7271}
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
set sel [atomselect top "index 4572 to 4607 5940 to 5975 6048 to 6083 7164 to 7199 7236 to 7271 7308 to 7343 8172 to 8207 8316 to 8351 9324 to 9359"]
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
# add_mol_label 128 4572 4607 1     ;# red
# add_mol_label 166 5940 5975 3     ;# orange
# add_mol_label 169 6048 6083 11    ;# purple
# add_mol_label 200 7164 7199 7     ;# green
# add_mol_label 204 7308 7343 4     ;# yellow
# add_mol_label 228 8172 8207 10    ;# cyan
# add_mol_label 232 8316 8351 9     ;# pink
# add_mol_label 260 9324 9359 14    ;# ochre
# add_mol_label 202 7236 7271 2     ;# gray (QM)

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/picene/input/viewpoint_selected.tcl"}} {
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {rotate_matrix} {[molinfo top get rotate_matrix]}"
    puts $fp "molinfo top set {center_matrix} {[molinfo top get center_matrix]}"
    puts $fp "molinfo top set {scale_matrix} {[molinfo top get scale_matrix]}"
    puts $fp "molinfo top set {global_matrix} {[molinfo top get global_matrix]}"
    close $fp
    puts "Viewpoint saved to: $filename"
}

proc load_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/picene/input/viewpoint_selected.tcl"}} {
    if {[file exists $filename]} {
        source $filename
        puts "Viewpoint loaded from: $filename"
    } else {
        puts "Viewpoint file not found: $filename"
    }
}

# Auto-load viewpoint if exists
if {[file exists "/home/takahashi/python/q4mdkit/examples/picene/input/viewpoint_selected.tcl"]} {
    load_viewpoint
}

puts ""
puts "=== Selected Molecules ==="
puts "  128: red"
puts "  166: orange"
puts "  169: purple"
puts "  200: green"
puts "  204: yellow"
puts "  228: cyan"
puts "  232: pink"
puts "  260: blue"
puts "  202: gray (QM molecule)"
puts ""
puts "Toggle representations:"
puts "  mol showrep top 0 off/on  - Molecule 128 (red)"
puts "  mol showrep top 1 off/on  - Molecule 166 (orange)"
puts "  mol showrep top 2 off/on  - Molecule 169 (purple)"
puts "  mol showrep top 3 off/on  - Molecule 200 (green)"
puts "  mol showrep top 4 off/on  - Molecule 204 (yellow)"
puts "  mol showrep top 5 off/on  - Molecule 228 (cyan)"
puts "  mol showrep top 6 off/on  - Molecule 232 (pink)"
puts "  mol showrep top 7 off/on  - Molecule 260 (blue)"
puts "  mol showrep top 8 off/on  - Molecule 202 (QM)"
puts ""
puts "Commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
puts "  graphics top delete all - Hide labels"
