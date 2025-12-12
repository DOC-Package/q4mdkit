#!/usr/bin/env python3
"""
Generate PDB file for visualizing QM/MM system.
Highlight QM region using B-factor/occupancy values.
"""

import numpy as np
from pathlib import Path

def read_intlist_from_file(filename):
    """Read a list of integers from a file."""
    with open(filename, 'r') as f:
        return [int(x) for x in f.read().strip().split()]


def read_gro_atoms(grofile):
    """
    Read GRO file and return atom information.
    
    Returns:
        atoms: list of dicts with resnum, resname, atomname, atomnum, x, y, z (in Angstrom)
        box_info: dict with a, b, c, alpha, beta, gamma
        title: str
    """
    with open(grofile, 'r') as f:
        lines = f.readlines()
    
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    atoms = []
    for i in range(2, 2 + natoms):
        line = lines[i]
        atoms.append({
            'resnum': int(line[0:5].strip()),
            'resname': line[5:10].strip(),
            'atomname': line[10:15].strip(),
            'atomnum': int(line[15:20].strip()),
            'x': float(line[20:28].strip()) * 10.0,  # nm -> Angstrom
            'y': float(line[28:36].strip()) * 10.0,
            'z': float(line[36:44].strip()) * 10.0,
        })
    
    # Box information
    box_line = lines[2 + natoms].strip().split()
    if len(box_line) >= 9:
        v1x = float(box_line[0]) * 10.0
        v2y = float(box_line[1]) * 10.0
        v3z = float(box_line[2]) * 10.0
        v2x = float(box_line[5]) * 10.0
        v3x = float(box_line[7]) * 10.0
        v3y = float(box_line[8]) * 10.0
        
        a = v1x
        b = np.sqrt(v2x**2 + v2y**2)
        c = np.sqrt(v3x**2 + v3y**2 + v3z**2)
        alpha = np.degrees(np.arccos((v2x*v3x + v2y*v3y) / (b * c))) if b*c > 0 else 90.0
        beta = np.degrees(np.arccos(v3x / c)) if c > 0 else 90.0
        gamma = np.degrees(np.arccos(v2x / b)) if b > 0 else 90.0
    else:
        a = float(box_line[0]) * 10.0
        b = float(box_line[1]) * 10.0
        c = float(box_line[2]) * 10.0
        alpha = beta = gamma = 90.0
    
    box_info = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}
    
    return atoms, box_info, title


def read_pdb_atoms(pdbfile):
    """
    Read PDB file and return atom information.
    
    Returns:
        atoms: list of dicts with resnum, resname, atomname, atomnum, x, y, z (in Angstrom)
        box_info: dict with a, b, c, alpha, beta, gamma
        title: str
    """
    atoms = []
    title = "PDB file"
    box_info = {'a': 0, 'b': 0, 'c': 0, 'alpha': 90, 'beta': 90, 'gamma': 90}
    
    with open(pdbfile, 'r') as f:
        for line in f:
            if line.startswith("TITLE"):
                title = line[10:].strip()
            elif line.startswith("CRYST1"):
                box_info['a'] = float(line[6:15].strip())
                box_info['b'] = float(line[15:24].strip())
                box_info['c'] = float(line[24:33].strip())
                box_info['alpha'] = float(line[33:40].strip())
                box_info['beta'] = float(line[40:47].strip())
                box_info['gamma'] = float(line[47:54].strip())
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append({
                    'atomnum': int(line[6:11].strip()),
                    'atomname': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'resnum': int(line[22:26].strip()),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                })
    
    return atoms, box_info, title


def read_structure_atoms(filepath):
    """
    Read GRO or PDB file based on extension.
    
    Returns:
        atoms: list of dicts
        box_info: dict
        title: str
    """
    ext = Path(filepath).suffix.lower()
    if ext == '.gro':
        return read_gro_atoms(filepath)
    elif ext == '.pdb':
        return read_pdb_atoms(filepath)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .gro or .pdb")

def write_pdb_with_qm_highlight(structure_file, qmatoms_file, output_pdb):
    """
    Read GRO or PDB file and convert to PDB.
    Highlight QM region using B-factor values.
    """
    # Read list of QM atoms
    qmatoms = set(read_intlist_from_file(qmatoms_file))
    
    # Read structure file
    atoms, box_info, title = read_structure_atoms(structure_file)
    natoms = len(atoms)
    
    # Write to PDB file
    with open(output_pdb, 'w') as f:
        f.write(f"TITLE     {title} - QM/MM visualization\n")
        f.write(f"REMARK    QM atoms: {len(qmatoms)} atoms\n")
        f.write(f"REMARK    MM atoms: {natoms - len(qmatoms)} atoms\n")
        f.write(f"REMARK    B-factor: 100.0 = QM region, 0.0 = MM region\n")
        f.write(f"REMARK    Occupancy: 1.0 = QM region, 0.5 = MM region\n")
        
        # Write atomic coordinates
        for atom in atoms:
            atom_index = atom['atomnum'] - 1  # 0-indexed
            if atom_index in qmatoms:
                bfactor = 100.0
                occupancy = 1.0
            else:
                bfactor = 0.0
                occupancy = 0.5
            
            element = atom['atomname'][0].upper()
            
            f.write(f"ATOM  {atom['atomnum']:5d} {atom['atomname']:4s} {atom['resname']:3s} A{atom['resnum']:4d}    "
                   f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}\n")
        
        # Box information
        f.write(f"CRYST1{box_info['a']:9.3f}{box_info['b']:9.3f}{box_info['c']:9.3f}"
               f"{box_info['alpha']:7.2f}{box_info['beta']:7.2f}{box_info['gamma']:7.2f} P 1           1\n")
        f.write("END\n")
    
    print(f"✓ PDB file created: {output_pdb}")
    print(f"  QM atoms: {len(qmatoms)}")
    print(f"  Total atoms: {natoms}")
    print(f"  QM region: B-factor=100.0, MM region: B-factor=0.0")

def write_vmd_script(pdb_file, script_file):
    """Generate VMD visualization script (standard QM/MM view)."""
    script_dir = Path(script_file).parent
    viewpoint_file = script_dir / "viewpoint.tcl"
    
    with open(script_file, 'w') as f:
        f.write(f"""# VMD visualization script for QM/MM system
# Usage: vmd -e {script_file}

# Load structure
mol new {pdb_file} type pdb waitfor all

# Color by B-factor (QM/MM region)
mol delrep 0 top

# MM region: thin lines
mol representation Lines 1.000000
mol color Beta
mol selection {{beta < 50}}
mol material Opaque
mol addrep top

# QM region: CPK balls
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol color Name
mol selection {{beta > 50}}
mol material Opaque
mol addrep top

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on QM region
set qm_sel [atomselect top "beta > 50"]
set center [measure center $qm_sel]
set box_size 30
molinfo top set center [list $center]

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {{rotate_matrix}} {{[molinfo top get rotate_matrix]}}"
    puts $fp "molinfo top set {{center_matrix}} {{[molinfo top get center_matrix]}}"
    puts $fp "molinfo top set {{scale_matrix}} {{[molinfo top get scale_matrix]}}"
    puts $fp "molinfo top set {{global_matrix}} {{[molinfo top get global_matrix]}}"
    close $fp
    puts "Viewpoint saved to: $filename"
}}

proc load_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    if {{[file exists $filename]}} {{
        source $filename
        puts "Viewpoint loaded from: $filename"
    }} else {{
        puts "Viewpoint file not found: $filename"
    }}
}}

# ============================================================
# Image save functions
# ============================================================
proc save_image {{{{filename "snapshot.tga"}}}} {{
    render snapshot $filename
    puts "Snapshot saved to: $filename"
}}

proc save_eps {{{{filename "snapshot.eps"}}}} {{
    render PostScript $filename
    puts "EPS (vector) saved to: $filename"
}}

proc save_tachyon {{{{filename "snapshot.tga"}}}} {{
    render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga"
    puts "High-quality image saved to: $filename.tga"
}}

# Auto-load viewpoint if exists
if {{[file exists "{viewpoint_file}"]}} {{
    load_viewpoint
    puts "Auto-loaded viewpoint from {viewpoint_file}"
}}

puts "QM region: CPK representation (colored by element)"
puts "MM region: Lines (colored by B-factor)"
puts ""
puts "Useful commands:"
puts "  mol showrep top 0 off    - Hide MM region"
puts "  mol showrep top 1 off    - Hide QM region"
puts ""
puts "Viewpoint commands:"
puts "  save_viewpoint           - Save current view to {viewpoint_file}"
puts "  load_viewpoint           - Load view from {viewpoint_file}"
""")
    
    print(f"✓ VMD script created: {script_file}")
    print(f"  Usage: vmd -e {script_file}")


def write_vmd_script_allqm(pdb_file, script_file):
    """Generate VMD script showing all atoms in QM-like style (CPK)."""
    script_dir = Path(script_file).parent
    viewpoint_file = script_dir / "viewpoint.tcl"
    
    with open(script_file, 'w') as f:
        f.write(f"""# VMD visualization script - All QM-like view
# Shows entire system with CPK representation (QM-style)
# Usage: vmd -e {script_file}

# Load structure
mol new {pdb_file} type pdb waitfor all

# Delete default representation
mol delrep 0 top

# All atoms: CPK balls (QM-like style)
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol color Name
mol selection {{all}}
mol material Opaque
mol addrep top

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on QM region
set qm_sel [atomselect top "beta > 50"]
if {{[$qm_sel num] > 0}} {{
    set center [measure center $qm_sel]
    molinfo top set center [list $center]
}}

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {{rotate_matrix}} {{[molinfo top get rotate_matrix]}}"
    puts $fp "molinfo top set {{center_matrix}} {{[molinfo top get center_matrix]}}"
    puts $fp "molinfo top set {{scale_matrix}} {{[molinfo top get scale_matrix]}}"
    puts $fp "molinfo top set {{global_matrix}} {{[molinfo top get global_matrix]}}"
    close $fp
    puts "Viewpoint saved to: $filename"
}}

proc load_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    if {{[file exists $filename]}} {{
        source $filename
        puts "Viewpoint loaded from: $filename"
    }} else {{
        puts "Viewpoint file not found: $filename"
    }}
}}

# ============================================================
# Image save functions
# ============================================================
proc save_image {{{{filename "snapshot.tga"}}}} {{
    render snapshot $filename
    puts "Snapshot saved to: $filename"
}}

proc save_eps {{{{filename "snapshot.eps"}}}} {{
    render PostScript $filename
    puts "EPS (vector) saved to: $filename"
}}

proc save_tachyon {{{{filename "snapshot.tga"}}}} {{
    render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga"
    puts "High-quality image saved to: $filename.tga"
}}

# Auto-load viewpoint if exists
if {{[file exists "{viewpoint_file}"]}} {{
    load_viewpoint
    puts "Auto-loaded viewpoint from {viewpoint_file}"
}}

puts "All atoms shown with CPK representation (QM-like style)"
puts "Colored by element name"
puts ""
puts "Viewpoint commands:"
puts "  save_viewpoint           - Save current view to {viewpoint_file}"
puts "  load_viewpoint           - Load view from {viewpoint_file}"
""")
    
    print(f"✓ VMD script (All-QM style) created: {script_file}")


def write_vmd_script_pointcharge(pdb_file, script_file):
    """Generate VMD script showing QM region with CPK and MM region as point charges."""
    script_dir = Path(script_file).parent
    viewpoint_file = script_dir / "viewpoint.tcl"
    
    with open(script_file, 'w') as f:
        f.write(f"""# VMD visualization script - QM + MM point charges
# QM region: CPK representation
# MM region: Small spheres (point charge representation)
# Usage: vmd -e {script_file}

# Load structure
mol new {pdb_file} type pdb waitfor all

# Delete default representation
mol delrep 0 top

# MM region: Small points (point charge representation)
mol representation Points 3.000000
mol color ColorID 8  ; # white/gray color for point charges
mol selection {{beta < 50}}
mol material Opaque
mol addrep top

# Alternative: Very small VDW spheres for point charges
mol representation VDW 0.100000 12.000000
mol color ColorID 2  ; # gray
mol selection {{beta < 50}}
mol material Transparent
mol addrep top

# QM region: CPK balls (normal atomic view)
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol color Name
mol selection {{beta > 50}}
mol material Opaque
mol addrep top

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on QM region
set qm_sel [atomselect top "beta > 50"]
if {{[$qm_sel num] > 0}} {{
    set center [measure center $qm_sel]
    molinfo top set center [list $center]
}}

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {{rotate_matrix}} {{[molinfo top get rotate_matrix]}}"
    puts $fp "molinfo top set {{center_matrix}} {{[molinfo top get center_matrix]}}"
    puts $fp "molinfo top set {{scale_matrix}} {{[molinfo top get scale_matrix]}}"
    puts $fp "molinfo top set {{global_matrix}} {{[molinfo top get global_matrix]}}"
    close $fp
    puts "Viewpoint saved to: $filename"
}}

proc load_viewpoint {{{{filename "{viewpoint_file}"}}}} {{
    if {{[file exists $filename]}} {{
        source $filename
        puts "Viewpoint loaded from: $filename"
    }} else {{
        puts "Viewpoint file not found: $filename"
    }}
}}

# ============================================================
# Image save functions
# ============================================================
proc save_image {{{{filename "snapshot.tga"}}}} {{
    render snapshot $filename
    puts "Snapshot saved to: $filename"
}}

proc save_eps {{{{filename "snapshot.eps"}}}} {{
    render PostScript $filename
    puts "EPS (vector) saved to: $filename"
}}

proc save_tachyon {{{{filename "snapshot.tga"}}}} {{
    render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga"
    puts "High-quality image saved to: $filename.tga"
}}

# Auto-load viewpoint if exists
if {{[file exists "{viewpoint_file}"]}} {{
    load_viewpoint
    puts "Auto-loaded viewpoint from {viewpoint_file}"
}}

puts "QM region: CPK representation (colored by element)"
puts "MM region: Point charges (small gray spheres)"
puts ""
puts "Useful commands:"
puts "  mol showrep top 0 off    - Hide MM points"
puts "  mol showrep top 1 off    - Hide MM small spheres"
puts "  mol showrep top 2 off    - Hide QM region"
puts ""
puts "To change point charge size:"
puts "  mol modstyle 0 top Points 5.0     - Larger points"
puts "  mol modstyle 1 top VDW 0.05 12    - Smaller spheres"
puts ""
puts "Viewpoint commands:"
puts "  save_viewpoint           - Save current view to {viewpoint_file}"
puts "  load_viewpoint           - Load view from {viewpoint_file}"
""")
    
    print(f"✓ VMD script (Point charge style) created: {script_file}")


def make_vmd_qmmm(structure_file, qmatoms_file, output_pdb=None, vmd_script=None):
    """
    Generate PDB and VMD scripts for QM/MM visualization.
    
    Creates three VMD scripts:
    1. Standard QM/MM view (QM=CPK, MM=Lines)
    2. All-QM style (entire system as CPK)
    3. Point charge style (QM=CPK, MM=point charges)
    
    Args:
        structure_file: Input GRO or PDB file
        qmatoms_file: File containing QM atom indices (0-indexed, space-separated)
        output_pdb: Output PDB file (default: {basename}_qmmm.pdb)
        vmd_script: Output VMD script basename (default: visualize_vmd.tcl)
    
    Returns:
        dict with 'pdb' and 'vmd_scripts' paths
    """
    structure_path = Path(structure_file)
    
    if output_pdb is None:
        output_pdb = str(structure_path.parent / f"{structure_path.stem}_qmmm.pdb")
    if vmd_script is None:
        vmd_script = str(structure_path.parent / "visualize_vmd.tcl")
    
    # Generate script filenames
    vmd_script_path = Path(vmd_script)
    vmd_script_allqm = str(vmd_script_path.parent / f"{vmd_script_path.stem}_allqm.tcl")
    vmd_script_pc = str(vmd_script_path.parent / f"{vmd_script_path.stem}_pointcharge.tcl")
    
    print("=" * 60)
    print("QM/MM Visualization File Generation")
    print("=" * 60)
    
    # Generate PDB file (highlight QM region)
    write_pdb_with_qm_highlight(structure_file, qmatoms_file, output_pdb)
    print()
    
    # Generate VMD scripts (3 versions)
    write_vmd_script(output_pdb, vmd_script)
    write_vmd_script_allqm(output_pdb, vmd_script_allqm)
    write_vmd_script_pointcharge(output_pdb, vmd_script_pc)
    print()
    
    print("=" * 60)
    print("Visualization methods:")
    print("=" * 60)
    print(f"1. Standard QM/MM:    vmd -e {vmd_script}")
    print(f"   (QM=CPK balls, MM=Lines)")
    print()
    print(f"2. All-QM style:      vmd -e {vmd_script_allqm}")
    print(f"   (All atoms as CPK balls)")
    print()
    print(f"3. Point charge:      vmd -e {vmd_script_pc}")
    print(f"   (QM=CPK balls, MM=point charges)")
    print()
    print(f"4. Manual: Open {output_pdb} with any viewer")
    print("   B-factor > 50 indicates QM region")
    print()
    print("=" * 60)
    
    return {
        'pdb': output_pdb,
        'vmd_script': vmd_script,
        'vmd_script_allqm': vmd_script_allqm,
        'vmd_script_pointcharge': vmd_script_pc
    }


def main():
    # File names
    grofile = "pentacene.gro"
    qmatoms_file = "qmatoms"
    output_pdb = "pentacene_qmmm.pdb"
    vmd_script = "visualize_vmd.tcl"
    
    make_vmd_qmmm(grofile, qmatoms_file, output_pdb, vmd_script)

if __name__ == "__main__":
    main()
