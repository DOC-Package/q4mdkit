#!/usr/bin/env python3
"""
Generate PDB file for visualizing QM/MM system.
Highlight QM region using B-factor/occupancy values.
"""

import numpy as np

def read_intlist_from_file(filename):
    """Read a list of integers from a file."""
    with open(filename, 'r') as f:
        return [int(x) for x in f.read().strip().split()]

def write_pdb_with_qm_highlight(grofile, qmatoms_file, output_pdb):
    """
    Read GRO file and convert to PDB.
    Highlight QM region using B-factor values.
    """
    # Read list of QM atoms
    qmatoms = set(read_intlist_from_file(qmatoms_file))
    
    # Read GRO file
    with open(grofile, 'r') as f:
        lines = f.readlines()
    
    # Title and number of atoms
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    # Write to PDB file
    with open(output_pdb, 'w') as f:
        f.write(f"TITLE     {title} - QM/MM visualization\n")
        f.write(f"REMARK    QM atoms: {len(qmatoms)} atoms\n")
        f.write(f"REMARK    MM atoms: {natoms - len(qmatoms)} atoms\n")
        f.write(f"REMARK    B-factor: 100.0 = QM region, 0.0 = MM region\n")
        f.write(f"REMARK    Occupancy: 1.0 = QM region, 0.5 = MM region\n")
        
        # Read and write atomic coordinates
        for i in range(2, 2 + natoms):
            line = lines[i]
            
            # Parse GRO file (fixed-width format)
            resnum = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomnum = int(line[15:20].strip())
            x = float(line[20:28].strip()) * 10.0  # nm -> Angstrom
            y = float(line[28:36].strip()) * 10.0
            z = float(line[36:44].strip()) * 10.0
            
            # Change values based on whether atom is in QM region
            atom_index = atomnum - 1  # 0-indexed
            if atom_index in qmatoms:
                bfactor = 100.0
                occupancy = 1.0
                element = atomname[0].upper()
            else:
                bfactor = 0.0
                occupancy = 0.5
                element = atomname[0].upper()
            
            # Write in PDB format
            f.write(f"ATOM  {atomnum:5d} {atomname:4s} {resname:3s} A{resnum:4d}    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}\n")
        
        # Box information (triclinic support)
        box_line = lines[2 + natoms].strip().split()
        
        if len(box_line) >= 9:
            # Triclinic: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            # Box vectors: v1 = (v1x, 0, 0)
            #              v2 = (v2x, v2y, 0)
            #              v3 = (v3x, v3y, v3z)
            v1x = float(box_line[0]) * 10.0  # nm -> Angstrom
            v2y = float(box_line[1]) * 10.0
            v3z = float(box_line[2]) * 10.0
            v1y = float(box_line[3]) * 10.0
            v1z = float(box_line[4]) * 10.0
            v2x = float(box_line[5]) * 10.0
            v2z = float(box_line[6]) * 10.0
            v3x = float(box_line[7]) * 10.0
            v3y = float(box_line[8]) * 10.0
            
            # Calculate lattice constants a, b, c and angles alpha, beta, gamma
            import numpy as np
            a = v1x
            b = np.sqrt(v2x**2 + v2y**2)
            c = np.sqrt(v3x**2 + v3y**2 + v3z**2)
            
            # Angle calculation
            alpha = np.degrees(np.arccos((v2x*v3x + v2y*v3y) / (b * c)))
            beta = np.degrees(np.arccos(v3x / c))
            gamma = np.degrees(np.arccos(v2x / b))
            
        else:
            # Orthorhombic
            a = float(box_line[0]) * 10.0
            b = float(box_line[1]) * 10.0
            c = float(box_line[2]) * 10.0
            alpha = 90.0
            beta = 90.0
            gamma = 90.0
        
        f.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n")
        f.write("END\n")
    
    print(f"✓ PDB file created: {output_pdb}")
    print(f"  QM atoms: {len(qmatoms)}")
    print(f"  Total atoms: {natoms}")
    print(f"  QM region: B-factor=100.0, MM region: B-factor=0.0")

def write_vmd_script(pdb_file, script_file):
    """Generate VMD visualization script."""
    with open(script_file, 'w') as f:
        f.write(f"""# VMD visualization script for QM/MM system
# Usage: vmd -e {script_file}

# Load structure
mol new {pdb_file} type pdb waitfor all

# Color by B-factor (QM/MM region)
mol delrep 0 top
mol representation Lines 1.000000
mol color Beta
mol selection {{beta < 50}}
mol material Opaque
mol addrep top

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

puts "QM region: CPK representation (colored by element)"
puts "MM region: Lines (colored by B-factor)"
puts ""
puts "Useful commands:"
puts "  mol showrep top 0 off    - Hide MM region"
puts "  mol showrep top 1 off    - Hide QM region"
""")
    
    print(f"✓ VMD script created: {script_file}")
    print(f"  Usage: vmd -e {script_file}")

def main():
    # File names
    grofile = "pentacene.gro"
    qmatoms_file = "qmatoms"
    output_pdb = "pentacene_qmmm.pdb"
    vmd_script = "visualize_vmd.tcl"
    
    print("=" * 60)
    print("QM/MM Visualization File Generation")
    print("=" * 60)
    
    # Generate PDB file (highlight QM region)
    write_pdb_with_qm_highlight(grofile, qmatoms_file, output_pdb)
    print()
    
    # Generate VMD script
    write_vmd_script(output_pdb, vmd_script)
    print()
    
    print("=" * 60)
    print("Visualization methods:")
    print("=" * 60)
    print(f"1. VMD:     vmd -e {vmd_script}")
    print(f"2. Browser: Run 'python generate_visualization.py' and")
    print(f"            open visualize_browser.html")
    print(f"3. Manual:  Open {output_pdb} with any visualization software")
    print("            B-factor > 50 indicates QM region")
    print()
    print("=" * 60)

if __name__ == "__main__":
    main()
