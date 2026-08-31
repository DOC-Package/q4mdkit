#!/bin/bash
set -e

echo "=== Step 1: Packing system with Packmol ==="
packmol < pack_system.inp

echo "=== Step 2: Building AMBER topology with tleap ==="
tleap -f make_system.in

echo "=== Step 3: Generating .box file from tleap PDB ==="
python3 -c "import sys; sys.path.insert(0, '../../../q4mdkit/prep'); from pdb2box import pdb_to_box; pdb_to_box('vg30_acetonitrile_leap.pdb', 'vg30_acetonitrile.box')"

echo "=== Step 4: Generating QM atom list ==="
python3 make_qmatoms.py vg30_acetonitrile_leap.pdb --resname VG3 -o qmatoms

echo "=== Done ==="
echo "Output files:"
echo "  vg30_acetonitrile.pdb"
echo "  vg30_acetonitrile.prmtop"
echo "  vg30_acetonitrile.inpcrd"
echo "  vg30_acetonitrile_leap.pdb"
echo "  vg30_acetonitrile.box"
echo "  qmatoms"
