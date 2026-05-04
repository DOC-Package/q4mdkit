#!/bin/bash
set -e

echo "=== Step 1: Packing system with Packmol ==="
packmol < pack_system.inp

echo "=== Step 2: Building AMBER topology with tleap ==="
tleap -f make_system.in

echo "=== Step 3: Generating .box file from tleap PDB ==="
python -c "from q4mdkit.prep.pdb2box import pdb_to_box; pdb_to_box('nile-red_acetone_leap.pdb', 'nile-red_acetone.box')"

echo "=== Step 4: Generating QM atom list ==="
python make_qmatoms.py nile-red_acetone_leap.pdb --resname QUC -o qmatoms

echo "=== Done ==="
echo "Output files:"
echo "  nile-red_acetone.pdb"
echo "  nile-red_acetone.prmtop"
echo "  nile-red_acetone.inpcrd"
echo "  nile-red_acetone_leap.pdb"
echo "  nile-red_acetone.box"
echo "  qmatoms"
