#!/bin/bash
set -e

echo "=== Step 1: Packing system with Packmol ==="
packmol < pack_system.inp

echo "=== Step 2: Building AMBER topology with tleap ==="
tleap -f make_system.in

echo "=== Step 3: Generating .box file from tleap PDB ==="
python -c "from q4mdkit.prep.pdb2box import pdb_to_box; pdb_to_box('pdi_acetonitrile_leap.pdb', 'pdi_acetonitrile.box')"

echo "=== Step 4: Generating QM atom list ==="
python make_qmatoms.py pdi_acetonitrile_leap.pdb --resname PDI -o qmatoms

echo "=== Done ==="
