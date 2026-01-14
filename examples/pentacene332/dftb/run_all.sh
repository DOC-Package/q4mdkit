#!/bin/bash
# Run geometry optimization and normal mode calculations for both neutral and cation

echo "======================================"
echo "Starting calculations for pentacene"
echo "======================================"

# Neutral calculations
echo ""
echo "--- Neutral pentacene ---"
cd neutral
echo "Running geometry optimization..."
bash run_opt.sh
echo "Running vibrational analysis..."
bash run_vib.sh
cd ..

# Cation calculations
echo ""
echo "--- Pentacene cation ---"
cd cation
echo "Running geometry optimization..."
bash run_opt.sh
echo "Running vibrational analysis..."
bash run_vib.sh
cd ..

echo ""
echo "======================================"
echo "All calculations completed!"
echo "======================================"
echo "Results are in:"
echo "  - neutral/: Neutral pentacene results"
echo "  - cation/: Cation pentacene results"
echo ""
echo "Key output files:"
echo "  - geom.out.xyz: Optimized geometry"
echo "  - hessian.out: Hessian matrix"
echo "  - vibrations.tag: Vibrational frequencies and modes"
