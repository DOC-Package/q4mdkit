#!/bin/bash
# 全計算の一括実行スクリプト

BASEDIR=$(dirname "$0")
cd "$BASEDIR"

echo "=============================================="
echo "Pentacene Aggregate DFTB+ Calculations"
echo "=============================================="
echo "原子数: 1296 (36分子)"
echo "中心分子: 分子11 (原子361-396)"
echo ""

# GENファイルの生成（まだない場合）
if [ ! -f "pentacene.gen" ]; then
    echo ">>> GEN形式ファイルの生成..."
    python pdb2gen.py
    echo ""
fi

# 中性計算
echo "=============================================="
echo "1. Neutral Optimization"
echo "=============================================="
cd neutral
chmod +x run_opt.sh run_vib.sh
./run_opt.sh
echo ""

echo "=============================================="
echo "2. Neutral Hessian"
echo "=============================================="
./run_vib.sh
cd ..
echo ""

# カチオン計算
echo "=============================================="
echo "3. Cation Optimization (CDFT)"
echo "=============================================="
cd cation
chmod +x run_opt.sh run_vib.sh
./run_opt.sh
echo ""

echo "=============================================="
echo "4. Cation Hessian (CDFT)"
echo "=============================================="
./run_vib.sh
cd ..
echo ""

echo "=============================================="
echo "All calculations completed!"
echo "=============================================="
