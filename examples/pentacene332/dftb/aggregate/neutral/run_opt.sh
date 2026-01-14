#!/bin/bash
# DFTB+による中性集合体の構造最適化

echo "=== Neutral Aggregate Optimization ==="
echo "原子数: 1296 (36分子)"
echo ""

# DFTB+の実行（最適化）
cp dftb_in_opt.hsd dftb_in.hsd
dftb+ > opt.log 2>&1

if [ $? -eq 0 ]; then
    echo "構造最適化完了"
    echo "結果: detailed.out, geom.out.gen"
else
    echo "エラー: 構造最適化に失敗しました"
    echo "opt.logを確認してください"
    exit 1
fi
