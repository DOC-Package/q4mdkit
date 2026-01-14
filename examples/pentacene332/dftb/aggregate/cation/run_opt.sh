#!/bin/bash
# DFTB+によるカチオン集合体の構造最適化（中心分子に電荷拘束）

echo "=== Cation Aggregate Optimization (CDFT) ==="
echo "原子数: 1296 (36分子)"
echo "中心分子: 分子11 (原子361-396)"
echo "電荷拘束: +1 on 中心分子"
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
