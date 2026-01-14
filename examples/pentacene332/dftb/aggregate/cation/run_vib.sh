#!/bin/bash
# DFTB+によるカチオン集合体のヘッシアン（振動数）計算

echo "=== Cation Aggregate Hessian Calculation (CDFT) ==="
echo "原子数: 1296 (36分子)"
echo "中心分子: 分子11 (原子361-396)"
echo "電荷拘束: +1 on 中心分子"
echo ""

# 最適化構造が存在するか確認
if [ ! -f "geom.out.gen" ]; then
    echo "エラー: geom.out.gen が見つかりません"
    echo "先に run_opt.sh を実行してください"
    exit 1
fi

# DFTB+の実行（ヘッシアン）
cp dftb_in_vib.hsd dftb_in.hsd
dftb+ > vib.log 2>&1

if [ $? -eq 0 ]; then
    echo "ヘッシアン計算完了"
    echo "結果: hessian.out, vibrations.tag"
else
    echo "エラー: ヘッシアン計算に失敗しました"
    echo "vib.logを確認してください"
    exit 1
fi
