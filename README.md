# OpenMM Crystal Project

## プロジェクト概要
OpenMMを使用した結晶構造のシミュレーションプロジェクト

## ディレクトリ構造
```
.
├── src/              # メインのソースコード
│   ├── main.py      # エントリーポイント
│   ├── config.py    # 設定管理
│   └── modules/     # 機能別モジュール
├── tests/           # テストコード
├── data/            # データファイル
│   ├── input/       # 入力データ
│   └── output/      # 出力結果
├── docs/            # ドキュメント
└── requirements.txt # 依存パッケージ
```

## 使い方
```bash
# 依存パッケージのインストール
pip install -r requirements.txt

# プログラムの実行
python src/main.py
```

## 開発
```bash
# テストの実行
python -m pytest tests/
```
