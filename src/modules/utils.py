"""
ユーティリティ関数
"""
import logging
from datetime import datetime
from pathlib import Path


def setup_logging(log_level: str = "INFO", log_format: str = None):
    """ロギングをセットアップ"""
    if log_format is None:
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format=log_format,
        handlers=[
            logging.StreamHandler(),
        ]
    )


def get_timestamp() -> str:
    """タイムスタンプを取得"""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def print_header(title: str):
    """ヘッダーを表示"""
    width = 60
    print("=" * width)
    print(f"{title:^{width}}")
    print("=" * width)


def print_dict(data: dict, indent: int = 0):
    """辞書を整形して表示"""
    for key, value in data.items():
        if isinstance(value, dict):
            print(" " * indent + f"{key}:")
            print_dict(value, indent + 2)
        else:
            print(" " * indent + f"{key}: {value}")
