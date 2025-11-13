"""
Utility functions
"""
import logging
from datetime import datetime


def setup_logging(level: str, format_string: str):
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=format_string,
    )
    # Adjust third-party log levels
    logging.getLogger("matplotlib").setLevel(logging.WARNING)


def get_timestamp() -> str:
    """Return a timestamp string (YYYYMMDD_HHMMSS)."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def print_header(text: str):
    """Print a centered header banner."""
    width = 60
    print("\n" + "=" * width)
    print(text.center(width))
    print("=" * width + "\n")


def print_dict(d: dict, indent: int = 0):
    """Print a dictionary with indentation."""
    for key, value in d.items():
        print(f"{'  ' * indent}{key}: {value}")
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
