"""
設定管理モジュール
プロジェクト全体で使用する設定を一元管理
"""
import os
from pathlib import Path

# プロジェクトのルートディレクトリ
PROJECT_ROOT = Path(__file__).parent.parent

# データディレクトリ
DATA_DIR = PROJECT_ROOT / "data"
INPUT_DIR = DATA_DIR / "input"
OUTPUT_DIR = DATA_DIR / "output"

# シミュレーションパラメータ
SIMULATION_CONFIG = {
    "temperature": 300,  # K
    "pressure": 1.0,     # bar
    "timestep": 2.0,     # fs
    "total_steps": 10000,
    "output_interval": 100,
}

# ログ設定
LOG_LEVEL = "INFO"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def ensure_directories():
    """必要なディレクトリが存在することを確認"""
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_output_path(filename: str) -> Path:
    """出力ファイルのパスを取得"""
    return OUTPUT_DIR / filename


def get_input_path(filename: str) -> Path:
    """入力ファイルのパスを取得"""
    return INPUT_DIR / filename
