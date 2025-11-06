"""
メインエントリーポイント
プログラムの実行開始点
"""
import logging
import sys
from pathlib import Path

# 親ディレクトリをパスに追加（モジュールインポートのため）
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    SIMULATION_CONFIG,
    LOG_LEVEL,
    LOG_FORMAT,
    ensure_directories,
    get_output_path,
)
from modules.utils import setup_logging, get_timestamp, print_header, print_dict
from modules.crystal_structure import CrystalStructure
from modules.simulation import Simulator
from modules.analysis import ResultAnalyzer

logger = logging.getLogger(__name__)


def main():
    """メイン実行関数"""
    # ロギング設定
    setup_logging(LOG_LEVEL, LOG_FORMAT)
    
    # ディレクトリ確認
    ensure_directories()
    
    print_header("OpenMM Crystal Simulation")
    
    # 1. 結晶構造の生成
    logger.info("結晶構造を生成中...")
    crystal = CrystalStructure(
        lattice_constant=0.4,  # nm
        num_cells=(3, 3, 3)
    )
    positions = crystal.generate_fcc_lattice()
    logger.info(f"原子数: {crystal.get_num_atoms()}")
    logger.info(f"ボックスサイズ: {crystal.get_box_vectors()}")
    
    # 2. シミュレーションのセットアップと実行
    logger.info("\n" + "="*60)
    simulator = Simulator(SIMULATION_CONFIG)
    simulator.setup()
    simulator.run()
    
    # 3. 結果の解析
    logger.info("\n" + "="*60)
    logger.info("結果を解析中...")
    results = simulator.get_results()
    analyzer = ResultAnalyzer(results)
    
    # 統計量を計算
    stats = analyzer.calculate_statistics()
    print("\n統計量:")
    print_dict(stats, indent=2)
    
    # 4. 結果の保存
    logger.info("\n" + "="*60)
    timestamp = get_timestamp()
    output_base = get_output_path(f"simulation_{timestamp}")
    
    analyzer.save_data(output_base)
    analyzer.plot_results(output_base)
    
    print_header("完了")
    logger.info(f"結果は {output_base.parent} に保存されました")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\n中断されました")
        sys.exit(0)
    except Exception as e:
        logger.error(f"エラーが発生しました: {e}", exc_info=True)
        sys.exit(1)
