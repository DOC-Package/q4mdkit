"""
解析モジュール
シミュレーション結果の解析と可視化
"""
import logging
from typing import Dict, Any
from pathlib import Path
import numpy as np

logger = logging.getLogger(__name__)


class ResultAnalyzer:
    """結果解析クラス"""
    
    def __init__(self, results: Dict[str, list]):
        """
        Args:
            results: シミュレーション結果
        """
        self.results = results
        
    def calculate_statistics(self) -> Dict[str, float]:
        """統計量を計算"""
        stats = {}
        
        if "energy" in self.results and len(self.results["energy"]) > 0:
            energies = np.array(self.results["energy"])
            stats["energy_mean"] = float(np.mean(energies))
            stats["energy_std"] = float(np.std(energies))
            stats["energy_min"] = float(np.min(energies))
            stats["energy_max"] = float(np.max(energies))
        
        if "temperature" in self.results and len(self.results["temperature"]) > 0:
            temps = np.array(self.results["temperature"])
            stats["temperature_mean"] = float(np.mean(temps))
            stats["temperature_std"] = float(np.std(temps))
        
        return stats
    
    def save_data(self, output_path: Path):
        """データをファイルに保存"""
        logger.info(f"データを {output_path} に保存中...")
        
        # CSVファイルとして保存
        csv_path = output_path.with_suffix(".csv")
        with open(csv_path, "w") as f:
            f.write("step,energy,temperature\n")
            for i in range(len(self.results["steps"])):
                step = self.results["steps"][i]
                energy = self.results["energy"][i]
                temp = self.results["temperature"][i]
                f.write(f"{step},{energy},{temp}\n")
        
        logger.info(f"保存完了: {csv_path}")
    
    def plot_results(self, output_path: Path):
        """結果をプロット（matplotlibが必要）"""
        try:
            import matplotlib.pyplot as plt
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            steps = self.results["steps"]
            
            # エネルギープロット
            ax1.plot(steps, self.results["energy"], "b-")
            ax1.set_xlabel("Step")
            ax1.set_ylabel("Energy (kJ/mol)")
            ax1.set_title("Energy vs Time")
            ax1.grid(True)
            
            # 温度プロット
            ax2.plot(steps, self.results["temperature"], "r-")
            ax2.set_xlabel("Step")
            ax2.set_ylabel("Temperature (K)")
            ax2.set_title("Temperature vs Time")
            ax2.grid(True)
            
            plt.tight_layout()
            
            plot_path = output_path.with_suffix(".png")
            plt.savefig(plot_path, dpi=300)
            logger.info(f"プロット保存: {plot_path}")
            
        except ImportError:
            logger.warning("matplotlibがインストールされていないため、プロットをスキップします")
