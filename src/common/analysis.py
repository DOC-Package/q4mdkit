"""
Analysis module
Analyze and visualize simulation results
"""

import logging
from typing import Dict, Any
from pathlib import Path
import numpy as np

logger = logging.getLogger(__name__)


class ResultAnalyzer:
    """Class for analyzing simulation results.

    Provides simple statistics, CSV export, and plotting utilities.
    """

    def __init__(self, results: Dict[str, list]):
        """Initialize with a results dictionary.

        Args:
            results: Dictionary containing simulation data (keys: 'steps', 'energy', 'temperature', ...)
        """
        self.results = results

    def calculate_statistics(self) -> Dict[str, float]:
        """Compute basic statistics for available quantities.

        Returns a dict with mean/std/min/max where applicable.
        """
        stats: Dict[str, float] = {}

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
        """Save results to a CSV file.

        The CSV will have columns: step,energy,temperature
        """
        logger.info(f"Saving data to {output_path}...")

        csv_path = output_path.with_suffix(".csv")
        with open(csv_path, "w") as f:
            f.write("step,energy,temperature\n")
            steps = self.results.get("steps", [])
            energies = self.results.get("energy", [])
            temps = self.results.get("temperature", [])
            for i in range(len(steps)):
                step = steps[i]
                energy = energies[i] if i < len(energies) else ""
                temp = temps[i] if i < len(temps) else ""
                f.write(f"{step},{energy},{temp}\n")

        logger.info(f"Saved: {csv_path}")

    def plot_results(self, output_path: Path):
        """Plot energy and temperature time series and save as PNG.

        Requires matplotlib.
        """
        try:
            import matplotlib.pyplot as plt

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

            steps = self.results.get("steps", [])

            # Energy plot
            ax1.plot(steps, self.results.get("energy", []), "b-")
            ax1.set_xlabel("Step")
            ax1.set_ylabel("Energy (kJ/mol)")
            ax1.set_title("Energy vs Time")
            ax1.grid(True)

            # Temperature plot
            ax2.plot(steps, self.results.get("temperature", []), "r-")
            ax2.set_xlabel("Step")
            ax2.set_ylabel("Temperature (K)")
            ax2.set_title("Temperature vs Time")
            ax2.grid(True)

            plt.tight_layout()

            plot_path = output_path.with_suffix(".png")
            plt.savefig(plot_path, dpi=300)
            logger.info(f"Plot saved: {plot_path}")

        except ImportError:
            logger.warning("matplotlib is not installed; skipping plotting")