"""
Base classes for simulation engines
Defines a common interface for OpenMM and pyCHARMM engines
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import numpy as np
from pathlib import Path


class BaseSimulator(ABC):
    """Abstract base class for simulators.

    All simulation engine implementations (OpenMM, pyCHARMM, etc.) should
    inherit from this class and implement the abstract methods.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Args:
            config: シミュレーション設定辞書
        """
        self.config = config
        self.temperature = config.get("temperature", 300)
        self.pressure = config.get("pressure", 1.0)
        self.timestep = config.get("timestep", 2.0)
        self.total_steps = config.get("total_steps", 10000)
        self.output_interval = config.get("output_interval", 100)
        
        self.current_step = 0
        self.results = {
            "steps": [],
            "energy": [],
            "temperature": [],
            "pressure": [],
        }
        self._is_setup = False
        
    @abstractmethod
    def setup_system(self, positions: np.ndarray, box_vectors: np.ndarray):
        """Set up the system (engine-specific).

        Args:
            positions: Atomic positions (N, 3)
            box_vectors: Box vectors (3, 3)
        """
        pass
    
    @abstractmethod
    def setup_integrator(self):
        """Set up the integrator (engine-specific)."""
        pass
    
    @abstractmethod
    def minimize_energy(self, max_iterations: int = 1000, tolerance: float = 10.0):
        """Perform energy minimization (engine-specific).

        Args:
            max_iterations: Maximum number of iterations
            tolerance: Convergence tolerance (kJ/mol/nm)
        """
        pass
    
    @abstractmethod
    def run_dynamics(self):
        """Run molecular dynamics (engine-specific)."""
        pass
    
    @abstractmethod
    def get_current_state(self) -> Dict[str, Any]:
        """Get the current state (engine-specific).

        Returns:
            State information (energy, temperature, positions, etc.)
        """
        pass
    
    @abstractmethod
    def save_trajectory(self, output_path: Path):
        """Save trajectory (engine-specific).

        Args:
            output_path: Output file path
        """
        pass
    
    # Common methods (shared across engines)
    
    def setup(self, positions: np.ndarray, box_vectors: np.ndarray):
        """High-level setup for the simulation.

        Args:
            positions: Atomic positions
            box_vectors: Box vectors
        """
        self.setup_system(positions, box_vectors)
        self.setup_integrator()
        self._is_setup = True
        
    def run(self):
        """Execute the simulation (common flow)."""
        if not self._is_setup:
            raise RuntimeError("Call setup() before running the simulation")
        
        # エネルギー最小化
        self.minimize_energy()
        
        # 分子動力学実行
        self.run_dynamics()
        
    def record_data(self, step: int, state: Dict[str, Any]):
        """Record data (common handling).

        Args:
            step: Step number
            state: State information
        """
        self.results["steps"].append(step)
        self.results["energy"].append(state.get("energy", 0.0))
        self.results["temperature"].append(state.get("temperature", 0.0))
        self.results["pressure"].append(state.get("pressure", 0.0))
        
    def get_results(self) -> Dict[str, list]:
        """Return simulation results.

        Returns:
            A dictionary with recorded simulation data.
        """
        return self.results
    
    def get_engine_name(self) -> str:
        """Return the engine name (derived from the class name)."""
        return self.__class__.__name__
