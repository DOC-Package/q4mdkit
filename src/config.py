"""
Configuration module
Centralized settings used across the project
"""
import os
from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent

# Data directories
DATA_DIR = PROJECT_ROOT / "data"
INPUT_DIR = DATA_DIR / "input"
OUTPUT_DIR = DATA_DIR / "output"

# Simulation parameters
SIMULATION_CONFIG = {
    "temperature": 300,  # K
    "pressure": 1.0,     # bar
    "timestep": 2.0,     # fs
    "total_steps": 10000,
    "output_interval": 100,
}

# Simulation engine selection
# Set to 'openmm' or 'pycharmm'
SIMULATION_ENGINE = "openmm"

# Logging configuration
LOG_LEVEL = "INFO"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def ensure_directories():
    """Ensure necessary directories exist."""
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_output_path(filename: str) -> Path:
    """Get the output file path."""
    return OUTPUT_DIR / filename


def get_input_path(filename: str) -> Path:
    """Get the input file path."""
    return INPUT_DIR / filename
