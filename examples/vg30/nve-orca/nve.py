#!/usr/bin/env python3
"""Continue the ORCA NVT final state with constant-energy QM/MM dynamics."""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from orca_md import run

if __name__ == "__main__":
    run("nve", SCRIPT_DIR)
