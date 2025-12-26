#!/usr/bin/env python3
"""
Run CDFTB analysis with online CDFTB-CI calculation.

This script performs:
1. CDFTB calculations for each fragment
2. CDFTB-CI calculation immediately after (on-the-fly)

Unlike the separate cdftb.py + cdftbci.py workflow:
- Uses a single working directory (reused for each frame)
- No H/S matrix storage needed
- Significantly reduced disk usage

Usage:
    python run_cdftb_with_ci.py [config_file]
    
Default config file: cdftb_with_ci_settings.yaml
"""

from pathlib import Path
import sys

# Add parent directory to path for module imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from qm4d4crystal.analysis.cdftb_with_ci import run_cdftb_with_ci


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run CDFTB with online CDFTB-CI calculation"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftb_with_ci_settings.yaml",
        help="Path to YAML configuration file (default: cdftb_with_ci_settings.yaml)"
    )
    args = parser.parse_args()
    
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        return 1
    
    run_cdftb_with_ci(config_path)
    return 0


if __name__ == "__main__":
    exit(main())
