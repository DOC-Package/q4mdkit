#!/usr/bin/env python3
"""
CDFTB-CI analysis script.

Perform CDFTB calculations with online CDFTB-CI calculation for each frame.

This script performs:
1. CDFTB calculations for each fragment (constrained DFT-B)
2. CDFTB-CI calculation immediately after (on-the-fly)
3. Spin population analysis

Features:
- Uses a single working directory (reused for each frame)
- No H/S matrix storage needed  
- Significantly reduced disk usage
- Optional charge continuation between frames

Usage:
    python run_cdftbci.py [config.yaml]

Default config file: cdftb_settings.yaml
"""

import argparse
import os
import signal
import sys
from pathlib import Path
from q4mdkit.analysis.cdftbci import run_cdftbci_analysis


def _install_group_cleanup():
    """Put this process in its own process group and ensure all descendants
    (forkserver, DFTB+ worker subprocesses, OpenMP threads owned by them) are
    terminated when the main script exits — whether via Ctrl-C, SIGTERM, or
    normal completion. Without this, multiprocessing's forkserver children can
    survive as orphans running DFTB+ at 99% CPU."""
    try:
        os.setpgrp()  # new process group; all children inherit it
    except OSError:
        return

    def _graceful_exit(signum, _frame):
        # sys.exit triggers the finally-block below, which kills the group
        sys.exit(128 + signum)

    signal.signal(signal.SIGINT, _graceful_exit)
    signal.signal(signal.SIGTERM, _graceful_exit)


def _kill_process_group():
    """Send SIGTERM to the whole process group, then SIGKILL after a short
    grace period for anything that ignored SIGTERM."""
    pgid = os.getpgrp()
    # Don't kill ourselves first; ignore the signal in this process
    try:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    except Exception:
        pass
    try:
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except Exception:
        pass
    # Brief grace, then force-kill stragglers (DFTB+ rarely honors SIGTERM mid-SCC)
    import time
    time.sleep(1.0)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except Exception:
        pass


def main():
    _install_group_cleanup()

    parser = argparse.ArgumentParser(
        description="Run CDFTB with online CDFTB-CI calculation"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftb_settings.yaml",
        help="Path to YAML configuration file (default: cdftb_settings.yaml)"
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)

    print(f"Running CDFTB-CI analysis with config: {config_path}")
    try:
        run_cdftbci_analysis(config_path)
    finally:
        _kill_process_group()


if __name__ == "__main__":
    main()
