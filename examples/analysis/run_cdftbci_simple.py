#!/usr/bin/env python3

import argparse
import os
import signal
import sys
from pathlib import Path

from q4mdkit.analysis.cdftbci_simple import run_cdftbci_simple_analysis


def _install_group_cleanup() -> None:
    try:
        os.setpgrp()
    except OSError:
        return

    def _graceful_exit(signum, _frame):
        sys.exit(128 + signum)

    signal.signal(signal.SIGINT, _graceful_exit)
    signal.signal(signal.SIGTERM, _graceful_exit)


def _kill_process_group() -> None:
    pgid = os.getpgrp()
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

    import time

    time.sleep(1.0)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except Exception:
        pass


def main() -> None:
    _install_group_cleanup()

    parser = argparse.ArgumentParser(
        description="Run minimal CDFTB-CI analysis with ODIN-based phase correction"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="cdftbci_simple.yaml",
        help="Path to YAML configuration file (default: cdftbci_simple.yaml)",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)

    print(f"Running simple CDFTB-CI analysis with config: {config_path}")
    try:
        run_cdftbci_simple_analysis(config_path)
    finally:
        _kill_process_group()


if __name__ == "__main__":
    main()