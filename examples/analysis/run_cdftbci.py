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


DEFAULT_PYTHON = Path("/home/takahashi/anaconda3/envs/mydftbplus/bin/python")


def _has_dftbplus_module() -> bool:
    try:
        import dftbplus  # noqa: F401
    except ModuleNotFoundError:
        return False
    return True


def _reexec_with_dftbplus_python_if_needed() -> None:
    if _has_dftbplus_module():
        return
    if Path(sys.executable).resolve() == DEFAULT_PYTHON.resolve():
        return
    if not DEFAULT_PYTHON.exists():
        return
    os.execv(str(DEFAULT_PYTHON), [str(DEFAULT_PYTHON), *sys.argv])


def _prepend_ld_library_path(path: Path) -> None:
    lib_dir = str(path.parent)
    current = os.environ.get("LD_LIBRARY_PATH", "")
    parts = [p for p in current.split(":") if p]
    if lib_dir not in parts:
        os.environ["LD_LIBRARY_PATH"] = ":".join([lib_dir, *parts])
        os.environ["Q4MDKIT_REEXEC_FOR_DFTB_LIBS"] = "1"


def _prepend_ld_preload(paths) -> None:
    current = os.environ.get("LD_PRELOAD", "")
    parts = [p for p in current.split(":") if p]
    changed = False
    for path in reversed([str(p) for p in paths if Path(p).exists()]):
        if path not in parts:
            parts.insert(0, path)
            changed = True
    if changed:
        os.environ["LD_PRELOAD"] = ":".join(parts)
        os.environ["Q4MDKIT_REEXEC_FOR_DFTB_LIBS"] = "1"


def _configure_dftb_runtime_environment(config_path: Path) -> None:
    """Expose DFTB+ dependent shared libraries before worker processes start."""
    try:
        import yaml
    except ModuleNotFoundError:
        return
    try:
        with config_path.open() as handle:
            data = yaml.safe_load(handle) or {}
    except OSError:
        return
    library_path = data.get("dftb", {}).get("library_path")
    if library_path:
        lib_path = Path(library_path)
        _prepend_ld_library_path(lib_path)
        _prepend_ld_preload(
            [
                Path("/lib/x86_64-linux-gnu/libgomp.so.1"),
                Path("/opt/intel/oneapi/mkl/2025.0/lib/libmkl_gnu_thread.so.2"),
                Path("/opt/intel/oneapi/mkl/2025.0/lib/libmkl_core.so.2"),
            ]
        )


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


def _terminate_process_group():
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


def main():
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

    _configure_dftb_runtime_environment(config_path)
    if os.environ.pop("Q4MDKIT_REEXEC_FOR_DFTB_LIBS", None) == "1":
        os.execv(sys.executable, [sys.executable, *sys.argv])

    _reexec_with_dftbplus_python_if_needed()
    _install_group_cleanup()

    from q4mdkit.analysis.cdftbci import run_cdftbci_analysis

    print(f"Running CDFTB-CI analysis with config: {config_path}")
    try:
        run_cdftbci_analysis(config_path)
    finally:
        _terminate_process_group()


if __name__ == "__main__":
    main()
