#!/usr/bin/env python3
"""Quick test of run_md.py with reduced steps"""
import sys
sys.path.insert(0, '/home/takahashi/python/mdcrystal/src')

from pathlib import Path
exec(compile(Path('/home/takahashi/python/mdcrystal/example/charmm/pentacene/run_md.py').read_text(), 'run_md.py', 'exec'))
