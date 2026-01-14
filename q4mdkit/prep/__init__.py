"""
Compatibility subpackage for `qm4d4crystal.prep` under `q4mdkit.prep`.

Re-exports the original `qm4d4crystal.prep` symbols.
"""

from importlib import import_module

_orig_pkg = import_module('qm4d4crystal.prep')

__all__ = [name for name in dir(_orig_pkg) if not name.startswith('_')]
for _name in __all__:
    globals()[_name] = getattr(_orig_pkg, _name)
