"""
Compatibility subpackage for `q4md4crystal.qmmm` under `q4mdkit.qmmm`.

This module re-exports the submodules from the original package so that
code using `from q4mdkit.qmmm import ...` or
`from q4mdkit.qmmm.qmmm_config import ...` will continue to work.
"""

from importlib import import_module

# Load the original package and use its __path__ so submodule imports
# like `q4mdkit.qmmm.qmmm_config` are resolved by the import system.
_orig_pkg = import_module('qm4d4crystal.qmmm')
__path__ = getattr(_orig_pkg, '__path__', [])

__all__ = [name for name in dir(_orig_pkg) if not name.startswith('_')]
for _name in __all__:
    globals()[_name] = getattr(_orig_pkg, _name)
