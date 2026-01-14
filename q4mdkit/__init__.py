"""
Compatibility alias package `q4mdkit` for `qm4d4crystal`.
This allows `import q4mdkit` while the original package directory
`qm4d4crystal/` remains in place during the rename transition.

The intent is to fully rename the package later by moving files,
but this alias is a safe first step.
"""

from importlib import import_module

# Import the original package and re-export its public symbols
_orig = import_module('qm4d4crystal')

# Re-export everything that doesn't start with '_' from the original package
__all__ = [name for name in dir(_orig) if not name.startswith('_')]
for _name in __all__:
    globals()[_name] = getattr(_orig, _name)

# Keep a reference to the original module
__original_module__ = _orig
