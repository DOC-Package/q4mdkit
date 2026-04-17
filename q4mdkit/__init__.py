# q4mdkit - QM/MM tools for molecular crystals

from . import analysis

# Optional modules (may require additional dependencies)
try:
    from . import prep
except ImportError:
    prep = None  # ase not available

try:
    from . import qmmm
except ImportError:
    qmmm = None  # ash not available

try:
    from . import mm
except ImportError:
    mm = None  # ash not available

__version__ = "0.1.0"
