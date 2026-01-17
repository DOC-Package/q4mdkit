"""
MD Configuration Module for MM simulations

This is a symbolic link / re-export of the qmmm.md_config module.
The MD settings (NVT, NPT, NVE) are the same regardless of whether
we use QM/MM or pure MM.
"""

# Re-export everything from qmmm.md_config
from q4mdkit.qmmm.md_config import MDConfig, get_md_config

__all__ = ['MDConfig', 'get_md_config']
