"""Canonical paths shared by the VG30 preparation and MD stages."""

from pathlib import Path


VG30_DIR = Path(__file__).resolve().parent
INPUT_DIR = VG30_DIR / "input"
OPT_MM_DIR = VG30_DIR / "opt-mm"

PACKED_PDB = INPUT_DIR / "vg30_acetonitrile.pdb"
AMBER_PRMTOP = INPUT_DIR / "vg30_acetonitrile.prmtop"
BOX_FILE = INPUT_DIR / "vg30_acetonitrile.box"
MINIMIZED_PDB = OPT_MM_DIR / "frag-minimized.pdb"
