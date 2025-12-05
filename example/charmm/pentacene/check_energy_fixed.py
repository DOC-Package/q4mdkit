#!/usr/bin/env python3
"""Check energy with proper PBC setup (no COOR ORIENT)"""
import os
try:
    import pycharmm
    from pycharmm import lingo
    cs = lingo.charmm_script
except Exception:
    from pycharmm import charmm_script as cs

def open_unit(unit: int, mode: str, path: str):
    path = os.path.abspath(path)
    cs(f"OPEN {mode.upper()} UNIT {unit:d} CARD NAME {path}")

def find_valid_fft(n):
    if n <= 1: return 2
    candidate = n
    while True:
        temp = candidate
        while temp % 2 == 0: temp //= 2
        while temp % 3 == 0: temp //= 3
        while temp % 5 == 0: temp //= 5
        if temp == 1: return candidate
        candidate += 1

cs("BOMLEV -2")
cs("PRNLEV 5")

# Load force field
base_rtf = "/home/takahashi/charmm/charmm/toppar/top_all36_cgenff.rtf"
base_prm = "/home/takahashi/charmm/charmm/toppar/par_all36_cgenff.prm"
if os.path.exists(base_rtf):
    cs(f"open read unit 1 card name {base_rtf}")
    cs("read rtf card unit 1")
    cs("close unit 1")
    cs(f"open read unit 1 card name {base_prm}")
    cs("read param card unit 1")
    cs("close unit 1")

import tempfile
with open("pentacene.str", 'r') as f:
    content = f.read()
content = content.replace('read param card flex append', 'read param card append')
with tempfile.NamedTemporaryFile(mode='w', suffix='.str', delete=False) as tmp:
    tmp_path = tmp.name
    tmp.write(content)
cs(f"stream {tmp_path}")
os.unlink(tmp_path)

cs("BOMLEV 0")

# Read PSF
open_unit(10, "READ", "pentacene.psf")
cs("READ PSF CARD UNIT 10")

# Read wrapped coordinates WITHOUT COOR ORIENT
print("="*70)
print("Reading pentacene.crd (wrapped) WITHOUT COOR ORIENT")
print("="*70)
open_unit(20, "READ", "pentacene.crd")
cs("READ COOR CARD UNIT 20")
# NO COOR ORIENT!
cs("COOR STAT")

# Read box
with open("pentacene.crd", "rt") as f:
    lines = f.readlines()
last_line = lines[-1].strip().split()
a, b, c, alpha, beta, gamma = map(float, last_line)
print(f"\nBox: {a:.3f} {b:.3f} {c:.3f} Å, {alpha:.2f}° {beta:.2f}° {gamma:.2f}°")

# Setup PBC
cs(f"CRYSTAL DEFINE TRIClinic {a:.3f} {b:.3f} {c:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f}")
cs("CRYSTAL BUILD NOPER 0")
cs("IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELECT ALL END")

# Setup PME
gx = find_valid_fft(max(32, int(round(a))))
gy = find_valid_fft(max(32, int(round(b))))
gz = find_valid_fft(max(32, int(round(c))))
print(f"PME grid: {gx} x {gy} x {gz}")

cs(
    "NBOND ATOM CDIEL SWITCH VATOM VDIST VSWITCH "
    f"CUTNB 14.0 CTOFNB 12.0 CTONNB 10.0 "
    f"EWALD PMEWALD KAPPA 0.34 ORDER 6 "
    f"FFTX {gx:d} FFTY {gy:d} FFTZ {gz:d}"
)

print("\n" + "="*70)
print("INITIAL ENERGY:")
print("="*70)
cs("ENERGY")

print("\nDone.")
