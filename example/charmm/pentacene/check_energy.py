#!/usr/bin/env python3
"""Check initial energy of the system"""
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

# Suppress warnings
cs("BOMLEV -2")
cs("PRNLEV 5")

# Load base force field
base_rtf = "/home/takahashi/charmm/charmm/toppar/top_all36_cgenff.rtf"
base_prm = "/home/takahashi/charmm/charmm/toppar/par_all36_cgenff.prm"

if os.path.exists(base_rtf) and os.path.exists(base_prm):
    print(f"Loading base CHARMM36 CGenFF...")
    cs(f"open read unit 1 card name {base_rtf}")
    cs("read rtf card unit 1")
    cs("close unit 1")
    cs(f"open read unit 1 card name {base_prm}")
    cs("read param card unit 1")
    cs("close unit 1")

# Load molecule parameters
print("Loading pentacene parameters...")
import tempfile
with open("pentacene.str", 'r') as f:
    content = f.read()
content_fixed = content.replace('read param card flex append', 'read param card append')
with tempfile.NamedTemporaryFile(mode='w', suffix='.str', delete=False) as tmp:
    tmp_path = tmp.name
    tmp.write(content_fixed)
cs(f"stream {tmp_path}")
os.unlink(tmp_path)

cs("BOMLEV 0")

# Read PSF and coordinates
print("\nReading PSF and coordinates...")
open_unit(10, "READ", "pentacene.psf")
cs("READ PSF CARD UNIT 10")

# Test with nowrap coordinates
open_unit(20, "READ", "pentacene.crd")
cs("READ COOR CARD UNIT 20")

# Read box info from CRD
with open("pentacene_nowrap.crd", "rt") as f:
    lines = f.readlines()
last_line = lines[-1].strip().split()
if len(last_line) == 6:
    a, b, c, alpha, beta, gamma = map(float, last_line)
    print(f"\nBox from CRD: {a:.3f} {b:.3f} {c:.3f} Å, {alpha:.2f}° {beta:.2f}° {gamma:.2f}°")
    
    # Setup crystal
    cs(f"CRYSTAL DEFINE TRIClinic {a:.3f} {b:.3f} {c:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f}")
    cs("CRYSTAL BUILD NOPER 0")
    cs("IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELECT ALL END")
    
    # Setup nonbond with PME
    print("\nSetting up nonbond interactions with PME...")
    from math import ceil
    
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
    
    gx = find_valid_fft(max(32, int(ceil(a))))
    gy = find_valid_fft(max(32, int(ceil(b))))
    gz = find_valid_fft(max(32, int(ceil(c))))
    print(f"PME grid: {gx} x {gy} x {gz}")
    
    cs(
        "NBOND ATOM CDIEL SWITCH VATOM VDIST VSWITCH "
        f"CUTNB 14.0 CTOFNB 12.0 CTONNB 10.0 "
        f"EWALD PMEWALD KAPPA 0.34 ORDER 6 "
        f"FFTX {gx:d} FFTY {gy:d} FFTZ {gz:d}"
    )
    
    print("\n" + "="*60)
    print("INITIAL ENERGY (pentacene.crd):")
    print("="*60)
    cs("ENERGY")
    
    # Get coordinate statistics
    print("\n" + "="*60)
    print("COORDINATE STATISTICS:")
    print("="*60)
    cs("COOR STAT")
    
    # Check distances between atoms
    print("\n" + "="*60)
    print("CHECKING SHORT DISTANCES (< 2.0 Å):")
    print("="*60)
    cs("COOR MINDIST CUTOFF 2.0")

print("\n" + "="*60)
print("Now testing with wrapped coordinates...")
print("="*60)

# Test with wrapped coordinates
open_unit(20, "READ", "pentacene.crd")
cs("READ COOR CARD UNIT 20")

with open("pentacene.crd", "rt") as f:
    lines = f.readlines()
last_line = lines[-1].strip().split()
if len(last_line) == 6:
    a, b, c, alpha, beta, gamma = map(float, last_line)
    print(f"\nBox from CRD: {a:.3f} {b:.3f} {c:.3f} Å, {alpha:.2f}° {beta:.2f}° {gamma:.2f}°")
    
    cs(f"CRYSTAL DEFINE TRIClinic {a:.3f} {b:.3f} {c:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f}")
    cs("CRYSTAL BUILD NOPER 0")
    cs("IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELECT ALL END")
    
    cs(
        "NBOND ATOM CDIEL SWITCH VATOM VDIST VSWITCH "
        f"CUTNB 14.0 CTOFNB 12.0 CTONNB 10.0 "
        f"EWALD PMEWALD KAPPA 0.34 ORDER 6 "
        f"FFTX {gx:d} FFTY {gy:d} FFTZ {gz:d}"
    )
    
    print("\n" + "="*60)
    print("INITIAL ENERGY (pentacene.crd - wrapped):")
    print("="*60)
    cs("ENERGY")
    
    print("\n" + "="*60)
    print("COORDINATE STATISTICS:")
    print("="*60)
    cs("COOR STAT")
    
    print("\n" + "="*60)
    print("CHECKING SHORT DISTANCES (< 2.0 Å):")
    print("="*60)
    cs("COOR MINDIST CUTOFF 2.0")

print("\nDone.")
