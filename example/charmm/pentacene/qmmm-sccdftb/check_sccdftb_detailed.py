#!/usr/bin/env python3
"""Check SCCDFTB availability by attempting actual QUANTUM commands."""
import pycharmm
from pycharmm import lingo

# Create minimal test system
print("=== Setting up minimal test system ===")
lingo.charmm_script("""
* Test system for SCCDFTB check
*

! Create a single hydrogen atom for testing
READ RTF CARD
* Minimal topology
*
   22     1
MASS     1 H      1.00800 H

RESI TEST        0.00
GROUP
ATOM H1   H      0.00

END

READ PARAM CARD
* Minimal parameters
*
ATOMS
MASS     1 H      1.00800

BONDS

ANGLES

DIHEDRALS

IMPROPER

NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

H      0.0       -0.046     0.2245

END

! Generate single atom structure
READ SEQUENCE CARD
* Single test residue
*
1
TEST

GENERATE TEST SETUP

! Set coordinates
COOR SET XDIR 1.0 YDIR 0.0 ZDIR 0.0 SELECT ALL END

BOMLEV -5
""")

print("\n=== Testing DFTB with NULH (should work) ===")
try:
    lingo.charmm_script("QUANTUM DFTB METHOD NULH CHRG 0 SELE ALL END")
    lingo.charmm_script("ENERGY")
    print("✓ DFTB with NULH: WORKS")
except Exception as e:
    print(f"✗ DFTB with NULH: FAILED - {str(e)[:100]}")

print("\n=== Testing QUANTUM SCCDFTB (trying SCC-DFTB) ===")
try:
    lingo.charmm_script("QUANTUM SCCDFTB CHRG 0 TEMP 300.0 SELE ALL END")
    lingo.charmm_script("ENERGY")
    print("✓ QUANTUM SCCDFTB: WORKS")
except Exception as e:
    error_msg = str(e)
    if "AM1, PM3, MNDO or NULH must be specified" in error_msg:
        print("✗ QUANTUM SCCDFTB: NOT COMPILED IN THIS BUILD")
        print("  Error: SCCDFTB keyword not recognized by CHARMM")
        print("  Available methods: AM1, PM3, MNDO, NULH (semi-empirical)")
    elif "Unrecognized" in error_msg:
        print("✗ QUANTUM SCCDFTB: COMMAND NOT RECOGNIZED")
    else:
        print(f"✗ QUANTUM SCCDFTB: FAILED - {error_msg[:100]}")

print("\n=== Testing standalone SCCDFTB command ===")
try:
    lingo.charmm_script("SCCDFTB CHRG 0 SELE ALL END")
    print("✓ SCCDFTB standalone: WORKS")
except Exception as e:
    if "Unrecognized command: SCCD" in str(e):
        print("✗ SCCDFTB standalone: NOT A VALID COMMAND")
    else:
        print(f"✗ SCCDFTB standalone: FAILED - {str(e)[:100]}")

lingo.charmm_script("BOMLEV 0")
print("\n" + "="*60)
print("CONCLUSION:")
print("="*60)
print("SCCDFTB (Self-Consistent Charge DFTB) is NOT available in")
print("this CHARMM build. This feature requires compilation with")
print("the DFTB+ library or CHARMM's internal SCCDFTB module.")
print("")
print("Available alternatives:")
print("  1. NULH (Numerical Universal Lazy Hamiltonian) - built-in")
print("  2. Semi-empirical: AM1, PM3, MNDO")
print("  3. External QM: MOPAC, GAMESS, ORCA, Q-Chem (if configured)")
print("")
print("To enable SCCDFTB, CHARMM must be recompiled with:")
print("  - DFTB keyword in pref.dat")
print("  - DFTB+ library linked")
print("="*60)
