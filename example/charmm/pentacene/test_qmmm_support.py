#!/usr/bin/env python3
"""
Test if CHARMM QM/MM capabilities are available.
"""
try:
    import pycharmm
    from pycharmm import lingo
    cs = lingo.charmm_script
except:
    from pycharmm import charmm_script as cs

print("Testing CHARMM QM/MM capabilities...")
print("="*60)

# Check CHARMM version and available features
cs("SYSTEM \"echo 'CHARMM version info:'\"")

# Try to get help on QUANTUM command
print("\nChecking QUANTUM command availability...")
try:
    cs("HELP QUANTUM")
    print("✓ QUANTUM command is available")
except Exception as e:
    print(f"✗ QUANTUM command may not be available: {e}")

print("\nNote: Full QM/MM functionality requires:")
print("  - CHARMM compiled with QM interface support")
print("  - External QM program (MOPAC, GAMESS, ORCA, etc.) installed")
print("  - Proper environment variables set for QM program")
print("="*60)
