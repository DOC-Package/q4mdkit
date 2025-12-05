#!/usr/bin/env python3
"""Check if SCCDFTB is available in CHARMM build."""
import pycharmm
from pycharmm import lingo

# Print CHARMM version banner
lingo.charmm_script("system \"echo '=== CHARMM Build Info ==='\"")

# Check for SCCDFTB support by trying to get help
print("\n=== Testing SCCDFTB command ===")
try:
    lingo.charmm_script("HELP SCCDFTB")
    print("SCCDFTB help found")
except:
    print("SCCDFTB command not recognized")

# Check QUANTUM command options
print("\n=== Testing QUANTUM command ===")
try:
    lingo.charmm_script("HELP QUANTUM")
    print("QUANTUM help retrieved")
except:
    print("QUANTUM help not available")

# Try to check what QM methods are compiled in
print("\n=== Available QM methods ===")
test_methods = ["DFTB", "SCCDFTB", "GAMESS", "MOPAC", "ORCA", "QCHEM", "CADPAC"]
for method in test_methods:
    try:
        # Just check if the keyword is recognized (will fail on missing selection, but that's OK)
        lingo.charmm_script(f"BOMLEV -5")
        result = lingo.charmm_script(f"QUANTUM {method} SELE NONE END")
        print(f"  {method}: potentially available")
    except Exception as e:
        if "Unrecognized" in str(e) or "not specified" in str(e):
            print(f"  {method}: NOT available")
        else:
            print(f"  {method}: check inconclusive ({str(e)[:50]})")

lingo.charmm_script("BOMLEV 0")
print("\nCheck complete.")
