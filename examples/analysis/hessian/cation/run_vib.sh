#!/bin/bash
# Normal mode calculation for neutral pentacene
# Run this after geometry optimization is complete

cp dftb_in_vib.hsd dftb_in.hsd
dftb+ > vib.log 2>&1

# Post-processing with modes
modes > modes.log 2>&1

echo "Vibrational analysis completed. Check vib.log and modes.log for details."
