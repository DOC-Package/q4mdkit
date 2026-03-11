#!/bin/bash
# Normal mode calculation for neutral pentacene
# Run this after geometry optimization is complete

cp dftb_in_vib.hsd dftb_in.hsd
screen -S dftb_vib -L -Logfile vib.log dftb+

# Post-processing with modes
screen -S modes -L -Logfile modes.log modes

echo "Vibrational analysis completed. Check vib.log and modes.log for details."
