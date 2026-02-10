#!/bin/bash
# Normal mode calculation for neutral pentacene
# Run this after geometry optimization is complete

# Post-processing with modes
screen -S modes -L -Logfile modes.log modes

echo "Vibrational analysis completed. Check vib.log and modes.log for details."
