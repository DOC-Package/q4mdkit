#!/bin/bash

# Session name
SESSION_NAME="ash_nve"

# Command to run (using conda environment ash)
# Output log to minimization_screen.log
COMMAND="python nve.py"

# Start screen session in detached mode
# Add read to prevent immediate closure after completion
screen -dmS $SESSION_NAME -L -Logfile "nve.log" bash -c "$COMMAND; echo '--------------------------------'; echo 'Job finished. Press Enter to close.'; read"

echo "=================================================="
echo "Job started in screen session: $SESSION_NAME"
echo "Log file: nve.log"
echo "=================================================="
