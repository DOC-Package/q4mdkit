#!/bin/bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Session name
SESSION_NAME="vg30_qmmm_npt"

# Command to run (using conda environment ash)
# Output log to minimization_screen.log
COMMAND="python3 npt.py"

# Start screen session in detached mode
# Add read to prevent immediate closure after completion
screen -dmS $SESSION_NAME -L -Logfile "npt.log" bash -c "$COMMAND; echo '--------------------------------'; echo 'Job finished. Press Enter to close.'; read"

echo "=================================================="
echo "Job started in screen session: $SESSION_NAME"
echo "Log file: npt.log"
echo "=================================================="
