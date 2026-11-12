#!/bin/bash

# Session name
SESSION_NAME="tddftb"

# Command to run (using conda environment ash)
# Output log to nvesc.log
COMMAND="python run_tddftb.py"

# Start screen session in detached mode
# Add read to prevent immediate closure after completion
screen -dmS $SESSION_NAME -L -Logfile "tddftb.log" bash -c "$COMMAND; echo '--------------------------------'; echo 'Job finished. Press Enter to close.'; read"

echo "=================================================="
echo "Job started in screen session: $SESSION_NAME"
echo "Log file: tddftb.log"
echo "=================================================="
