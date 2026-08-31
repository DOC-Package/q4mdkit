#!/bin/bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SESSION_NAME="vg30_qmmm_nvt"
#PYTHON="/home/takahashi/miniforge3/envs/q4mdkit/bin/python3"
#COMMAND="$PYTHON nvt.py"
COMMAND="python3 nvt.py"

screen -dmS "$SESSION_NAME" -L -Logfile "nvt.log" bash -c "$COMMAND; echo '--------------------------------'; echo 'Job finished. Press Enter to close.'; read"

echo "=================================================="
echo "Job started in screen session: $SESSION_NAME"
echo "Log file: nvt.log"
echo "=================================================="
