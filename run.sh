#!/bin/bash
# Helper script to run mdcrystal tools from any directory

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Add src directory to PYTHONPATH
export PYTHONPATH="${SCRIPT_DIR}/src:${PYTHONPATH}"

# Run the specified Python script with all arguments
python "$@"
