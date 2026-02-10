#!/bin/bash
# Geometry optimization for neutral pentacene

cp dftb_in_opt.hsd dftb_in.hsd
screen -S dftb_opt -L -Logfile opt.log dftb+

echo "Optimization completed. Check opt.log for details."
