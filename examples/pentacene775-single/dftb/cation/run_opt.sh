#!/bin/bash
# Geometry optimization for pentacene cation

cp dftb_in_opt.hsd dftb_in.hsd
dftb+ > opt.log 2>&1

echo "Optimization completed. Check opt.log for details."
