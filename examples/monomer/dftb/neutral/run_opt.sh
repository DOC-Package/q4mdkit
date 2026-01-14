#!/bin/bash
# Geometry optimization for neutral pentacene

cp dftb_in_opt.hsd dftb_in.hsd
dftb+ > opt.log 2>&1

echo "Optimization completed. Check opt.log for details."
