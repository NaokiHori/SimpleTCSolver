#!/bin/bash

## overall configurations
# curved channel (true) or planar channel (false)
# NOTE: for curved channels, the inner cylinder radius is always fixed to 1
export is_curved=true

## domain size
# domain lengths
export lx=1.0e+0
export ly=1.0e+0
export lz=2.0e+0
# number of cell centers
export glisize=32
export gljsize=8
export glksize=64

## where to write resulting NPY files
export dirname="output"

python3 main.py ${dirname}
