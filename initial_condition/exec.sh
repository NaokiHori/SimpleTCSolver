#!/bin/bash

## where to write resulting NPY files
export dirname="output"

## domain
# domain lengths
export lx=1.0e+0
export ly=1.5707963267948966
export lz=2.0e+0
# number of cell centers
export glisize=32
export gljsize=16
export glksize=32
# grid is uniformly distributed in x (true) or not (false)
export uniformx=false

./a.out
