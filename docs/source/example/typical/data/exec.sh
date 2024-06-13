#!/bin/bash

set -e
set -x

# set initial condition, model half cylinder
cd initial_condition
make output && make datadel
is_curved=true \
lx=1.0e+0 \
ly=0.39269908169872414 \
lz=2.0e+0 \
glisize=64 \
gljsize=8 \
glksize=32 \
python3 main.py output
cd ..

# build and run
make output && make datadel && make all
timemax=4.0e+2 \
wtimemax=6.0e+2 \
log_rate=1.0e+0 \
save_rate=1.0e+2 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=1.0e+2 \
coef_dt_adv=1. \
coef_dt_dif=1. \
Re=80. \
Pr=1. \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts
python3 \
  docs/source/example/typical/data/snapshot.py \
  $(find output/save -type d | sort | tail -n 1) \
  artifacts/snapshot.png
python3 \
  docs/source/example/typical/data/divergence.py \
  output/log/divergence.dat \
  artifacts/divergence.png
python3 \
  docs/source/example/typical/data/balance_main.py \
  $(find output/save -type d | sort | tail -n 1) \
  output/log/dissipation.dat \
  output/log/injection.dat \
  artifacts/balance_main.png
python3 \
  docs/source/example/typical/data/balance_dif.py \
  output/log/dissipation.dat \
  output/log/injection.dat \
  artifacts/balance_dif.png

