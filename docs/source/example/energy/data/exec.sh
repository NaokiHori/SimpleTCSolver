#!/bin/bash

set -e
set -x

# set initial condition, model 1/4 cylinder
cd initial_condition
make output && make datadel
is_curved=true \
lx=1.0e+0 \
ly=0.7853981633974483 \
lz=2.0e+0 \
glisize=32 \
gljsize=16 \
glksize=32 \
python main.py output
cd ..

# build and run
make output && make datadel && make all
timemax=10. \
wtimemax=6.0e+2 \
log_rate=1.0e+0 \
save_rate=2.0e+2 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=1.0e+2 \
coef_dt_adv=1. \
coef_dt_dif=1. \
Re=1.e16 \
Pr=1.e0 \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# stash last flow field
mv $(find output/save -type d | sort | tail -n 1) dirname_ic

# try different dt
for factor in 0.1 0.2 0.4 0.8
do
  make datadel
  timemax=15. \
  wtimemax=6.0e+2 \
  log_rate=1.0e-1 \
  save_rate=5.0e+1 \
  save_after=0.0e+0 \
  stat_rate=1.0e-1 \
  stat_after=1.0e+3 \
  coef_dt_adv=${factor} \
  coef_dt_dif=0.95 \
  Re=1.e16 \
  Pr=1.e0 \
  mpirun -n 2 --oversubscribe ./a.out dirname_ic
  mv output/log/energy.dat ./energy-${factor}.dat
done

# post process
mkdir artifacts
python3 \
  docs/source/example/energy/data/energy.py \
  . \
  artifacts/energy1.png \
  artifacts/energy2.png

