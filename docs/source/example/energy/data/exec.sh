#!/bin/bash

set -x

# set initial condition, model 1/4 cylinder
cd initial_conditions
make all && make output
dirname=output \
lx=1.0e+0 \
ly=0.7853981633974483 \
lz=2.0e+0 \
glisize=32 \
gljsize=16 \
glksize=32 \
uniformx=false \
./a.out
# go back to root
cd ..

# pre process
make all && make output

# create seed
dirname_ic=initial_conditions/output \
timemax=1.0e+2 \
wtimemax=6.0e+2 \
log_rate=1.0e+0 \
save_rate=2.0e+2 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=1.0e+2 \
implicitx=false \
implicity=false \
implicitz=false \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Re=1.e8 \
mpirun -n 2 --oversubscribe ./a.out

# stash last flow field
dirname_ic=$(find output/save -type d | sort | tail -n 1)
mv ${dirname_ic} dirname_ic

# try different dt
for factor in 0.1 0.2 0.4 0.8
do
  make datadel
  dirname_ic=dirname_ic \
  timemax=105.   \
  wtimemax=6.0e+2 \
  log_rate=1.0e-1 \
  save_rate=5.0e+1 \
  save_after=0.0e+0 \
  stat_rate=1.0e-1 \
  stat_after=1.0e+3 \
  implicitx=false \
  implicity=false \
  implicitz=false \
  coef_dt_adv=${factor} \
  coef_dt_dif=0.95 \
  Re=1.e16 \
  mpirun -n 2 --oversubscribe ./a.out
  mv output/log/energy.dat ./energy-${factor}.dat
done

# post process
mkdir artifacts
python \
  docs/source/example/energy/data/energy.py \
  . \
  artifacts/energy1.png \
  artifacts/energy2.png
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci.txt
echo "Date :" $(date) >> artifacts/ci.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci.txt

set +x
