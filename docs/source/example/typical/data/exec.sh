#!/bin/bash

set -e
set -x

# set initial condition, model half cylinder
cd initial_condition
make output && make datadel && make all
dirname=output \
lx=1.0e+0 \
ly=1.5707963267948966 \
lz=2.0e+0 \
glisize=32 \
gljsize=16 \
glksize=32 \
uniformx=false \
./a.out
cd ..

# pre process
make output && make datadel && make all

# try high Reynolds number to trigger instability
timemax=1.0e+2 \
wtimemax=6.0e+2 \
log_rate=5.0e-1 \
save_rate=5.0e+1 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=1.0e+2 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Re=160. \
mpirun -n 4 --oversubscribe ./a.out initial_condition/output

# try low Reynolds number to obtain result
dirname_ic=$(find output/save -type d | sort | tail -n 1)
timemax=2.0e+2 \
wtimemax=6.0e+2 \
log_rate=5.0e-1 \
save_rate=5.0e+1 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=1.0e+2 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Re=80. \
mpirun -n 8 --oversubscribe ./a.out ${dirname_ic}

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
  docs/source/example/typical/data/nusselt.py \
  output/log/nusselt.dat \
  artifacts/nusselt.png
mv ${dirname_ic} artifacts/
mv output/log artifacts/
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci.txt
echo "Date :" $(date) >> artifacts/ci.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci.txt

set +x
