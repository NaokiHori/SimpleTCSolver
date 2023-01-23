#!/bin/bash

# restart simulation or not
export restart_sim=false
# if restart_sim=true, restart_dir is mandatory
# normally information exist under "output/save/stepxxxxxxxxxx"
# export restart_dir=output/save/stepxxxxxxxxxx

## durations
# maximum duration (in free-fall time)
export timemax=1.0e+3
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=2.0e-1
# logging after (in free-fall time)
export log_after=0.0e+0
# save rate (in free-fall time)
export save_rate=5.0e+1
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=1.0e+2

## domain
# domain lengths
# distance between inner and outer cylinders
# NOTE: inner cylinder radius is fixed to 1
export lx=1.0e+0
# azimuthal length (full cylinder: 6.28...)
export ly=1.5707963267948966
# axial length
export lz=2.0e+0
# number of cell centers
export glisize=32
export gljsize=16
export glksize=32
# gather grid points close to the walls in x direction or not
export use_stretched_grid=true

## treatment of diffusive terms
export implicitx=true
export implicity=false
export implicitz=false

## safety factors to decide time step size
## for advective and diffusive terms
export coef_dt_adv=0.95
export coef_dt_dif=0.95

## physical parameters
export Re=80.

mpirun -n 2 --oversubscribe ./a.out
