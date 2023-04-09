#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "statistics.h"
#include "save.h"
#include "logging.h"
#include "decide_dt.h"
#include "config.h"
#include "fileio.h"

/**
 * @brief main function
 * @return : error code
 */
int main(void){
  // launch MPI, start timer 
  MPI_Init(NULL, NULL);
  double wtimes[2] = {0., 0.};
  wtimes[0] = common_get_wtime();
  // load environment variables 
  // abort when variable is missing
  if(0 != config.construct()){
    goto abort;
  }
  // find name of directory where IC is stored 
  const char *dirname_ic = config.get.dirname_ic();
  // initialise time step and time units 
  int step = 0;
  double time = 0.;
  if(0 != fileio_r_serial(dirname_ic, "step", 0, NULL, NPY_INT, sizeof(   int), &step)){
    goto abort;
  }
  if(0 != fileio_r_serial(dirname_ic, "time", 0, NULL, NPY_DBL, sizeof(double), &time)){
    goto abort;
  }
  // initialise structures 
  domain_t *domain = domain_init(dirname_ic);
  if(NULL == domain){
    goto abort;
  }
  fluid_t *fluid = fluid_init(dirname_ic, domain);
  if(NULL == fluid){
    goto abort;
  }
  // initialise auxiliary objects 
  logging.init(
      domain,
      time,
      config.get.log_rate()
  );
  save.init(
      domain,
      time,
      config.get.save_rate(),
      config.get.save_after()
  );
  statistics.init(
      domain,
      time,
      config.get.stat_rate(),
      config.get.stat_after()
  );
  /* main loop */
  for(;;){
    // decide time step size
    const double dt = decide_dt(domain, fluid);
    // Runge-Kutta iterations
    for(int rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
      // update velocity 
      fluid_compute_rhs(domain, rkstep, fluid);
      fluid_update_velocity(domain, rkstep, dt, fluid);
      // compute scalar potential 
      const int retval = fluid_compute_potential(domain, rkstep, dt, fluid);
      if(0 != retval) goto abort;
      // correct velocity field to satisfy mass conservation
      fluid_correct_velocity(domain, rkstep, dt, fluid);
      // update pressure
      fluid_update_pressure(domain, rkstep, dt, fluid);
    }
    // update step and simulation / wall time 
    step += 1;
    time += dt;
    wtimes[1] = common_get_wtime();
    // terminate if one of the following conditions is met 
    if(config.get.timemax() < time){
      // the simulation is finished
      break;
    }
    if(config.get.wtimemax() < wtimes[1] - wtimes[0]){
      // wall time limit is reached
      break;
    }
    // compute and output log regularly 
    if(logging.get_next_time() < time){
      logging.check_and_output(domain, step, time, dt, wtimes[1] - wtimes[0], fluid);
    }
    // save flow fields regularly 
    if(save.get_next_time() < time){
      save.output(domain, step, time, fluid);
    }
    // collect statistics regularly 
    if(statistics.get_next_time() < time){
      statistics.collect(domain, fluid);
    }
  }
  // finalisation
  // save flow fields and statistics 
  save.output(domain, step, time, fluid);
  statistics.output(domain, step, time);
  // finalise structures 
  fluid_finalise(fluid);
  domain_finalise(domain);
  logging.finalise();
  save.finalise();
  statistics.finalise();
  config.destruct();
  // finalise MPI 
abort:
  MPI_Finalize();
  return 0;
}

