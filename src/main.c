#include <stdio.h>
#include <mpi.h>
#include "common.h"
#include "param.h"
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
 * @param[in] argc : number of arguments (expect 1)
 * @param[in] argv : name of the directory
 *                     where an initial condition is contained
 * @return         : error code
 */
int main(int argc, char *argv[]){
  // launch MPI, start timer
  MPI_Init(NULL, NULL);
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  double wtimes[2] = {0., 0.};
  wtimes[0] = common_get_wtime();
  // find name of directory where IC is stored
  if(2 != argc){
    if(0 == myrank) printf("directory name should be given as input\n");
    goto abort;
  }
  const char * restrict dirname_ic = argv[1];
  if(0 == myrank) printf("load from: %s\n", dirname_ic);
  // initialise time step and time units
  int step = 0;
  double time = 0.;
  if(0 != fileio_r_serial(dirname_ic, "step", 0, NULL, NPY_INT, sizeof(   int), &step)) goto abort;
  if(0 != fileio_r_serial(dirname_ic, "time", 0, NULL, NPY_DBL, sizeof(double), &time)) goto abort;
  if(0 == myrank) printf("step: %d, time: % .7e\n", step, time);
  // initialise structures
  domain_t * restrict domain = NULL;
  fluid_t  * restrict  fluid = NULL;
  if(0 != domain_init(dirname_ic, &domain))        goto abort;
  if(0 !=  fluid_init(dirname_ic, domain, &fluid)) goto abort;
  // initialise auxiliary objects
  if(0 !=    logging.init(domain, time)) goto abort;
  if(0 !=       save.init(domain, time)) goto abort;
  if(0 != statistics.init(domain, time)) goto abort;
  // check termination conditions
  double  timemax = 0.;
  double wtimemax = 0.;
  if(0 != config.get_double("timemax", &timemax))   goto abort;
  if(0 != config.get_double("wtimemax", &wtimemax)) goto abort;
  if(0 == myrank) printf("timemax: % .7e, wtimemax: % .7e\n", timemax, wtimemax);
  /* main loop */
  for(;;){
    // decide time step size
    double dt = 0.;
    if(0 != decide_dt(domain, fluid, &dt)) goto abort;
    // Runge-Kutta iterations
    // max iteration, should be three
    const int rkstepmax = sizeof(param_rkcoefs) / sizeof(rkcoef_t);
    for(int rkstep = 0; rkstep < rkstepmax; rkstep++){
      // update velocity
      fluid_compute_rhs(domain, rkstep, fluid);
      fluid_update_velocity(domain, rkstep, dt, fluid);
      // compute scalar potential
      if(0 != fluid_compute_potential(domain, rkstep, dt, fluid)) goto abort;
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
    // the simulation is finished
    if(timemax < time) break;
    // wall time limit is reached
    if(wtimemax < wtimes[1] - wtimes[0]) break;
    // compute and output log regularly
    if(logging.get_next_time() < time){
      logging.check_and_output(domain, step, time, dt, wtimes[1] - wtimes[0], fluid);
    }
    // save flow fields regularly
    if(save.get_next_time() < time){
      char * restrict dirname = NULL;
      save.prepare(domain, step, time, &dirname);
      domain_save(dirname, domain);
      fluid_save(dirname, domain, fluid);
    }
    // collect statistics regularly
    if(statistics.get_next_time() < time){
      statistics.collect(domain, fluid);
    }
  }
  // finalisation
  // save final flow fields
  {
    char * restrict dirname = NULL;
    save.prepare(domain, step, time, &dirname);
    domain_save(dirname, domain);
    fluid_save(dirname, domain, fluid);
  }
  // save collected statistics
  statistics.output(domain, step, time);
  // finalise structures
  fluid_finalise(fluid);
  domain_finalise(domain);
  logging.finalise();
  save.finalise();
  statistics.finalise();
  // finalise MPI
abort:
  MPI_Finalize();
  return 0;
}

