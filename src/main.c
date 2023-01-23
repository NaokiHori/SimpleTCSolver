#include <stdbool.h>
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
#include "tdm.h"


/**
 * @brief main function
 * @return : error code
 */
int main(void){
  /* ! launch MPI, start timer ! 3 ! */
  MPI_Init(NULL, NULL);
  double wtimes[2] = {0., 0.};
  wtimes[0] = common_get_wtime();
  // load environmental variables
  config.load();
  /* ! initialise time step and current (simulation) time ! 15 ! */
  int step;
  double time;
  {
    const bool restart_sim = config.get_bool("restart_sim");
    if(restart_sim){
      // directory name from which information are loaded
      const char *dirname = config.get_string("restart_dir");
      fileio_r_0d_serial(dirname, "step", NPYIO_INT,    sizeof(   int), &step);
      fileio_r_0d_serial(dirname, "time", NPYIO_DOUBLE, sizeof(double), &time);
    }else{
      // start from 0
      step = 0;
      time = 0.;
    }
  }
  /* ! initialise structures ! 3 ! */
  domain_t     *domain     = domain_init();
  fluid_t      *fluid      = fluid_init(domain);
  statistics_t *statistics = statistics_init(domain);
  /* main loop */
  for(;;){
    /* ! decide time step size and diffusive term treatment ! 1 ! */
    const double dt = decide_dt(domain, fluid);
    /* ! integrate mass, momentum, and internal energy balances in time ! 17 ! */
    for(int rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
      /* ! update velocity by integrating momentum equation ! 2 ! */
      fluid_compute_rhs(domain, rkstep, fluid);
      fluid_update_velocity(domain, rkstep, dt, fluid);
      /* ! compute scalar potential ! 1 ! */
      fluid_compute_potential(domain, rkstep, dt, fluid);
      /* ! correct velocity to be solenoidal ! 1 ! */
      fluid_correct_velocity(domain, rkstep, dt, fluid);
      /* ! update pressure ! 1 ! */
      fluid_update_pressure(domain, rkstep, dt, fluid);
    }
    /* ! step and time are incremented ! 3 ! */
    step += 1;
    time += dt;
    wtimes[1] = common_get_wtime();
    /* ! output log ! 3 ! */
    if(logging_next < time){
      logging(domain, step, time, dt, wtimes[1]-wtimes[0], fluid);
    }
    /* ! save flow fields ! 3 ! */
    if(save_next < time){
      save(domain, step, time, fluid);
    }
    /* ! collect statistics ! 3 ! */
    if(stat_next < time){
      statistics_collect(domain, time, fluid, statistics);
    }
    /* ! terminate when the simulation is finished ! 3 ! */
    if(time > config.get_double("timemax")){
      break;
    }
    /* ! terminate when wall time limit is reached ! 3 ! */
    if(wtimes[1]-wtimes[0] > config.get_double("wtimemax")){
      break;
    }
  }
  /* ! save restart file and statistics at last ! 2 ! */
  save(domain, step, time, fluid);
  statistics_output(domain, step, time, statistics);
  /* ! finalise structures ! 3 ! */
  statistics_finalise(statistics);
  fluid_finalise(fluid);
  domain_finalise(domain);
  // free internal memories used by tdm
  tdm_cleanup();
  // free memory used to store environmental variables
  config.unload();
  /* ! finalise MPI ! 1 ! */
  MPI_Finalize();
  return 0;
}

