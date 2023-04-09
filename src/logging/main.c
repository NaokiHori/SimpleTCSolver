#include <stdio.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "logging.h"
#include "internal.h"

static double g_rate = 0.;
static double g_next = 0.;

static int get_ndigits(int num){
  /*
   * E.g., num =    3 -> return 1
   * E.g., num =   13 -> return 2
   * E.g., num = 1234 -> return 4
   * N.B. negative values are not considered
   */
  if(num < 0) return 0;
  int retval = 1;
  while(num /= 10){
    retval++;
  }
  return retval;
}

/**
 * @brief constructor - schedule logging
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 * @param[in] rate   : output rate
 */
static void init(const domain_t *domain, const double time, const double rate){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const double next = rate * ceil(
      fmax(DBL_EPSILON, time) / rate
  );
  g_rate = rate;
  g_next = next;
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    printf("logging    initialised, next logging: % .3e\n", g_next);
  }
}

/**
 * @brief destructor
 */
static void finalise(void){
}

/**
 * @brief show current step, time, time step size, diffusive treatments
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fname  : file name to which the log is written
 * @param[in] time   : current simulation time
 * @param[in] step   : current time step
 * @param[in] dt     : time step size
 * @param[in] wtime  : current wall time
 */
static void show_progress(const char fname[], const domain_t *domain, const double time, const int step, const double dt, const double wtime){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      const double  timemax = config.get. timemax();
      const double wtimemax = config.get.wtimemax();
      // compute number of digits, just for pretty print
      const int ndigits_step     = step < pow(10, 9) ? 10 : 16;
      const int ndigits_timemax  = get_ndigits((int) timemax) + 2;
      const int ndigits_wtimemax = get_ndigits((int)wtimemax) + 2;
      // show progress to standard output and file
      // output to stdout and file
#define MPRINT(...) { \
      fprintf(fp,     __VA_ARGS__); \
      fprintf(stdout, __VA_ARGS__); \
}
      MPRINT("step %*d, time %*.1f, dt %.2e, elapsed %*.1f [sec]\n",
          ndigits_step, step,
          ndigits_timemax, time,
          dt,
          ndigits_wtimemax, wtime
      );
#undef MPRINT
      fileio_fclose(fp);
    }
  }
}

/**
 * @brief output log files to be monitored during simulation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @param[in] time   : current simulation time
 * @param[in] dt     : time step size
 * @param[in] wtime  : current wall time
 * @param[in] fluid  : velocity
 */
static void check_and_output(const domain_t *domain, const int step, const double time, const double dt, const double wtime, const fluid_t *fluid){
  show_progress           ("output/log/progress.dat",   domain, time, step, dt, wtime);
  logging_check_divergence("output/log/divergence.dat", domain, time, fluid);
  logging_check_momentum  ("output/log/momentum.dat",   domain, time, fluid);
  logging_check_energy    ("output/log/energy.dat",     domain, time, fluid);
  logging_check_nusselt   ("output/log/nusselt.dat",    domain, time, fluid);
  /* logging_check_nusselt   ("output/log/nusselt.dat",    domain, time, fluid); */
  g_next += g_rate;
}

/**
 * @brief getter of a member: g_next
 * @return : g_next
 */
static double get_next_time(void){
  return g_next;
}

const logging_t logging = {
  .init             = init,
  .finalise         = finalise,
  .check_and_output = check_and_output,
  .get_next_time    = get_next_time,
};

