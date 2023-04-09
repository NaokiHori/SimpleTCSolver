#include <stdio.h>
#include <math.h>
#include <float.h>
#include "array.h"
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "save.h"
#include "fileio.h"
#include "config.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"

static double g_rate = 0.;
static double g_next = 0.;

/**
 * @brief constructor - schedule saving flow fields
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 * @param[in] rate   : output rate
 * @param[in] after  : information is saved after this time
 */
static void init(const domain_t *domain, const double time, const double rate, const double after){
  const double next = rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / rate
  );
  g_rate = rate;
  g_next = next;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  if(myrank == 0){
    printf("save       initialised, next output : % .3e\n", g_next);
  }
}

/**
 * @brief destructor
 */
static void finalise(void){
}

/**
 * @brief output flow fields etc. to files
 *          which contain essential information to restart
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : time step
 * @param[in] time   : current time units
 * @param[in] fluid  : velocity and pressure
 */
static void output(const domain_t *domain, const int step, const double time, const fluid_t *fluid){
  // get communicator to identify the main process
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  // create directory from main process
  const char prefix[] = {"output/save/step"};
  char *dirname = fileio_create_dirname(prefix, step);
  if(NULL == dirname){
    return;
  }
  if(0 == myrank){
    // create directory
    // this should be invoked only by the main process
    fileio_mkdir(dirname);
  }
  // other processes wait for the completion
  MPI_Barrier(MPI_COMM_WORLD);
  // save current time step and time units
  if(0 == myrank){
    fileio_w_serial(dirname, "step", 0, NULL, NPY_INT, sizeof   (int), &step);
    fileio_w_serial(dirname, "time", 0, NULL, NPY_DBL, sizeof(double), &time);
    config.output(dirname);
  }
  // save coordinates
  if(0 == myrank){
    domain_save(dirname, domain);
  }
  // ux
  array_dump(comm_cart, dirname, "ux", fluid->ux);
  // uy
  array_dump(comm_cart, dirname, "uy", fluid->uy);
  // uz
  array_dump(comm_cart, dirname, "uz", fluid->uz);
  // p
  array_dump(comm_cart, dirname, "p", fluid->p);
  common_free(dirname);
  // schedule next saving event
  g_next += g_rate;
}

/**
 * @brief getter of a member: g_next
 * @return : g_next
 */
static double get_next_time(void){
  return g_next;
}

const save_t save = {
  .init          = init,
  .finalise      = finalise,
  .output        = output,
  .get_next_time = get_next_time,
};

