#include <stdio.h>
#include <string.h>
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

// parameters deciding directory name
static const char dirname_prefix[] = {"output/save/step"};
static const int dirname_ndigits = 10;

// name of directory
static char * restrict g_dirname = NULL;
static size_t g_dirname_nchars = 0;

static double g_rate = 0.;
static double g_next = 0.;

/**
 * @brief constructor - schedule saving flow fields
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 */
static int init(const domain_t * restrict domain, const double time){
  // fetch timings
  if(0 != config.get_double("save_rate", &g_rate)){
    return 1;
  }
  double after = 0.;
  if(0 != config.get_double("save_after", &after)){
    return 1;
  }
  // schedule next event
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / g_rate
  );
  // allocate directory name
  g_dirname_nchars =
    + strlen(dirname_prefix)
    + dirname_ndigits;
  g_dirname = common_calloc(g_dirname_nchars + 2, sizeof(char));
  // report
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(myrank == 0){
    printf("SAVE\n");
    printf("\tnext: % .3e\n", g_next);
    printf("\trate: % .3e\n", g_rate);
    fflush(stdout);
  }
  return 0;
}

/**
 * @brief destructor
 */
static void finalise(void){
  common_free(g_dirname);
}

/**
 * @brief prepare place to output flow fields and save auxiliary data
 * @param[in]  domain  : information related to MPI domain decomposition
 * @param[in]  step    : time step
 * @param[in]  time    : current time units
 * @param[out] dirname : name of created directory
 */
static int prepare(const domain_t * restrict domain, const int step, const double time, char * restrict * dirname){
  // set directory name
  snprintf(g_dirname, g_dirname_nchars + 1, "%s%0*d", dirname_prefix, dirname_ndigits, step);
  *dirname = g_dirname;
  // get communicator to identify the main process
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory and save time step / time units
  if(0 == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio_mkdir(g_dirname);
    // save current time step and time units
    fileio_w_serial(g_dirname, "step", 0, NULL, NPY_INT, sizeof(   int), &step);
    fileio_w_serial(g_dirname, "time", 0, NULL, NPY_DBL, sizeof(double), &time);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // schedule next saving event
  g_next += g_rate;
  return 0;
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
  .prepare       = prepare,
  .get_next_time = get_next_time,
};

