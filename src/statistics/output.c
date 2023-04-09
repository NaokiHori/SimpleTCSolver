#include "common.h"
#include "config.h"
#include "domain.h"
#include "statistics.h"
#include "fileio.h"
#include "internal.h"
#include "arrays/ux1.h"
#include "arrays/ux2.h"
#include "arrays/uy1.h"
#include "arrays/uy2.h"
#include "arrays/uz1.h"
#include "arrays/uz2.h"

/**
 * @brief save structures which contains collected statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @param[in] time   : current time units
 */
void output(const domain_t *domain, const int step, const double time){
  const int num = g_st_int->num;
  if(0 == num){
    // no reason to save statistics,
    //   since they are not collected
    return;
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  // create directory from main process
  const char prefix[] = {"output/stat/step"};
  char *dirname = fileio_create_dirname(prefix, step);
  if(NULL == dirname){
    return;
  }
  if(0 == myrank){
    fileio_mkdir(dirname);
  }
  // wait for the completion of mkdir
  MPI_Barrier(comm_cart);
  /* save parameters */
  if(0 == myrank){
    fileio_w_serial(dirname,  "num", 0, NULL, NPY_INT, sizeof(   int),  &num);
    fileio_w_serial(dirname, "step", 0, NULL, NPY_INT, sizeof(   int), &step);
    fileio_w_serial(dirname, "time", 0, NULL, NPY_DBL, sizeof(double), &time);
    config.output(dirname);
  }
  // save domain info (coordinates)
  if(0 == myrank){
    domain_save(dirname, domain);
  }
  // save collected statistics
  array_dump(comm_cart, dirname, "ux1", g_st_int->ux1);
  array_dump(comm_cart, dirname, "ux2", g_st_int->ux2);
  array_dump(comm_cart, dirname, "uy1", g_st_int->uy1);
  array_dump(comm_cart, dirname, "uy2", g_st_int->uy2);
  array_dump(comm_cart, dirname, "uz1", g_st_int->uz1);
  array_dump(comm_cart, dirname, "uz2", g_st_int->uz2);
  // don't forget to free memory holding directory name
  common_free(dirname);
}

