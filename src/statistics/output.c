#include <string.h>
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

static const char dirname_prefix[] = {"output/stat/step"};
static const int dirname_ndigits = 10;
static const bool reduction = true;

/**
 * @brief reduce N-dimensional (xy, xyz) array to 1D (x) vector
 * @param[in] info     : information related to MPI domain decomposition
 * @param[in] dirname  : name of directory
 * @param[in] dsetname : name of dataset
 * @param[in] array    : N-dimensional array to be reduced
 */
static int reduce_and_write(const sdecomp_info_t * restrict info, const char dirname[restrict], const char dsetname[restrict], const array_t * restrict array){
  const int isize = array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1];
  const int jsize = array->mysizes[1];
  const int ksize = array->mysizes[2];
  const double * restrict data = array->data;
  double * restrict vec = common_calloc(isize, sizeof(double));
  for(int k = 0; k < ksize; k++){
    for(int j = 0; j < jsize; j++){
      for(int i = 0; i < isize; i++){
        vec[i] += data[k * jsize * isize + j * isize + i];
      }
    }
  }
  int myrank = 0;
  sdecomp.get_comm_rank(info, &myrank);
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, vec, isize, MPI_DOUBLE, MPI_SUM, comm_cart);
  if(0 == myrank){
    fileio_w_serial(dirname, dsetname, 1, (size_t [1]){isize}, NPY_DBL, sizeof(double) * isize, vec);
  }
  common_free(vec);
  return 0;
}

/**
 * @brief save structures which contains collected statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @param[in] time   : current time units
 */
void output(const domain_t *domain, const int step, const double time){
  // when no statistics are collected (num is 0),
  //   no reason to save, so abort
  const int num = g_st_int->num;
  if(0 == num){
    return;
  }
  // set directory name
  // allocate directory name
  const size_t dirname_nchars = strlen(dirname_prefix) + dirname_ndigits;
  char *dirname = common_calloc(dirname_nchars + 2, sizeof(char));
  snprintf(dirname, dirname_nchars + 1, "%s%0*d", dirname_prefix, dirname_ndigits, step);
  // get communicator to identify the main process
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory and save scalars from main process
  if(0 == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio_mkdir(dirname);
    // save scalars
    fileio_w_serial(dirname,  "num",        0, NULL, NPY_INT, sizeof(   int),  &num);
    fileio_w_serial(dirname, "step",        0, NULL, NPY_INT, sizeof(   int), &step);
    fileio_w_serial(dirname, "time",        0, NULL, NPY_DBL, sizeof(double), &time);
    fileio_w_serial(dirname, "diffusivity", 0, NULL, NPY_DBL, sizeof(double), &g_st_int->diffusivity);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // save domain info (coordinates)
  domain_save(dirname, domain);
  // save collected statistics
  if(reduction){
    reduce_and_write(domain->info, dirname, "ux1", g_st_int->ux1);
    reduce_and_write(domain->info, dirname, "ux2", g_st_int->ux2);
    reduce_and_write(domain->info, dirname, "uy1", g_st_int->uy1);
    reduce_and_write(domain->info, dirname, "uy2", g_st_int->uy2);
    reduce_and_write(domain->info, dirname, "uz1", g_st_int->uz1);
    reduce_and_write(domain->info, dirname, "uz2", g_st_int->uz2);
  }else{
    MPI_Comm comm_cart = MPI_COMM_NULL;
    sdecomp.get_comm_cart(domain->info, &comm_cart);
    array_dump(comm_cart, dirname, "ux1", g_st_int->ux1);
    array_dump(comm_cart, dirname, "ux2", g_st_int->ux2);
    array_dump(comm_cart, dirname, "uy1", g_st_int->uy1);
    array_dump(comm_cart, dirname, "uy2", g_st_int->uy2);
    array_dump(comm_cart, dirname, "uz1", g_st_int->uz1);
    array_dump(comm_cart, dirname, "uz2", g_st_int->uz2);
  }
  common_free(dirname);
}

