#include <stdbool.h>
#include "common.h"
#include "fileio.h"
#include "domain.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "internal.h"

/**
 * @brief load members in domain_t
 * @param[in]  dirname : name of directory from which data is loaded
 * @param[out] domain  : global domain sizes and resolutions
 * @return             : error code
 */
int domain_load(const char dirname[], domain_t *domain){
  bool   * restrict uniformx = &(domain->uniformx);
  int    * restrict *glsizes = &(domain->glsizes);
  double * restrict *lengths = &(domain->lengths);
  double * restrict *xf      = &(domain->xf);
  double * restrict *xc      = &(domain->xc);
  if(0 != fileio_r_serial(dirname, "uniformx", 0, NULL, NPY_BOL, sizeof(bool), uniformx)){
    return 1;
  }
  *glsizes = common_calloc(NDIMS, sizeof(   int));
  *lengths = common_calloc(NDIMS, sizeof(double));
  if(0 != fileio_r_serial(dirname, "glisize", 0, NULL, NPY_INT, sizeof(int), &(*glsizes)[0])){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "gljsize", 0, NULL, NPY_INT, sizeof(int), &(*glsizes)[1])){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "glksize", 0, NULL, NPY_INT, sizeof(int), &(*glsizes)[2])){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "lx", 0, NULL, NPY_DBL, sizeof(double), &(*lengths)[0])){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "ly", 0, NULL, NPY_DBL, sizeof(double), &(*lengths)[1])){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "lz", 0, NULL, NPY_DBL, sizeof(double), &(*lengths)[2])){
    return 1;
  }
  // xf
  {
    const size_t shape[1] = {(size_t)XF_NITEMS_0((*glsizes)[0])};
    *xf = common_calloc(shape[0], sizeof(double));
    const size_t size = sizeof(double) * shape[0];
    if(0 != fileio_r_serial(dirname, "xf", 1, shape, NPY_DBL, size, *xf)){
      return 1;
    }
  }
  // xc
  {
    const size_t shape[1] = {(size_t)XC_NITEMS_0((*glsizes)[0])};
    *xc = common_calloc(shape[0], sizeof(double));
    const size_t size = sizeof(double) * shape[0];
    if(0 != fileio_r_serial(dirname, "xc", 1, shape, NPY_DBL, size, *xc)){
      return 1;
    }
  }
  return 0;
}

/**
 * @brief save members in domain_t
 * @param[in] dirname : name of directory to which data is saved
 * @param[in] domain  : global domain sizes and resolutions
 * @return            : error code
 */
int domain_save(const char dirname[], const domain_t *domain){
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // serial operation
  if(0 != myrank) return 0;
  const bool   uniformx = domain->uniformx;
  const int    *glsizes = domain->glsizes;
  const double *lengths = domain->lengths;
  const double *xf      = domain->xf;
  const double *xc      = domain->xc;
  fileio_w_serial(dirname, "uniformx", 0, NULL, NPY_BOL, sizeof(bool), &uniformx);
  fileio_w_serial(dirname, "glisize", 0, NULL, NPY_INT, sizeof(int), &(glsizes[0]));
  fileio_w_serial(dirname, "gljsize", 0, NULL, NPY_INT, sizeof(int), &(glsizes[1]));
  fileio_w_serial(dirname, "glksize", 0, NULL, NPY_INT, sizeof(int), &(glsizes[2]));
  fileio_w_serial(dirname, "lx", 0, NULL, NPY_DBL, sizeof(double), &(lengths[0]));
  fileio_w_serial(dirname, "ly", 0, NULL, NPY_DBL, sizeof(double), &(lengths[1]));
  fileio_w_serial(dirname, "lz", 0, NULL, NPY_DBL, sizeof(double), &(lengths[2]));
  // xf
  {
    const size_t shape[1] = {(size_t)XF_NITEMS_0(glsizes[0])};
    const size_t size = sizeof(double) * shape[0];
    fileio_w_serial(dirname, "xf", 1, shape, NPY_DBL, size, xf);
  }
  // xc
  {
    const size_t shape[1] = {(size_t)XC_NITEMS_0(glsizes[0])};
    const size_t size = sizeof(double) * shape[0];
    fileio_w_serial(dirname, "xc", 1, shape, NPY_DBL, size, xc);
  }
  return 0;
}

