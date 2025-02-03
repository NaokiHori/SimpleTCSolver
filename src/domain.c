#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <fftw3.h>
#include "timer.h"
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "fileio.h"
#include "array_macros/domain/xf.h"
#include "array_macros/domain/xc.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"

/**
 * @brief load members in domain_t
 * @param[in]  dirname : name of directory from which data is loaded
 * @param[out] domain  : global domain sizes and resolutions
 * @return             : error code
 */
static int domain_load(
    const char dirname[],
    domain_t * domain
){
  bool * is_curved = &domain->is_curved;
  size_t * glsizes = domain->glsizes;
  double * restrict lengths = domain->lengths;
  double * restrict * restrict xf = &domain->xf;
  double * restrict * restrict xc = &domain->xc;
  if(0 != fileio.r_serial(dirname, "is_curved", 0, NULL, fileio.npy_boolean, sizeof(bool), is_curved)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths)){
    return 1;
  }
  *xf = memory_calloc(glsizes[0] + 1, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), *xf)){
    return 1;
  }
  *xc = memory_calloc(glsizes[0] + 2, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), *xc)){
    return 1;
  }
  return 0;
}

/**
 * @brief save members in domain_t
 * @param[in] dirname : name of directory to which data is saved
 * @param[in] domain  : global domain sizes and resolutions
 * @return            : error code
 */
int domain_save(
    const char dirname[],
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // since this is a serial operation,
  //   other processes are not involved
  if(root != myrank){
    return 0;
  }
  const bool is_curved = domain->is_curved;
  const size_t * glsizes = domain->glsizes;
  const double * lengths = domain->lengths;
  const double * xf      = domain->xf;
  const double * xc      = domain->xc;
  fileio.w_serial(dirname, "is_curved", 0, NULL, fileio.npy_boolean, sizeof(bool), &is_curved);
  fileio.w_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes);
  fileio.w_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths);
  fileio.w_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), xf);
  fileio.w_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), xc);
  return 0;
}

static double * allocate_and_init_hxxf(
    const int isize,
    const double * xc
){
  assert(XC_NADDS[0] == 1);
  assert(XC_NADDS[1] == 1);
  assert(HXXF_NADDS[0] == 0);
  assert(HXXF_NADDS[1] == 1);
  double * hxxf = memory_calloc(isize + HXXF_NADDS[0] + HXXF_NADDS[1], sizeof(double));
  // x scale factors at x cell faces
  for(int i = 1; i <= isize + 1; i++){
    HXXF(i  ) = XC(i  ) - XC(i-1);
  }
  return hxxf;
}

static double * allocate_and_init_hxxc(
    const int isize,
    const double * xf
){
  assert(XF_NADDS[0] == 0);
  assert(XF_NADDS[1] == 1);
  assert(HXXC_NADDS[0] == 1);
  assert(HXXC_NADDS[1] == 1);
  double * hxxc = memory_calloc(isize + HXXC_NADDS[0] + HXXC_NADDS[1], sizeof(double));
  // x scale factors at x cell centers
  HXXC(        0) = 0.5 * XF(        2) - 0.5 * XF(    1);
  for(int i = 1; i <= isize; i++){
    HXXC(i  ) = XF(i+1) - XF(i  );
  }
  HXXC(isize + 1) = 0.5 * XF(isize + 1) - 0.5 * XF(isize);
  return hxxc;
}

static double * allocate_and_init_hyxf(
    const bool is_curved,
    const int isize,
    const double * xf,
    const double ly,
    const size_t ny
){
  assert(HYXF_NADDS[0] == 0);
  assert(HYXF_NADDS[1] == 1);
  double * hyxf = memory_calloc(isize + HYXF_NADDS[0] + HYXF_NADDS[1], sizeof(double));
  // y scale factors at x cell faces
  for(int i = 1; i <= isize + 1; i++){
    HYXF(i  ) = is_curved ? XF(i  ) * ly / ny : ly / ny;
  }
  return hyxf;
}

static double * allocate_and_init_hyxc(
    const bool is_curved,
    const int isize,
    const double * xc,
    const double ly,
    const size_t ny
){
  assert(HYXC_NADDS[0] == 1);
  assert(HYXC_NADDS[1] == 1);
  double * hyxc = memory_calloc(isize + HYXC_NADDS[0] + HYXC_NADDS[1], sizeof(double));
  // y scale factors at x cell centers
  for(int i = 0; i <= isize + 1; i++){
    HYXC(i  ) = is_curved ? XC(i  ) * ly / ny : ly / ny;
  }
  return hyxc;
}

static double * allocate_and_init_jdxf(
    const int isize,
    const double * hxxf,
    const double * hyxf
){
  assert(JDXF_NADDS[0] == 0);
  assert(JDXF_NADDS[1] == 1);
  double * jdxf = memory_calloc(isize + JDXF_NADDS[0] + JDXF_NADDS[1], sizeof(double));
  // jacobian determinants at x cell faces
  for(int i = 1; i <= isize + 1; i++){
    JDXF(i  ) = HXXF(i  ) * HYXF(i  );
  }
  return jdxf;
}

static double * allocate_and_init_jdxc(
    const int isize,
    const double * hxxc,
    const double * hyxc
){
  assert(JDXC_NADDS[0] == 1);
  assert(JDXC_NADDS[1] == 1);
  double * jdxc = memory_calloc(isize + JDXC_NADDS[0] + JDXC_NADDS[1], sizeof(double));
  // jacobian determinants at x cell centers
  for(int i = 0; i <= isize + 1; i++){
    JDXC(i  ) = HXXC(i  ) * HYXC(i  );
  }
  return jdxc;
}

static void report(
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("DOMAIN\n");
    if (domain->is_curved) {
      puts("\tCylindrical domain");
    } else {
      puts("\tCartesian domain");
    }
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tglsizes[%u]: %zu\n", dim, domain->glsizes[dim]);
    }
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tlengths[%u]: % .7e\n", dim, domain->lengths[dim]);
    }
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial conditions are stored
 * @param[out] domain     : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int domain_init(
    const char dirname_ic[],
    domain_t * domain
){
  sdecomp_info_t ** info    = &domain->info;
  bool * is_curved = &domain->is_curved;
  size_t * restrict glsizes = domain->glsizes;
  size_t * restrict mysizes = domain->mysizes;
  size_t * restrict offsets = domain->offsets;
  double * restrict lengths = domain->lengths;
  double * restrict * xf    = &domain->xf;
  double * restrict * xc    = &domain->xc;
  double * restrict * jdxf  = &domain->jdxf;
  double * restrict * jdxc  = &domain->jdxc;
  double * restrict * hxxf  = &domain->hxxf;
  double * restrict * hxxc  = &domain->hxxc;
  double * restrict * hyxf  = &domain->hyxf;
  double * restrict * hyxc  = &domain->hyxc;
  // load spatial information
  if(0 != domain_load(dirname_ic, domain)){
    return 1;
  }
  // scale factors
  *hxxf = allocate_and_init_hxxf(glsizes[0], *xc);
  *hxxc = allocate_and_init_hxxc(glsizes[0], *xf);
  *hyxf = allocate_and_init_hyxf(*is_curved, glsizes[0], *xf, lengths[1], glsizes[1]);
  *hyxc = allocate_and_init_hyxc(*is_curved, glsizes[0], *xc, lengths[1], glsizes[1]);
  // Jacobian determinants at x cell faces and centers
  *jdxf = allocate_and_init_jdxf(glsizes[0], *hxxf, *hyxf);
  *jdxc = allocate_and_init_jdxc(glsizes[0], *hxxc, *hyxc);
  if(0 != sdecomp.construct(
        MPI_COMM_WORLD,
        NDIMS,
        (size_t [NDIMS]){0, 0},
        (bool [NDIMS]){false, true},
        info
  )) return 1;
  // local array sizes and offsets
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], mysizes + dim);
    sdecomp.get_pencil_offset(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], offsets + dim);
  }
  report(domain);
  return 0;
}

