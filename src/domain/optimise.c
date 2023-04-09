#include <stdio.h>
#include <float.h>
#include <fftw3.h>
#include <mpi.h>
#include "common.h"
#include "sdecomp.h"
#include "domain.h"
#include "internal.h"

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
 * @brief return optimised MPI domain decomposition information
 * @return : structure sdecomp_t containing domain info
 */
sdecomp_info_t *optimise_sdecomp_init(const bool uniformx, const int *glsizes){
  // periodicity in each dimension
  const int periods[NDIMS] = {0, 1, 1};
  // number of processes in each dimension,
  //   which is to be optimised
  int dims_optimum[NDIMS] = {0, 0, 0};
  double wtime_optimum = DBL_MAX;
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  // factorise and decide dims
  for(int n = 1; n <= nprocs; n++){
    if(nprocs % n == 0){
      const int dims[NDIMS] = {1, n, nprocs / n};
      // sanitise, for all dimension, dims should not exceed the number of grid points
      bool valid = true;
      {
        const int glsizes_[NDIMS] = {
          glsizes[0]        ,
          glsizes[1] / 2 + 1, // after FFTed
          glsizes[2]        ,
        };
        for(int dim1 = 0; dim1 < NDIMS; dim1++){
          const int glsize = glsizes_[dim1];
          for(int dim0 = 0; dim0 < NDIMS; dim0++){
            const int np = dims[dim0];
            if(np > glsize){
              valid = false;
            }
          }
        }
      }
      if(!valid) continue;
      // execute transposes which are used to solve Poisson equation
      //   and check how long they take in total
      // for the time being only transposes when solving Poisson equation are considered,
      //   i.e. implicity and implicitz also request transposes, which are neglected for now
      sdecomp_info_t *info = NULL;
      {
        const size_t dims_[NDIMS] = {
          (size_t)dims[0],
          (size_t)dims[1],
          (size_t)dims[2],
        };
        const bool periods_[NDIMS] = {
          (bool)periods[0],
          (bool)periods[1],
          (bool)periods[2],
        };
        sdecomp.construct(MPI_COMM_WORLD, NDIMS, dims_, periods_, &info);
      }
      // initialise pencils and rotations
      const size_t r_dsize = sizeof(      double);
      const size_t c_dsize = sizeof(fftw_complex);
      const size_t r_sizes[NDIMS] = {
        (size_t)glsizes[0],
        (size_t)glsizes[1],
        (size_t)glsizes[2],
      };
      const size_t c_sizes[NDIMS] = {
        (size_t)glsizes[0]        ,
        (size_t)glsizes[1] / 2 + 1, // after FFTed
        (size_t)glsizes[2]        ,
      };
      double       *r_x1pcnl = NULL;
      double       *r_y1pcnl = NULL;
      fftw_complex *c_y1pcnl = NULL;
      fftw_complex *c_z1pcnl = NULL;
      fftw_complex *c_x2pcnl = NULL;
      sdecomp_transpose_plan_t *r_x1_to_y1 = NULL;
      sdecomp_transpose_plan_t *r_y1_to_x1 = NULL;
      sdecomp_transpose_plan_t *c_y1_to_z1 = NULL;
      sdecomp_transpose_plan_t *c_z1_to_y1 = NULL;
      sdecomp_transpose_plan_t *c_z1_to_x2 = NULL;
      sdecomp_transpose_plan_t *c_x2_to_z1 = NULL;
      sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_sizes, r_dsize, &r_x1_to_y1);
      sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_sizes, r_dsize, &r_y1_to_x1);
      sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, c_sizes, c_dsize, &c_y1_to_z1);
      sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, c_sizes, c_dsize, &c_z1_to_y1);
      sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_X2PENCIL, c_sizes, c_dsize, &c_z1_to_x2);
      sdecomp.transpose.construct(info, SDECOMP_X2PENCIL, SDECOMP_Z1PENCIL, c_sizes, c_dsize, &c_x2_to_z1);
      size_t r_x1sizes[NDIMS] = {0};
      size_t r_y1sizes[NDIMS] = {0};
      size_t c_y1sizes[NDIMS] = {0};
      size_t c_z1sizes[NDIMS] = {0};
      size_t c_x2sizes[NDIMS] = {0};
      for(int dim = 0; dim < NDIMS; dim++){
        sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_sizes[dim], r_x1sizes + dim);
        sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_sizes[dim], r_y1sizes + dim);
        sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_sizes[dim], c_y1sizes + dim);
        sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, c_sizes[dim], c_z1sizes + dim);
        sdecomp.get_pencil_mysize(info, SDECOMP_X2PENCIL, dim, c_sizes[dim], c_x2sizes + dim);
      }
      r_x1pcnl = common_calloc(r_x1sizes[0] * r_x1sizes[1] * r_x1sizes[2], r_dsize);
      r_y1pcnl = common_calloc(r_y1sizes[0] * r_y1sizes[1] * r_y1sizes[2], r_dsize);
      c_y1pcnl = common_calloc(c_y1sizes[0] * c_y1sizes[1] * c_y1sizes[2], c_dsize);
      c_z1pcnl = common_calloc(c_z1sizes[0] * c_z1sizes[1] * c_z1sizes[2], c_dsize);
      c_x2pcnl = common_calloc(c_x2sizes[0] * c_x2sizes[1] * c_x2sizes[2], c_dsize);
      // execute transpose, repeat for "niter" times
      const int niter = 4;
      double wtimes[2] = {0., 0.};
      wtimes[0] = common_get_wtime();
      for(int iter = 0; iter < niter; iter++){
        sdecomp.transpose.execute(r_x1_to_y1, r_x1pcnl, r_y1pcnl);
        sdecomp.transpose.execute(c_y1_to_z1, c_y1pcnl, c_z1pcnl);
        if(!uniformx){
          sdecomp.transpose.execute(c_z1_to_x2, c_z1pcnl, c_x2pcnl);
          sdecomp.transpose.execute(c_x2_to_z1, c_x2pcnl, c_z1pcnl);
        }
        sdecomp.transpose.execute(c_z1_to_y1, c_z1pcnl, c_y1pcnl);
        sdecomp.transpose.execute(r_y1_to_x1, r_y1pcnl, r_x1pcnl);
      }
      wtimes[1] = common_get_wtime();
      // clean-up tentative transpose plans and buffers
      sdecomp.transpose.destruct(r_x1_to_y1);
      sdecomp.transpose.destruct(c_y1_to_z1);
      sdecomp.transpose.destruct(c_z1_to_x2);
      sdecomp.transpose.destruct(c_x2_to_z1);
      sdecomp.transpose.destruct(c_z1_to_y1);
      sdecomp.transpose.destruct(r_y1_to_x1);
      common_free(r_x1pcnl);
      common_free(r_y1pcnl);
      common_free(c_y1pcnl);
      common_free(c_z1pcnl);
      common_free(c_x2pcnl);
      // clean-up current sdecomp config
      sdecomp.destruct(info);
      // check time
      double wtime = (wtimes[1] - wtimes[0]) / niter;
      if(wtime < wtime_optimum){
        dims_optimum[0] = dims[0];
        dims_optimum[1] = dims[1];
        dims_optimum[2] = dims[2];
        wtime_optimum = wtime;
      }
      int myrank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      if(myrank == 0){
        const int nd = get_ndigits(nprocs);
        printf("dims: [%*d, %*d, %*d]: % .7e [sec]\n", nd, dims[0], nd, dims[1], nd, dims[2], wtime);
      }
    }
  }
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    const int nd = get_ndigits(nprocs);
    printf("Conclusive domain decomposition: [%*d, %*d, %*d]\n", nd, dims_optimum[0], nd, dims_optimum[1], nd, dims_optimum[2]);
  }
  // create sdecomp which will be used in the main run
  {
    // in size_t
    const size_t dims_optimum_[NDIMS] = {
      (size_t)dims_optimum[0],
      (size_t)dims_optimum[1],
      (size_t)dims_optimum[2],
    };
    // in bool
    const bool periods_[NDIMS] = {
      (bool)periods[0],
      (bool)periods[1],
      (bool)periods[2],
    };
    sdecomp_info_t *info = NULL;
    sdecomp.construct(MPI_COMM_WORLD, NDIMS, dims_optimum_, periods_, &info);
    return info;
  }
}
