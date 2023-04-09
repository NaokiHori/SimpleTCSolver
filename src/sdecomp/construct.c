/*
 * Copyright 2022 Naoki Hori
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// https://github.com/NaokiHori/SimpleDecomp

// Construct sdecomp_info_t
// Allocate and initialise (decompose domain) sdecomp_info_t,
//   and return a pointer to it

#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <mpi.h>
#include "sdecomp.h"
#include "internal.h"


static int check_dims(const char error_label[], const MPI_Comm comm_default, const size_t ndims, const size_t *dims, bool *decomp_automatically){
  // if at least one of dims[dim] is non-zero,
  //   user specifies the decomposition pattern
  *decomp_automatically = true;
  for(size_t dim = 0; dim < ndims; dim++){
    if(0 != dims[dim]){
      *decomp_automatically = false;
      break;
    }
  }
  // dims are all-zero,
  //   indicating that I can use automatic domain decomposition
  //   (i.e. no need to check dims further)
  if(*decomp_automatically) return 0;
  // now I know user specifies how to decompose the domain,
  //   so I need to check whether the input is valid
  // check first element of "dims", which should be 1
  if(1 != dims[0]){
    // user tries to decompose in x direction,
    //   which is not allowed
    SDECOMP_ERROR(
        "dims contain at least one non-zero value, but dims[0] is not 1 (%zu is given)\n",
        error_label, dims[0]
    );
    return 1;
  }
  // check overflow
  // threshold: guaranteed to be larger than 65535,
  //   which is usually large enough (until 2040s or so...?)
  // NOTE: this is not the total number of processes,
  //   but number of processes in one direction
  for(size_t dim = 0; dim < ndims; dim++){
    const unsigned int threshold = UINT_MAX;
    if(threshold <= dims[dim]){
      SDECOMP_ERROR(
          "dims[%zu] (%zu: number of processes in one dimension) exceeds threshold (%u)\n",
          error_label, dim, dims[dim], threshold
      );
      return 1;
    }
  }
  // check number of total processes match
  int nprocs = 0;
  MPI_Comm_size(comm_default, &nprocs);
  size_t nprocs_user = 1;
  for(size_t dim = 0; dim < ndims; dim++){
    // compute product of dims, while checking overflow
    const size_t threshold = SIZE_MAX;
    if(threshold / nprocs_user <= dims[dim]){
      SDECOMP_ERROR(
          "product of dims exceeds threshold (%zu)\n",
          error_label, threshold
      );
      return 1;
    }
    nprocs_user *= dims[dim];
  }
  if((size_t)nprocs != nprocs_user){
    SDECOMP_ERROR(
        "Number of processes in the default communicator (%d) and the product of dims (%zu) do not match\n",
        error_label, nprocs, nprocs_user
    );
    return 1;
  }
  return 0;
}

static int create_new_communicator(const char error_label[], const MPI_Comm comm_default, const size_t ndims, const size_t *dims, const bool *periods, const bool decomp_automatically, MPI_Comm *comm_cart){
  // number of total processes participating in this decomposition
  int nprocs = 0;
  MPI_Comm_size(comm_default, &nprocs);
  // prepare parameters which will be passed to MPI APIs
  // NOTE: need to convert to int from size_t or bool
  int *dims_    = sdecomp_internal_calloc(error_label, ndims, sizeof(int));
  int *periods_ = sdecomp_internal_calloc(error_label, ndims, sizeof(int));
  if(NULL ==    dims_) return 1;
  if(NULL == periods_) return 1;
  if(decomp_automatically){
    for(size_t dim = 0; dim < ndims; dim++){
      // force 1st dimension NOT decomposed (assign 1)
      dims_[dim] = dim == 0 ? 1 : 0;
    }
    // let MPI library decompose domain
    MPI_Dims_create(nprocs, (int)ndims, dims_);
  }else{
    // copy user-specified value
    for(size_t dim = 0; dim < ndims; dim++){
      dims_[dim] = (int)(dims[dim]);
    }
  }
  for(size_t dim = 0; dim < ndims; dim++){
    // copy user-specified value
    periods_[dim] = (int)(periods[dim]);
  }
  // MPI rank in comm_default is no longer important,
  //   so I allow to override it
  const int reorder = 1;
  // create communicator
  *comm_cart = MPI_COMM_NULL;
  MPI_Cart_create(
      comm_default,
      (int)ndims,
      dims_,
      periods_,
      reorder,
      comm_cart
  );
  // clean-up integer buffers
  sdecomp_internal_free(   dims_);
  sdecomp_internal_free(periods_);
  return 0;
}

/**
 * @brief construct a structure sdecomp_info_t
 * @param[in] comm_default : MPI communicator which contains all processes
 *                             participating in the decomposition
 *                             (normally MPI_COMM_WORLD)
 * @param[in] ndims        : number of dimensions of the target domain
 * @param[in] dims         : number of processes in each dimension
 * @param[in] periods      : periodicities in each dimension
 * @param[out] info        : (success) a pointer to sdecomp_info_t
 *                           (failure) NULL pointer
 * @return                 : (success) 0
 *                           (failure) non-zero value
 */
int sdecomp_internal_construct(const MPI_Comm comm_default, const size_t ndims, const size_t *dims, const bool *periods, sdecomp_info_t **info){
  const char error_label[] = {"sdecomp.construct"};
  *info = NULL;
  if(0 != sdecomp_internal_sanitise_null(error_label,    "dims",    dims)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "periods", periods)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label,    "info",    info)) return 1;
  if(0 != sdecomp_internal_sanitise_ndims(error_label, ndims))             return 1;
  if(0 != sdecomp_internal_sanitise_comm(error_label, comm_default))       return 1;
  // check argument "dims"
  //   1. check whether user specifies how to decompose the domain
  //   2. if it is, whether the input is valid
  bool decomp_automatically = false;
  if(0 != check_dims(error_label, comm_default, ndims, dims, &decomp_automatically)) return 1;
  // create new communicator (x1 pencil)
  MPI_Comm comm_cart = MPI_COMM_NULL;
  if(0 != create_new_communicator(
        error_label,
        comm_default,
        ndims,
        dims,
        periods,
        decomp_automatically,
        &comm_cart
  )) return 1;
  // create sdecomp_info_t
  *info = sdecomp_internal_calloc(error_label, 1, sizeof(sdecomp_info_t));
  if(NULL == *info) return 1;
  // assign members
  (*info)->ndims = ndims;
  (*info)->comm_cart = comm_cart;
  return 0;
}

