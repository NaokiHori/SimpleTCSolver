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

#include <string.h>
#include <limits.h>
#include <mpi.h>
#include "sdecomp.h"
#include "../internal.h"
#include "internal.h"


#define SDECOMP_INTERNAL_NDIMS 2

static int convert_glsizes_to_sizes(const char error_label[], const sdecomp_pencil_t pencil, const size_t glsizes[SDECOMP_INTERNAL_NDIMS], size_t sizes[SDECOMP_INTERNAL_NDIMS]){
  // "glsizes" tell size of domain in each PHYSICAL direction (x, y)
  // convert them in terms of MEMORY order
  //                   <-- contiguous sparse -->
  // from x1 pencils : (glsizes[0],  glsizes[1])
  // from y1 pencils : (glsizes[1],  glsizes[0])
  if(SDECOMP_X1PENCIL == pencil){
    sizes[0] = glsizes[0];
    sizes[1] = glsizes[1];
    return 0;
  }else if(SDECOMP_Y1PENCIL == pencil){
    sizes[0] = glsizes[1];
    sizes[1] = glsizes[0];
    return 0;
  }else{
    sizes[0] = 0;
    sizes[1] = 0;
    SDECOMP_ERROR(
        "should not reach here (%s:%d)\n",
        error_label, __FILE__, __LINE__
    );
    return 1;
  }
}

/**
 * @brief initialise transpose plan for 2d
 * @param[in]  info            : struct contains information of process distribution
 * @param[in]  pencil_bef      : type of pencil before being rotated
 * @param[in]  pencil_aft      : type of pencil after  being rotated
 * @param[in]  glsizes         : global array size in each dimension
 * @param[in]  size_of_element : size of each element, e.g., sizeof(double)
 * @param[out] plan            : (success) a pointer to the created plan
 *                               (failure) undefined
 * @return                     : (success) 0
 *                               (failure) non-zero value
 */
int sdecomp_internal_transpose_init_2d(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil_bef, sdecomp_pencil_t pencil_aft, const size_t *glsizes, const size_t size_of_element, sdecomp_transpose_plan_t **plan){
  *plan = NULL;
  if(0 != sdecomp_internal_sanitise_pencil_pair_2d(error_label, pencil_bef, pencil_aft)) return 1;
  // create 2d communicator,
  //   in which comm_2d collective comm. will be called
  // NOTE: since now the domain is 2D, this is identical to
  //   the default Cartesian communicator comm_cart
  //   (= sdecomp.get_comm_cart(info)).
  // It is duplicated here to simplify the destruction process
  const MPI_Comm comm_cart = info->comm_cart;
  MPI_Comm comm_2d = MPI_COMM_NULL;
  MPI_Comm_dup(comm_cart, &comm_2d);
  // number of total process and my position
  int nprocs_2d = 0;
  int myrank_2d = 0;
  {
#define SDECOMP_INTERNAL_LOCAL_NDIMS 2 // 2d communicator
    int    dims[SDECOMP_INTERNAL_LOCAL_NDIMS] = {0};
    int periods[SDECOMP_INTERNAL_LOCAL_NDIMS] = {0};
    int  coords[SDECOMP_INTERNAL_LOCAL_NDIMS] = {0};
    MPI_Cart_get(comm_2d, SDECOMP_INTERNAL_LOCAL_NDIMS, dims, periods, coords);
    // I am interested in the memory-sparse direction
    // since comm_2d is defined in PHYSICAL (not memory) order,
    //   1: physical y direction is what I should extract here
    nprocs_2d =   dims[1];
    myrank_2d = coords[1];
#undef SDECOMP_INTERNAL_LOCAL_NDIMS
  }
  // get number of grid points in all directions
  //   in memory order (NOT physical order x, y)
  // NOTE: should never fail as long as
  //   glsizes is sanitised at the entrypoint
  size_t sizes[SDECOMP_INTERNAL_NDIMS] = {0};
  if(0 != convert_glsizes_to_sizes(error_label, pencil_bef, glsizes, sizes)) return 1;
  // check I can safely decompose the domain into chunks beforehand,
  //   i.e. local chunk sizes are positive (0 is not accepted)
  // I do it here for early return (before allocating buffers)
  for(int rank = 0; rank < nprocs_2d; rank++){
    size_t dummy = 0;
    if(
           0 != sdecomp_internal_kernel_get_mysize(error_label, sizes[0], nprocs_2d, rank, &dummy)
        || 0 != sdecomp_internal_kernel_get_mysize(error_label, sizes[1], nprocs_2d, rank, &dummy)
    ) return 1;
  }
  // allocate plan and its members
  if(0 != sdecomp_internal_transpose_allocate(error_label, nprocs_2d, plan)) return 1;
  (*plan)->comm_2d = comm_2d;
  // consider communication between my (myrank_2d-th) and your (yrrank_2d-th) pencils
  // base datatype having contiguous size_of_elemerror_label, ent bytes
  MPI_Datatype basetype = MPI_BYTE;
  MPI_Type_contiguous((int)size_of_element, basetype, &basetype);
  // no need to commit, but just in case
  MPI_Type_commit(&basetype);
  // NOTE: params for "get_mysize" and "get_offset" are already sanitised above
  //       and thus error check is skipped
  for(int yrrank_2d = 0; yrrank_2d < nprocs_2d; yrrank_2d++){
    // send
    {
      int *count         = &((*plan)->scounts[yrrank_2d]);
      int *displ         = &((*plan)->sdispls[yrrank_2d]);
      MPI_Datatype *type = &((*plan)->stypes [yrrank_2d]);
      size_t chunk_isize = 0;
      size_t chunk_ioffs = 0;
      size_t chunk_jsize = 0;
      sdecomp_internal_kernel_get_mysize(error_label, sizes[0], nprocs_2d, yrrank_2d, &chunk_isize);
      sdecomp_internal_kernel_get_offset(error_label, sizes[0], nprocs_2d, yrrank_2d, &chunk_ioffs);
      sdecomp_internal_kernel_get_mysize(error_label, sizes[1], nprocs_2d, myrank_2d, &chunk_jsize);
      // column data type 
      MPI_Type_create_hvector(
          (int)(chunk_jsize),
          1,
          (MPI_Aint)(size_of_element * sizes[0]),
          basetype,
          type
      );
      // repeat in x 
      MPI_Type_create_hvector(
          (int)(chunk_isize),
          1,
          (MPI_Aint)(size_of_element),
          *type,
          type
      );
      MPI_Type_commit(type);
      *count = 1;
      *displ = (int)(size_of_element * chunk_ioffs);
    }
    // recv
    {
      int *count         = &((*plan)->rcounts[yrrank_2d]);
      int *displ         = &((*plan)->rdispls[yrrank_2d]);
      MPI_Datatype *type = &((*plan)->rtypes [yrrank_2d]);
      size_t chunk_isize = 0;
      size_t chunk_jsize = 0;
      size_t chunk_joffs = 0;
      sdecomp_internal_kernel_get_mysize(error_label, sizes[0], nprocs_2d, myrank_2d, &chunk_isize);
      sdecomp_internal_kernel_get_mysize(error_label, sizes[1], nprocs_2d, yrrank_2d, &chunk_jsize);
      sdecomp_internal_kernel_get_offset(error_label, sizes[1], nprocs_2d, yrrank_2d, &chunk_joffs);
      // row data type 
      MPI_Type_create_hvector(
          (int)(chunk_isize),
          (int)(chunk_jsize),
          (MPI_Aint)(size_of_element * sizes[1]),
          basetype,
          type
      );
      MPI_Type_commit(type);
      *count = 1;
      *displ = (int)(size_of_element * chunk_joffs);
    }
  }
  // based plan is no longer used, so free here
  MPI_Type_free(&basetype);
  return 0;
}

#undef SDECOMP_INTERNAL_NDIMS
