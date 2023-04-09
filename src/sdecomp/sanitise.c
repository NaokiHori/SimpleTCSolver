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

// utility functions to check arguments
//   return 0 if the argument is   valid
//   return 1 if the argument is invalid and dump an error message

#include <stdbool.h>
#include <limits.h>
#include "sdecomp.h"
#include "internal.h"


// NULL checker
int sdecomp_internal_sanitise_null(const char error_label[], const char ptr_name[], const void *ptr){
  if(NULL != ptr) return 0;
  SDECOMP_ERROR(
      "NULL pointer is passed (%s)\n",
      error_label, ptr_name
  );
  return 1;
}

// check ndims, which should be 2 or 3
int sdecomp_internal_sanitise_ndims(const char error_label[], const size_t ndims){
  if(2 == ndims || 3 == ndims) return 0;
  SDECOMP_ERROR(
      "ndims (%zu is given) should be one of 2 or 3\n",
      error_label, ndims
  );
  return 1;
}

// check comm, which should NOT be MPI_COMM_NULL
int sdecomp_internal_sanitise_comm(const char error_label[], const MPI_Comm comm){
  if(MPI_COMM_NULL != comm) return 0;
  SDECOMP_ERROR(
      "MPI_COMM_NULL is given\n",
      error_label
  );
  return 1;
}

// check the given pencil is valid
//   2D: SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL
//   3D: SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL
//       SDECOMP_X2PENCIL, SDECOMP_Y2PENCIL, SDECOMP_Z2PENCIL
int sdecomp_internal_sanitise_pencil(const char error_label[], const size_t ndims, const sdecomp_pencil_t pencil){
  bool is_valid = false;
  if(2 == ndims){
    is_valid
      =  SDECOMP_X1PENCIL == pencil
      || SDECOMP_Y1PENCIL == pencil;
  }else{
    is_valid
      =  SDECOMP_X1PENCIL == pencil
      || SDECOMP_Y1PENCIL == pencil
      || SDECOMP_Z1PENCIL == pencil
      || SDECOMP_X2PENCIL == pencil
      || SDECOMP_Y2PENCIL == pencil
      || SDECOMP_Z2PENCIL == pencil;
  }
  if(is_valid) return 0;
  SDECOMP_ERROR(
      "invalid pencil: %u\n",
      error_label, pencil
  );
  return 1;
}

// check the given direction is valid
//   2D: SDECOMP_XDIR, SDECOMP_YDIR
//   3D: SDECOMP_XDIR, SDECOMP_YDIR, SDECOMP_ZDIR
int sdecomp_internal_sanitise_dir(const char error_label[], const size_t ndims, const sdecomp_dir_t dir){
  bool is_valid = false;
  if(2 == ndims){
    is_valid
      =  SDECOMP_XDIR == dir
      || SDECOMP_YDIR == dir;
  }else{
    is_valid
      =  SDECOMP_XDIR == dir
      || SDECOMP_YDIR == dir
      || SDECOMP_ZDIR == dir;
  }
  if(is_valid) return 0;
  SDECOMP_ERROR(
      "invalid dir: %u\n",
      error_label, dir
  );
  return 1;
}

// check number of process, which should be positive
// I do this since MPI takes integer
int sdecomp_internal_sanitise_nprocs(const char error_label[], const int nprocs){
  if(nprocs > 0) return 0;
  SDECOMP_ERROR(
      "total number of processes is non-positive (%d)\n",
      error_label, nprocs
  );
  return 1;
}

// check my MPI rank, which should be non-negative
// I do this since MPI takes integer
int sdecomp_internal_sanitise_myrank(const char error_label[], const int myrank){
  if(myrank >= 0) return 0;
  SDECOMP_ERROR(
      "MPI rank is negative (%d)\n",
      error_label, myrank
  );
  return 1;
}

// check number of global grid points,
//   which should be smaller than INT_MAX (for size_t -> int conversion purpose)
int sdecomp_internal_sanitise_glsize(const char error_label[], const size_t glsize){
  const int maxglsize = INT_MAX;
  if((size_t)maxglsize > glsize) return 0;
  SDECOMP_ERROR(
      "glsize %zu exceeds %d\n",
      error_label, glsize, maxglsize
  );
  return 1;
}

// reject too big array element
int sdecomp_internal_sanitise_size_of_element(const char error_label[], const size_t size_of_element){
  const size_t threshold = USHRT_MAX;
  if(threshold > size_of_element) return 0;
  SDECOMP_ERROR(
      "size_of_element = %zu should be smaller than %zu\n",
      error_label, size_of_element, threshold
  );
  return 1;
}

// check the given pencil pair is valid (2D)
int sdecomp_internal_sanitise_pencil_pair_2d(const char error_label[], const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  const size_t ndims = 2;
  if(0 != sdecomp_internal_sanitise_pencil(error_label, ndims, pencil_bef)) return 1;
  if(0 != sdecomp_internal_sanitise_pencil(error_label, ndims, pencil_aft)) return 1;
  const sdecomp_pencil_t total = 2;
  const sdecomp_pencil_t value = (total + pencil_aft - pencil_bef) % total;
  // valid pencil pairs
  if(1 == value) return 0;
  // invalid pencil pair
  SDECOMP_ERROR(
      "pair of pencil_bef (%u) pencil_aft (%u) is not valid\n",
      error_label, pencil_bef, pencil_aft
  );
  return 1;
}

// check the given pencil pair is valid (3D)
int sdecomp_internal_sanitise_pencil_pair_3d(const char error_label[], const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft, bool *is_forward){
  const size_t ndims = 3;
  if(0 != sdecomp_internal_sanitise_pencil(error_label, ndims, pencil_bef)) return 1;
  if(0 != sdecomp_internal_sanitise_pencil(error_label, ndims, pencil_aft)) return 1;
  const sdecomp_pencil_t total = 6;
  const sdecomp_pencil_t value = (total + pencil_aft - pencil_bef) % total;
  // valid pencil pairs, clockwise (forward)
  if(1 == value){
    *is_forward = true;
    return 0;
  }else if(total - 1 == value){
    *is_forward = false;
    return 0;
  }
  // invalid pencil pair
  SDECOMP_ERROR(
      "pair of pencil_bef (%u) pencil_aft (%u) is not valid\n",
      error_label, pencil_bef, pencil_aft
  );
  return 1;
}

