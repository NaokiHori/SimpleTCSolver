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

#include <stdint.h>
#include "internal.h"

static int compare_glsize_and_nprocs(const char error_label[], const size_t glsize, const int nprocs){
  // I try to split glsize elements with nprocs nprocesses
  // then it does not make sense if glsize is less than nprocs
  if((int)glsize < nprocs){
    SDECOMP_ERROR(
        "glsize (%zu) should be equal to or more than nprocs (%d)\n",
        error_label, glsize, nprocs
    );
    return 1;
  }
  // check overflow
  if(SIZE_MAX - glsize < (size_t)nprocs){
    SDECOMP_ERROR(
        "sum of glsize (%zu) and nprocs (%d) exceeds SIZE_MAX (%zu)\n",
        error_label, glsize, nprocs, SIZE_MAX
    );
    return 1;
  }
  return 0;
}

static int compare_nprocs_and_myrank(const char error_label[], const int nprocs, const int myrank){
  // myrank is an ID in the communicator which holds this process
  // then it is strange if myrank is equal to or more than nprocs
  if(myrank >= nprocs){
    SDECOMP_ERROR(
        "myrank (%d) should be smaller than nprocs (%d)\n",
        error_label, myrank, nprocs
    );
    return 1;
  }
  return 0;
}

/* compute number of grid points taken care of by the process */
int sdecomp_internal_kernel_get_mysize(const char error_label[], const size_t glsize, const int nprocs, const int myrank, size_t *mysize){
  if(0 != sdecomp_internal_sanitise_null(error_label, "mysize", mysize)) return 1;
  if(0 != sdecomp_internal_sanitise_nprocs(error_label, nprocs)) return 1;
  if(0 != sdecomp_internal_sanitise_myrank(error_label, myrank)) return 1;
  if(0 != sdecomp_internal_sanitise_glsize(error_label, glsize)) return 1;
  if(0 != compare_glsize_and_nprocs(error_label, glsize, nprocs)) return 1;
  if(0 != compare_nprocs_and_myrank(error_label, nprocs, myrank)) return 1;
  // example: glsize: 10, nprocs: 3 (3 processes in total)
  //   myrank = 0 -> mysize = 3
  //   myrank = 1 -> mysize = 3
  //   myrank = 2 -> mysize = 4
  // NOTE: sum of "mysize"s is "glsize"
  *mysize = (glsize + (size_t)myrank) / (size_t)nprocs;
  return 0;
}

/* compute offset of the grid point taken care of by the process */
int sdecomp_internal_kernel_get_offset(const char error_label[], const size_t glsize, const int nprocs, const int myrank, size_t *offset){
  if(0 != sdecomp_internal_sanitise_null(error_label, "offset", offset)) return 1;
  if(0 != sdecomp_internal_sanitise_nprocs(error_label, nprocs)) return 1;
  if(0 != sdecomp_internal_sanitise_myrank(error_label, myrank)) return 1;
  if(0 != sdecomp_internal_sanitise_glsize(error_label, glsize)) return 1;
  if(0 != compare_glsize_and_nprocs(error_label, glsize, nprocs)) return 1;
  if(0 != compare_nprocs_and_myrank(error_label, nprocs, myrank)) return 1;
  // example: glsize: 10, nprocs: 3 (3 processes in total)
  //   myrank = 0 -> offset = 0
  //   myrank = 1 -> offset = 3
  //   myrank = 2 -> offset = 6
  *offset = 0;
  for(int i = 0; i < myrank; i++){
    size_t val = 0;
    if(0 != sdecomp_internal_kernel_get_mysize(error_label, glsize, nprocs, i, &val)) return 1;
    *offset += val;
  }
  return 0;
}

