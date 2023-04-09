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

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include "sdecomp.h"
#include "../internal.h"
#include "internal.h"


/**
 * @brief initialise transpose plan
 * @param[in]  info            : struct contains information of process distribution
 * @param[in]  pencil_bef      : type of pencil before rotated
 * @param[in]  pencil_aft      : type of pencil after  rotated
 * @param[in]  glsizes         : global array size in each dimension
 * @param[in]  size_of_element : size of each element, e.g. sizeof(double)
 * @param[out] plan            : (success) a pointer to the created plan (struct)
 *                               (failure) undefined
 * @return                     : (success) 0
 *                               (failure) non-zero value
 */
int sdecomp_internal_transpose_construct(const sdecomp_info_t *info, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft, const size_t *glsizes, const size_t size_of_element, sdecomp_transpose_plan_t **plan){
  *plan = NULL;
  const char error_label[] = {"sdecomp.transpose.construct"};
  if(0 != sdecomp_internal_sanitise_null(error_label,    "info",    info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "glsizes", glsizes)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label,    "plan",    plan)) return 1;
  const size_t ndims = info->ndims;
  for(size_t dim = 0; dim < ndims; dim++){
    if(0 != sdecomp_internal_sanitise_glsize(error_label, glsizes[dim])) return 1;
  }
  if(0 != sdecomp_internal_sanitise_size_of_element(error_label, size_of_element)) return 1;
  if(2 == ndims){
    if(0 != sdecomp_internal_transpose_init_2d(error_label, info, pencil_bef, pencil_aft, glsizes, size_of_element, plan)) return 1;
  }else{
    if(0 != sdecomp_internal_transpose_init_3d(error_label, info, pencil_bef, pencil_aft, glsizes, size_of_element, plan)) return 1;
  }
  return 0;
}

/**
 * @brief execute transpose
 * @param[in]  plan    : transpose plan initialised by constructor
 * @param[in]  sendbuf : pointer to the input  buffer
 * @param[out] recvbuf : pointer to the output buffer
 * @return             : (success) 0
 *                       (failure) non-zero value
 */
int sdecomp_internal_transpose_execute(sdecomp_transpose_plan_t * restrict plan, const void * restrict sendbuf, void * restrict recvbuf){
  const char error_label[] = {"sdecomp.transpose.execute"};
  if(0 != sdecomp_internal_sanitise_null(error_label,    "plan",    plan)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "sendbuf", sendbuf)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "recvbuf", recvbuf)) return 1;
  const int * restrict scounts = plan->scounts;
  const int * restrict rcounts = plan->rcounts;
  const int * restrict sdispls = plan->sdispls;
  const int * restrict rdispls = plan->rdispls;
  const MPI_Datatype * restrict stypes = plan->stypes;
  const MPI_Datatype * restrict rtypes = plan->rtypes;
  MPI_Comm restrict comm_2d = plan->comm_2d;
  MPI_Alltoallw(
      sendbuf, scounts, sdispls, stypes,
      recvbuf, rcounts, rdispls, rtypes,
      comm_2d
  );
  return 0;
}

/**
 * @brief finalise transpose plan
 * @param[in,out] plan : transpose plan to be cleaned-up
 * @return             : (success) 0
 *                       (failure) non-zero value
 */
int sdecomp_internal_transpose_destruct(sdecomp_transpose_plan_t *plan){
  const char error_label[] = {"sdecomp.transpose.destruct"};
  if(0 != sdecomp_internal_sanitise_null(error_label, "plan", plan)) return 1;
  // type-free data types used by all-to-allw
  int nprocs_2d = 0;
  MPI_Comm_size(plan->comm_2d, &nprocs_2d);
  for(int n = 0; n < nprocs_2d; n++){
    MPI_Type_free(&plan->stypes[n]);
    MPI_Type_free(&plan->rtypes[n]);
  }
  // free communicator used by all-to-allw
  MPI_Comm_free(&(plan->comm_2d));
  // free struct itself and members
  sdecomp_internal_transpose_deallocate(plan);
  return 0;
}

