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

#include "sdecomp.h"
#include "../internal.h"
#include "internal.h"

int sdecomp_internal_transpose_allocate(const char error_label[], const int nprocs_2d, sdecomp_transpose_plan_t **plan){
  *plan = sdecomp_internal_calloc(error_label, 1, sizeof(sdecomp_transpose_plan_t));
  int *scounts = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(int));
  int *rcounts = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(int));
  int *sdispls = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(int));
  int *rdispls = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(int));
  MPI_Datatype  *stypes = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(MPI_Datatype));
  MPI_Datatype  *rtypes = sdecomp_internal_calloc(error_label, (size_t)nprocs_2d, sizeof(MPI_Datatype));
  if(NULL ==   *plan) return 1;
  if(NULL == scounts) return 1;
  if(NULL == rcounts) return 1;
  if(NULL == sdispls) return 1;
  if(NULL == rdispls) return 1;
  if(NULL ==  stypes) return 1;
  if(NULL ==  rtypes) return 1;
  (*plan)->scounts = scounts;
  (*plan)->rcounts = rcounts;
  (*plan)->sdispls = sdispls;
  (*plan)->rdispls = rdispls;
  (*plan)->stypes  = stypes;
  (*plan)->rtypes  = rtypes;
  return 0;
}

int sdecomp_internal_transpose_deallocate(sdecomp_transpose_plan_t *plan){
  sdecomp_internal_free(plan->scounts);
  sdecomp_internal_free(plan->rcounts);
  sdecomp_internal_free(plan->sdispls);
  sdecomp_internal_free(plan->rdispls);
  sdecomp_internal_free(plan->stypes);
  sdecomp_internal_free(plan->rtypes);
  sdecomp_internal_free(plan);
  return 0;
}

